use sha2::{Digest, Sha256};
use std::env;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::time::UNIX_EPOCH;

// GENOMEIDX1: the only path you must set — absolute path to your GRCh38 reference FASTA.
// All other tools are bundled under software/ and scripts/ in the repository root.
const GENOMEIDX1: &str = "";            // e.g. "/data/ref/Homo_sapiens.GRCh38.dna.toplevel.fna"
const SOFTWAREDIR: &str = "software";   // relative to repo root; contains VarScan.v2.3.9.jar
const SCRIPTSDIR: &str = "scripts";     // relative to repo root; contains fpfilter.pl

const BAM_SUFFIX: &str = "_final.bam";
const FILE_PAIRS_LIST: &str = "sample_pairs.csv";
const PAIRS_SUFFIX: &str = "_final.bam";

const TARGET_BED: &str = "";
const VCF_SAMPLE_LIST: &str = "";

const MIN_COVERAGE: i32 = 20;
const MIN_COVERAGE_NORMAL: i32 = 10;
const MIN_COVERAGE_TUMOR: i32 = 20;
const MIN_BASE_QUAL: i32 = 20;
const MIN_VAR_FREQ: f64 = 0.10;
const MIN_FREQ_FOR_HOM: f64 = 0.75;
const NORMAL_PURITY: f64 = 1.0;
const TUMOR_PURITY: f64 = 1.0;
const P_VALUE: f64 = 0.99;
const SOMATIC_P_VALUE: f64 = 0.05;
const STRAND_FILTER: i32 = 1;

const MIN_TUMOR_FREQ: f64 = 0.10;
const MAX_NORMAL_FREQ: f64 = 0.05;
const PROCESS_P_VALUE: f64 = 0.05;

const CNV_MIN_COVERAGE: i32 = 20;
const CNV_P_VALUE: f64 = 0.01;
const MIN_SEGMENT_SIZE: i32 = 10;
const MAX_SEGMENT_SIZE: i32 = 100;

const CNV_AMP_THRESHOLD: f64 = 0.25;
const CNV_DEL_THRESHOLD: f64 = 0.25;
const CNV_RECENTER_UP: i32 = 0;
const CNV_RECENTER_DOWN: i32 = 0;

const BRC_MAP_QUAL: i32 = 10;
const BRC_BASE_QUAL: i32 = 15;
const MAX_PARALLEL_JOBS: usize = 30;

type AppResult<T> = Result<T, String>;

#[derive(Clone)]
struct Paths {
    bam_dir: PathBuf,
    resume_state_dir: PathBuf,
    flagstat_dir: PathBuf,
    mpileup_dir: PathBuf,
    somatic_dir: PathBuf,
    copy_number_dir: PathBuf,
    snp_var_dir: PathBuf,
    indel_var_dir: PathBuf,
    readcount_dir: PathBuf,
    filtered_dir: PathBuf,
    summary_file: PathBuf,
}

#[derive(Clone)]
struct Args {
    from_stage: u8,
    to_stage: u8,
    resume: bool,
    dry_run: bool,
}

fn now_string() -> String {
    let output = Command::new("date")
        .arg("+%Y-%m-%d %H:%M:%S")
        .output();
    match output {
        Ok(o) if o.status.success() => String::from_utf8_lossy(&o.stdout).trim().to_string(),
        _ => "unknown-time".to_string(),
    }
}

fn log_message(msg: &str) {
    println!("[{}] {}", now_string(), msg);
}

fn create_directory(path: &Path) -> AppResult<()> {
    if !path.exists() {
        fs::create_dir_all(path).map_err(|e| format!("mkdir failed {}: {}", path.display(), e))?;
        let mut perm = fs::metadata(path)
            .map_err(|e| format!("metadata failed {}: {}", path.display(), e))?
            .permissions();
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            perm.set_mode(0o755);
            fs::set_permissions(path, perm)
                .map_err(|e| format!("chmod failed {}: {}", path.display(), e))?;
        }
        log_message(&format!("Created: {}", path.display()));
    }
    Ok(())
}

fn check_file_exists(path: &Path) -> AppResult<()> {
    if !path.is_file() {
        return Err(format!("Required file not found: {}", path.display()));
    }
    Ok(())
}

fn check_command_exists(cmd: &str) -> AppResult<()> {
    let status = Command::new("sh")
        .arg("-lc")
        .arg(format!("command -v {} >/dev/null 2>&1", cmd))
        .status()
        .map_err(|e| format!("failed checking command {}: {}", cmd, e))?;
    if status.success() {
        Ok(())
    } else {
        Err(format!("Required command not found on PATH: {}", cmd))
    }
}

fn get_sample_name(entry: &str) -> String {
    if entry.ends_with(PAIRS_SUFFIX) {
        entry.trim_end_matches(PAIRS_SUFFIX).to_string()
    } else {
        entry.to_string()
    }
}

fn get_bam_path(paths: &Paths, entry: &str) -> PathBuf {
    let mut p = paths.bam_dir.clone();
    p.push(format!("{}{}", get_sample_name(entry), BAM_SUFFIX));
    p
}

fn check_bam_index(bam: &Path) -> AppResult<()> {
    let bam_bai = PathBuf::from(format!("{}.bai", bam.display()));
    let alt_bai = if let Some(stem) = bam.to_string_lossy().strip_suffix(".bam") {
        PathBuf::from(format!("{}.bai", stem))
    } else {
        PathBuf::from("")
    };
    if bam_bai.is_file() || alt_bai.is_file() {
        Ok(())
    } else {
        Err(format!(
            "BAM index not found for: {} (expected {} or {})",
            bam.display(),
            bam_bai.display(),
            alt_bai.display()
        ))
    }
}

fn clean_pair_field(s: &str) -> String {
    s.chars().filter(|c| *c != ' ' && *c != '\r').collect()
}

fn read_pairs() -> AppResult<Vec<(String, String)>> {
    let file = File::open(FILE_PAIRS_LIST).map_err(|e| format!("open {}: {}", FILE_PAIRS_LIST, e))?;
    let reader = BufReader::new(file);
    let mut pairs = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("read {}: {}", FILE_PAIRS_LIST, e))?;
        let mut parts = line.splitn(3, ',');
        let p1 = clean_pair_field(parts.next().unwrap_or(""));
        let p2 = clean_pair_field(parts.next().unwrap_or(""));
        if p1.is_empty() || p2.is_empty() {
            log_message(&format!("WARNING: Incomplete pair (empty field) — skipping: '{}','{}'", p1, p2));
            continue;
        }
        pairs.push((p1, p2));
    }
    Ok(pairs)
}

fn glob_match(pattern: &str, name: &str) -> bool {
    if pattern == "*" {
        return true;
    }
    if !pattern.contains('*') {
        return pattern == name;
    }
    let parts: Vec<&str> = pattern.split('*').collect();
    let mut idx = 0usize;
    let starts_with_star = pattern.starts_with('*');
    let ends_with_star = pattern.ends_with('*');

    for (i, part) in parts.iter().enumerate() {
        if part.is_empty() {
            continue;
        }
        if i == 0 && !starts_with_star {
            if !name[idx..].starts_with(part) {
                return false;
            }
            idx += part.len();
            continue;
        }
        if i == parts.len() - 1 && !ends_with_star {
            if !name.ends_with(part) {
                return false;
            }
            continue;
        }
        if let Some(pos) = name[idx..].find(part) {
            idx += pos + part.len();
        } else {
            return false;
        }
    }
    true
}

fn list_files_matching(dir: &Path, pattern: &str) -> AppResult<Vec<PathBuf>> {
    if !dir.is_dir() {
        return Ok(Vec::new());
    }
    let mut files = Vec::new();
    for entry in fs::read_dir(dir).map_err(|e| format!("read_dir {}: {}", dir.display(), e))? {
        let entry = entry.map_err(|e| format!("read_dir entry {}: {}", dir.display(), e))?;
        let path = entry.path();
        if !path.is_file() {
            continue;
        }
        let name = path
            .file_name()
            .and_then(OsStr::to_str)
            .ok_or_else(|| format!("non-utf8 filename in {}", dir.display()))?;
        if glob_match(pattern, name) {
            files.push(path);
        }
    }
    files.sort();
    Ok(files)
}

fn stage_output_exists(stage: u8, paths: &Paths) -> AppResult<bool> {
    let yes = match stage {
        1 => !list_files_matching(&paths.flagstat_dir, "*.flagstats")?.is_empty(),
        2 => !list_files_matching(&paths.mpileup_dir, "*.mpileup")?.is_empty(),
        3 => !list_files_matching(&paths.somatic_dir, "*.snp.vcf")?.is_empty(),
        4 => !list_files_matching(&paths.somatic_dir, "*.hc.vcf")?.is_empty(),
        5 => !list_files_matching(&paths.copy_number_dir, "*.copynumber")?.is_empty(),
        6 => !list_files_matching(&paths.copy_number_dir, "*.copynumber.called")?.is_empty(),
        7 => !list_files_matching(&paths.snp_var_dir, "*.var")?.is_empty(),
        8 => !list_files_matching(&paths.readcount_dir, "*.readcount")?.is_empty(),
        9 => !list_files_matching(&paths.filtered_dir, "*.fpfilter.vcf")?.is_empty(),
        10 => paths.summary_file.is_file(),
        _ => false,
    };
    Ok(yes)
}

fn sha256_of_string(s: &str) -> AppResult<String> {
    let mut hasher = Sha256::new();
    hasher.update(s.as_bytes());
    Ok(format!("{:x}", hasher.finalize()))
}

fn sha256_of_file_or_missing(path: &Path) -> AppResult<String> {
    if !path.is_file() {
        return Ok("MISSING".to_string());
    }
    let file = File::open(path).map_err(|e| format!("open {}: {}", path.display(), e))?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();
    let mut buf = [0u8; 65536];
    loop {
        let n = reader
            .read(&mut buf)
            .map_err(|e| format!("read {}: {}", path.display(), e))?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn modified_epoch(meta: &fs::Metadata) -> AppResult<f64> {
    let t = meta
        .modified()
        .map_err(|e| format!("modified time error: {}", e))?
        .duration_since(UNIX_EPOCH)
        .map_err(|e| format!("time before epoch: {}", e))?;
    Ok(t.as_secs_f64())
}

fn hash_manifest(dir: &Path, pattern: &str) -> AppResult<String> {
    if !dir.is_dir() {
        return Ok("MISSING_DIR".to_string());
    }
    let files = list_files_matching(dir, pattern)?;
    if files.is_empty() {
        return Ok("EMPTY".to_string());
    }

    let mut lines = Vec::new();
    for p in files {
        let name = p
            .file_name()
            .and_then(OsStr::to_str)
            .ok_or_else(|| format!("non-utf8 filename {}", p.display()))?;
        let meta = fs::metadata(&p).map_err(|e| format!("metadata {}: {}", p.display(), e))?;
        let mtime = modified_epoch(&meta)?;
        lines.push(format!("{}|{}|{}", name, meta.len(), mtime));
    }
    lines.sort();
    sha256_of_string(&lines.join("\n"))
}

fn hash_pairs_cleaned() -> AppResult<String> {
    if !Path::new(FILE_PAIRS_LIST).is_file() {
        return Ok("MISSING".to_string());
    }
    let pairs = read_pairs()?;
    if pairs.is_empty() {
        return Ok("EMPTY".to_string());
    }
    let normalized = pairs
        .iter()
        .map(|(t, n)| format!("{},{}", t, n))
        .collect::<Vec<_>>()
        .join("\n");
    sha256_of_string(&normalized)
}

fn hash_paired_bam_metadata(paths: &Paths) -> AppResult<String> {
    if !Path::new(FILE_PAIRS_LIST).is_file() {
        return Ok("MISSING".to_string());
    }

    let pairs = read_pairs()?;
    if pairs.is_empty() {
        return Ok("EMPTY".to_string());
    }

    let mut lines = Vec::new();
    for (entry1, entry2) in pairs {
        let normal_bam = get_bam_path(paths, &entry2);
        let tumor_bam = get_bam_path(paths, &entry1);
        for bam in [&normal_bam, &tumor_bam] {
            if bam.is_file() {
                let m = fs::metadata(bam).map_err(|e| format!("metadata {}: {}", bam.display(), e))?;
                let modified = m
                    .modified()
                    .map_err(|e| format!("modified {}: {}", bam.display(), e))?
                    .duration_since(UNIX_EPOCH)
                    .map_err(|e| format!("time {}: {}", bam.display(), e))?
                    .as_secs();
                lines.push(format!("{}|{}|{}", bam.display(), m.len(), modified));
            } else {
                lines.push(format!("{}|MISSING", bam.display()));
            }
        }
    }
    sha256_of_string(&lines.join("\n"))
}

fn file_sha(path_str: &str) -> AppResult<String> {
    if path_str.is_empty() {
        return Ok("MISSING".to_string());
    }
    sha256_of_file_or_missing(Path::new(path_str))
}

fn compute_stage_hash(stage: u8, paths: &Paths) -> AppResult<String> {
    let exe = env::current_exe().map_err(|e| format!("current_exe: {}", e))?;
    let mut payload = String::new();
    payload.push_str(&format!("exe_sha256={}\n", sha256_of_file_or_missing(&exe)?));
    payload.push_str(&format!("pairs_sha256={}\n", hash_pairs_cleaned()?));
    payload.push_str(&format!("bam_sha256={}\n", hash_paired_bam_metadata(paths)?));
    payload.push_str(&format!("stage={}\n", stage));

    match stage {
        1 => {
            payload.push_str(&format!("bam_dir={}\n", paths.bam_dir.display()));
            payload.push_str(&format!("bam_suffix={}\n", BAM_SUFFIX));
        }
        2 => {
            payload.push_str(&format!("reference_sha256={}\n", file_sha(GENOMEIDX1)?));
            payload.push_str(&format!("target_bed_sha256={}\n", file_sha(TARGET_BED)?));
            payload.push_str(&format!("BRC_MAP_QUAL={}\n", BRC_MAP_QUAL));
            payload.push_str(&format!("MIN_BASE_QUAL={}\n", MIN_BASE_QUAL));
        }
        3 => {
            payload.push_str(&format!(
                "jar_sha256={}\n",
                file_sha(&format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR))?
            ));
            payload.push_str(&format!("mpileup_manifest={}\n", hash_manifest(&paths.mpileup_dir, "*.mpileup")?));
            payload.push_str(&format!("vcf_sample_sha={}\n", file_sha(VCF_SAMPLE_LIST)?));
            payload.push_str(&format!("MIN_COVERAGE={}\n", MIN_COVERAGE));
            payload.push_str(&format!("MIN_COVERAGE_NORMAL={}\n", MIN_COVERAGE_NORMAL));
            payload.push_str(&format!("MIN_COVERAGE_TUMOR={}\n", MIN_COVERAGE_TUMOR));
            payload.push_str(&format!("MIN_VAR_FREQ={}\n", MIN_VAR_FREQ));
            payload.push_str(&format!("MIN_FREQ_FOR_HOM={}\n", MIN_FREQ_FOR_HOM));
            payload.push_str(&format!("NORMAL_PURITY={}\n", NORMAL_PURITY));
            payload.push_str(&format!("TUMOR_PURITY={}\n", TUMOR_PURITY));
            payload.push_str(&format!("P_VALUE={}\n", P_VALUE));
            payload.push_str(&format!("SOMATIC_P_VALUE={}\n", SOMATIC_P_VALUE));
            payload.push_str(&format!("STRAND_FILTER={}\n", STRAND_FILTER));
            payload.push_str(&format!("BRC_MAP_QUAL={}\n", BRC_MAP_QUAL));
            payload.push_str(&format!("MIN_BASE_QUAL={}\n", MIN_BASE_QUAL));
        }
        4 => {
            payload.push_str(&format!(
                "jar_sha256={}\n",
                file_sha(&format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR))?
            ));
            payload.push_str(&format!("snp_manifest={}\n", hash_manifest(&paths.somatic_dir, "*.snp.vcf")?));
            payload.push_str(&format!("indel_manifest={}\n", hash_manifest(&paths.somatic_dir, "*.indel.vcf")?));
            payload.push_str(&format!("MIN_TUMOR_FREQ={}\n", MIN_TUMOR_FREQ));
            payload.push_str(&format!("MAX_NORMAL_FREQ={}\n", MAX_NORMAL_FREQ));
            payload.push_str(&format!("PROCESS_P_VALUE={}\n", PROCESS_P_VALUE));
        }
        5 => {
            payload.push_str(&format!(
                "jar_sha256={}\n",
                file_sha(&format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR))?
            ));
            payload.push_str(&format!("mpileup_manifest={}\n", hash_manifest(&paths.mpileup_dir, "*.mpileup")?));
            payload.push_str(&format!("flagstats_manifest={}\n", hash_manifest(&paths.flagstat_dir, "*.flagstats")?));
            payload.push_str(&format!("CNV_MIN_COVERAGE={}\n", CNV_MIN_COVERAGE));
            payload.push_str(&format!("MIN_BASE_QUAL={}\n", MIN_BASE_QUAL));
            payload.push_str(&format!("BRC_MAP_QUAL={}\n", BRC_MAP_QUAL));
            payload.push_str(&format!("MIN_SEGMENT_SIZE={}\n", MIN_SEGMENT_SIZE));
            payload.push_str(&format!("MAX_SEGMENT_SIZE={}\n", MAX_SEGMENT_SIZE));
            payload.push_str(&format!("CNV_P_VALUE={}\n", CNV_P_VALUE));
        }
        6 => {
            payload.push_str(&format!(
                "jar_sha256={}\n",
                file_sha(&format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR))?
            ));
            payload.push_str(&format!("cnv_manifest={}\n", hash_manifest(&paths.copy_number_dir, "*.copynumber")?));
            payload.push_str(&format!("CNV_AMP_THRESHOLD={}\n", CNV_AMP_THRESHOLD));
            payload.push_str(&format!("CNV_DEL_THRESHOLD={}\n", CNV_DEL_THRESHOLD));
            payload.push_str(&format!("CNV_RECENTER_UP={}\n", CNV_RECENTER_UP));
            payload.push_str(&format!("CNV_RECENTER_DOWN={}\n", CNV_RECENTER_DOWN));
        }
        7 => {
            payload.push_str(&format!("hc_manifest={}\n", hash_manifest(&paths.somatic_dir, "*.hc.vcf")?));
        }
        8 => {
            payload.push_str(&format!("reference_sha256={}\n", file_sha(GENOMEIDX1)?));
            payload.push_str(&format!("snp_var_manifest={}\n", hash_manifest(&paths.snp_var_dir, "*.var")?));
            payload.push_str(&format!("indel_var_manifest={}\n", hash_manifest(&paths.indel_var_dir, "*.var")?));
            payload.push_str(&format!("BRC_MAP_QUAL={}\n", BRC_MAP_QUAL));
            payload.push_str(&format!("BRC_BASE_QUAL={}\n", BRC_BASE_QUAL));
        }
        9 => {
            payload.push_str(&format!("fpfilter_sha256={}\n", file_sha(&format!("{}/fpfilter.pl", SCRIPTSDIR))?));
            payload.push_str(&format!("readcount_manifest={}\n", hash_manifest(&paths.readcount_dir, "*.readcount")?));
            payload.push_str(&format!("somatic_manifest={}\n", hash_manifest(&paths.somatic_dir, "*.vcf")?));
            payload.push_str(&format!("MIN_VAR_FREQ={}\n", MIN_VAR_FREQ));
            payload.push_str(&format!("BRC_BASE_QUAL={}\n", BRC_BASE_QUAL));
        }
        10 => {
            payload.push_str(&format!("somatic_manifest={}\n", hash_manifest(&paths.somatic_dir, "*.vcf")?));
            payload.push_str(&format!("copy_manifest={}\n", hash_manifest(&paths.copy_number_dir, "*.copynumber*")?));
            payload.push_str(&format!("filtered_manifest={}\n", hash_manifest(&paths.filtered_dir, "*.fpfilter.vcf")?));
        }
        _ => {}
    }

    sha256_of_string(&payload)
}

fn stage_marker_path(paths: &Paths, stage: u8) -> PathBuf {
    paths.resume_state_dir.join(format!("stage_{}.sha256", stage))
}

fn write_stage_marker(paths: &Paths, stage: u8) -> AppResult<()> {
    let marker = stage_marker_path(paths, stage);
    let hash = compute_stage_hash(stage, paths)?;
    fs::write(&marker, format!("{}\n", hash)).map_err(|e| format!("write marker {}: {}", marker.display(), e))
}

fn resume_stage_match(paths: &Paths, stage: u8) -> AppResult<bool> {
    let marker = stage_marker_path(paths, stage);
    if !marker.is_file() {
        return Ok(false);
    }
    let expected = fs::read_to_string(&marker)
        .map_err(|e| format!("read marker {}: {}", marker.display(), e))?
        .trim()
        .to_string();
    let current = compute_stage_hash(stage, paths)?;
    if expected == current {
        stage_output_exists(stage, paths)
    } else {
        Ok(false)
    }
}

fn spawn_command(program: &str, args: &[String], stdout_file: Option<&Path>) -> AppResult<Child> {
    let mut cmd = Command::new(program);
    cmd.args(args);
    if let Some(out) = stdout_file {
        let f = File::create(out).map_err(|e| format!("create {}: {}", out.display(), e))?;
        cmd.stdout(Stdio::from(f));
    }
    cmd.stderr(Stdio::inherit());
    cmd.spawn().map_err(|e| format!("spawn {} failed: {}", program, e))
}

fn wait_all(children: &mut Vec<Child>) -> AppResult<()> {
    let mut failed = false;
    for child in children.iter_mut() {
        let status = child.wait().map_err(|e| format!("wait failed: {}", e))?;
        if !status.success() {
            failed = true;
        }
    }
    children.clear();
    if failed {
        Err("one or more subprocesses failed".to_string())
    } else {
        Ok(())
    }
}

fn push_child(children: &mut Vec<Child>, child: Child, limit: usize) -> AppResult<()> {
    children.push(child);
    if children.len() >= limit {
        wait_all(children)?;
    }
    Ok(())
}

fn setup_directories(paths: &Paths) -> AppResult<()> {
    log_message("Setting up output directory structure...");
    for d in [
        &paths.flagstat_dir,
        &paths.mpileup_dir,
        &paths.somatic_dir,
        &paths.copy_number_dir,
        &paths.snp_var_dir,
        &paths.indel_var_dir,
        &paths.readcount_dir,
        &paths.filtered_dir,
        &paths.resume_state_dir,
    ] {
        create_directory(d)?;
    }
    log_message("Directory setup complete");
    Ok(())
}

fn generate_flagstats(paths: &Paths) -> AppResult<()> {
    log_message("=== STAGE 1: Generating BAM flagstats ===");
    if !paths.bam_dir.is_dir() {
        return Err(format!("BAM directory not found: {}", paths.bam_dir.display()));
    }

    let mut children = Vec::new();
    let mut bams = list_files_matching(&paths.bam_dir, &format!("*{}", BAM_SUFFIX))?;
    bams.sort();
    for bam in bams {
        let sample = bam
            .file_name()
            .and_then(OsStr::to_str)
            .unwrap_or_default()
            .trim_end_matches(BAM_SUFFIX)
            .to_string();
        log_message(&format!("Flagstats: {}", sample));
        let out = paths.flagstat_dir.join(format!("{}.flagstats", sample));
        let args = vec!["flagstat".to_string(), bam.to_string_lossy().to_string()];
        let child = spawn_command("samtools", &args, Some(&out))?;
        push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
    }
    wait_all(&mut children)?;
    log_message("Flagstats complete");
    Ok(())
}

fn generate_mpileup(paths: &Paths) -> AppResult<()> {
    log_message("=== STAGE 2: Generating paired mpileup (normal first, tumour second) ===");
    check_file_exists(Path::new(FILE_PAIRS_LIST))?;
    check_file_exists(Path::new(GENOMEIDX1))?;
    if !TARGET_BED.is_empty() {
        check_file_exists(Path::new(TARGET_BED))?;
    }

    let pairs = read_pairs()?;
    let mut children = Vec::new();

    for (entry1, entry2) in pairs {
        let samplet = get_sample_name(&entry1);
        let samplen = get_sample_name(&entry2);
        let tumor_bam = get_bam_path(paths, &entry1);
        let normal_bam = get_bam_path(paths, &entry2);
        check_file_exists(&normal_bam)?;
        check_file_exists(&tumor_bam)?;

        let out = paths.mpileup_dir.join(format!("{}_{}.mpileup", samplen, samplet));
        log_message(&format!("Mpileup: Normal={}  Tumour={}", samplen, samplet));

        let mut args = vec![
            "mpileup".to_string(),
            "-B".to_string(),
            "-q".to_string(),
            BRC_MAP_QUAL.to_string(),
            "-Q".to_string(),
            MIN_BASE_QUAL.to_string(),
            "-F".to_string(),
            "0x400".to_string(),
        ];
        if !TARGET_BED.is_empty() {
            args.push("-l".to_string());
            args.push(TARGET_BED.to_string());
        }
        args.push("-f".to_string());
        args.push(GENOMEIDX1.to_string());
        args.push(normal_bam.to_string_lossy().to_string());
        args.push(tumor_bam.to_string_lossy().to_string());

        let child = spawn_command("samtools", &args, Some(&out))?;
        push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
    }

    wait_all(&mut children)?;
    log_message("Mpileup complete");
    Ok(())
}

fn run_varscan_somatic(paths: &Paths) -> AppResult<()> {
    log_message("=== STAGE 3: VarScan somatic variant calling ===");
    check_file_exists(Path::new(FILE_PAIRS_LIST))?;
    check_file_exists(Path::new(&format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR)))?;

    let pairs = read_pairs()?;
    let mut children = Vec::new();

    for (entry1, entry2) in pairs {
        let samplet = get_sample_name(&entry1);
        let samplen = get_sample_name(&entry2);
        let mpileup_file = paths.mpileup_dir.join(format!("{}_{}.mpileup", samplen, samplet));
        if !mpileup_file.is_file() {
            log_message(&format!("WARNING: Mpileup not found: {} — skipping pair", mpileup_file.display()));
            continue;
        }

        log_message(&format!("VarScan somatic: Normal={}  Tumour={}", samplen, samplet));
        let mut args = vec![
            "-jar".to_string(),
            format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR),
            "somatic".to_string(),
            mpileup_file.to_string_lossy().to_string(),
            paths.somatic_dir.join(format!("{}_{}", samplen, samplet)).to_string_lossy().to_string(),
            "--mpileup".to_string(),
            "1".to_string(),
            "--min-coverage".to_string(),
            MIN_COVERAGE.to_string(),
            "--min-coverage-normal".to_string(),
            MIN_COVERAGE_NORMAL.to_string(),
            "--min-coverage-tumor".to_string(),
            MIN_COVERAGE_TUMOR.to_string(),
            "--min-var-freq".to_string(),
            MIN_VAR_FREQ.to_string(),
            "--min-freq-for-hom".to_string(),
            MIN_FREQ_FOR_HOM.to_string(),
            "--normal-purity".to_string(),
            NORMAL_PURITY.to_string(),
            "--tumor-purity".to_string(),
            TUMOR_PURITY.to_string(),
            "--p-value".to_string(),
            P_VALUE.to_string(),
            "--somatic-p-value".to_string(),
            SOMATIC_P_VALUE.to_string(),
            "--strand-filter".to_string(),
            STRAND_FILTER.to_string(),
            "--min-map-qual".to_string(),
            BRC_MAP_QUAL.to_string(),
            "--min-base-qual".to_string(),
            MIN_BASE_QUAL.to_string(),
            "--output-vcf".to_string(),
            "1".to_string(),
        ];

        if !VCF_SAMPLE_LIST.is_empty() && Path::new(VCF_SAMPLE_LIST).is_file() {
            args.push("--vcf-sample-list".to_string());
            args.push(VCF_SAMPLE_LIST.to_string());
        }

        let child = spawn_command("java", &args, None)?;
        push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
    }

    wait_all(&mut children)?;
    log_message("VarScan somatic calling complete");
    Ok(())
}

fn process_somatic_variants(paths: &Paths) -> AppResult<()> {
    log_message("=== STAGE 4: processSomatic — classifying into Somatic / Germline / LOH ===");
    let mut children = Vec::new();

    for snp_file in list_files_matching(&paths.somatic_dir, "*.snp.vcf")? {
        log_message(&format!("processSomatic SNP: {}", snp_file.file_name().and_then(OsStr::to_str).unwrap_or("")));
        let args = vec![
            "-jar".to_string(),
            format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR),
            "processSomatic".to_string(),
            snp_file.to_string_lossy().to_string(),
            "--min-tumor-freq".to_string(),
            MIN_TUMOR_FREQ.to_string(),
            "--max-normal-freq".to_string(),
            MAX_NORMAL_FREQ.to_string(),
            "--p-value".to_string(),
            PROCESS_P_VALUE.to_string(),
        ];
        let child = spawn_command("java", &args, None)?;
        push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
    }

    for indel_file in list_files_matching(&paths.somatic_dir, "*.indel.vcf")? {
        log_message(&format!("processSomatic INDEL: {}", indel_file.file_name().and_then(OsStr::to_str).unwrap_or("")));
        let args = vec![
            "-jar".to_string(),
            format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR),
            "processSomatic".to_string(),
            indel_file.to_string_lossy().to_string(),
            "--min-tumor-freq".to_string(),
            MIN_TUMOR_FREQ.to_string(),
            "--max-normal-freq".to_string(),
            MAX_NORMAL_FREQ.to_string(),
            "--p-value".to_string(),
            PROCESS_P_VALUE.to_string(),
        ];
        let child = spawn_command("java", &args, None)?;
        push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
    }

    wait_all(&mut children)?;
    log_message("processSomatic complete");
    Ok(())
}

fn check_stage4_output(paths: &Paths) -> AppResult<()> {
    let hc = list_files_matching(&paths.somatic_dir, "*.hc.vcf")?;
    if hc.is_empty() {
        return Err("processSomatic (Stage 4) produced no .hc.vcf files".to_string());
    }
    log_message(&format!("Stage 4 guard passed: {} .hc.vcf file(s)", hc.len()));
    Ok(())
}

fn parse_primary_mapped_count(flagstats: &Path) -> AppResult<Option<u64>> {
    let file = File::open(flagstats).map_err(|e| format!("open {}: {}", flagstats.display(), e))?;
    for line in BufReader::new(file).lines() {
        let line = line.map_err(|e| format!("read {}: {}", flagstats.display(), e))?;
        if line.contains("primary mapped") {
            let n = line.split_whitespace().next().unwrap_or("");
            if let Ok(v) = n.parse::<u64>() {
                return Ok(Some(v));
            }
        }
    }
    Ok(None)
}

fn run_varscan_copynumber(paths: &Paths) -> AppResult<()> {
    log_message("=== STAGE 5: VarScan copy number analysis ===");
    check_file_exists(Path::new(FILE_PAIRS_LIST))?;

    let pairs = read_pairs()?;
    let mut children = Vec::new();

    for (entry1, entry2) in pairs {
        let samplet = get_sample_name(&entry1);
        let samplen = get_sample_name(&entry2);

        let normal_flags = paths.flagstat_dir.join(format!("{}.flagstats", samplen));
        let tumor_flags = paths.flagstat_dir.join(format!("{}.flagstats", samplet));

        let mut dataratio = 1.0f64;
        if normal_flags.is_file() && tumor_flags.is_file() {
            let n_mapped = parse_primary_mapped_count(&normal_flags)?;
            let t_mapped = parse_primary_mapped_count(&tumor_flags)?;
            if let (Some(n), Some(t)) = (n_mapped, t_mapped) {
                if t != 0 {
                    dataratio = n as f64 / t as f64;
                } else {
                    log_message(&format!("WARNING: t_mapped=0 for {} — using ratio=1.0", samplet));
                }
            } else {
                log_message(&format!("WARNING: Could not parse primary mapped for {}/{} — ratio=1.0", samplen, samplet));
            }
        } else {
            log_message(&format!("WARNING: Flagstats missing for {} or {} — ratio=1.0", samplen, samplet));
        }

        log_message(&format!("Data ratio {}/{}: {:.6}", samplen, samplet, dataratio));

        let mpileup_file = paths.mpileup_dir.join(format!("{}_{}.mpileup", samplen, samplet));
        if !mpileup_file.is_file() {
            log_message(&format!("WARNING: Mpileup not found: {} — skipping pair", mpileup_file.display()));
            continue;
        }

        let args = vec![
            "-jar".to_string(),
            format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR),
            "copynumber".to_string(),
            mpileup_file.to_string_lossy().to_string(),
            paths.copy_number_dir.join(format!("{}_{}", samplen, samplet)).to_string_lossy().to_string(),
            "--mpileup".to_string(),
            "1".to_string(),
            "--min-coverage".to_string(),
            CNV_MIN_COVERAGE.to_string(),
            "--min-base-qual".to_string(),
            MIN_BASE_QUAL.to_string(),
            "--min-map-qual".to_string(),
            BRC_MAP_QUAL.to_string(),
            "--min-segment-size".to_string(),
            MIN_SEGMENT_SIZE.to_string(),
            "--max-segment-size".to_string(),
            MAX_SEGMENT_SIZE.to_string(),
            "--p-value".to_string(),
            CNV_P_VALUE.to_string(),
            "--data-ratio".to_string(),
            format!("{:.6}", dataratio),
        ];
        let child = spawn_command("java", &args, None)?;
        push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
    }

    wait_all(&mut children)?;
    log_message("VarScan copy number complete");
    Ok(())
}

fn run_copy_caller(paths: &Paths) -> AppResult<()> {
    log_message("=== STAGE 6: VarScan copyCaller ===");
    let mut children = Vec::new();

    for cnv_file in list_files_matching(&paths.copy_number_dir, "*.copynumber")? {
        let base = cnv_file
            .file_name()
            .and_then(OsStr::to_str)
            .unwrap_or_default()
            .trim_end_matches(".copynumber")
            .to_string();
        log_message(&format!("copyCaller: {}", base));

        let args = vec![
            "-jar".to_string(),
            format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR),
            "copyCaller".to_string(),
            cnv_file.to_string_lossy().to_string(),
            "--output-file".to_string(),
            paths.copy_number_dir.join(format!("{}.copynumber.called", base)).to_string_lossy().to_string(),
            "--output-homdel-file".to_string(),
            paths.copy_number_dir.join(format!("{}.copynumber.homdel", base)).to_string_lossy().to_string(),
            "--amp-threshold".to_string(),
            CNV_AMP_THRESHOLD.to_string(),
            "--del-threshold".to_string(),
            CNV_DEL_THRESHOLD.to_string(),
            "--recenter-up".to_string(),
            CNV_RECENTER_UP.to_string(),
            "--recenter-down".to_string(),
            CNV_RECENTER_DOWN.to_string(),
        ];

        let child = spawn_command("java", &args, None)?;
        push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
    }

    wait_all(&mut children)?;
    log_message("copyCaller complete");
    Ok(())
}

fn prepare_filter_input(paths: &Paths) -> AppResult<()> {
    log_message("=== STAGE 7: Preparing VAR position files for bam-readcount ===");

    for hc_file in list_files_matching(&paths.somatic_dir, "*.snp.*.hc.vcf")? {
        let base = hc_file
            .file_name()
            .and_then(OsStr::to_str)
            .unwrap_or_default()
            .trim_end_matches(".vcf")
            .to_string();
        let out_var = paths.snp_var_dir.join(format!("{}.var", base));
        log_message(&format!("SNP VAR: {}", base));

        let in_f = File::open(&hc_file).map_err(|e| format!("open {}: {}", hc_file.display(), e))?;
        let mut out_f = File::create(&out_var).map_err(|e| format!("create {}: {}", out_var.display(), e))?;
        for line in BufReader::new(in_f).lines() {
            let line = line.map_err(|e| format!("read {}: {}", hc_file.display(), e))?;
            if line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() >= 2 {
                writeln!(out_f, "{}\t{}\t{}", cols[0], cols[1], cols[1])
                    .map_err(|e| format!("write {}: {}", out_var.display(), e))?;
            }
        }
    }

    for hc_file in list_files_matching(&paths.somatic_dir, "*.indel.Somatic.hc.vcf")? {
        let base = hc_file
            .file_name()
            .and_then(OsStr::to_str)
            .unwrap_or_default()
            .trim_end_matches(".vcf")
            .to_string();
        let out_var = paths.indel_var_dir.join(format!("{}.var", base));
        log_message(&format!("INDEL VAR: {}", base));

        let in_f = File::open(&hc_file).map_err(|e| format!("open {}: {}", hc_file.display(), e))?;
        let mut out_f = File::create(&out_var).map_err(|e| format!("create {}: {}", out_var.display(), e))?;
        for line in BufReader::new(in_f).lines() {
            let line = line.map_err(|e| format!("read {}: {}", hc_file.display(), e))?;
            if line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() >= 5 {
                let pos: u64 = cols[1].parse().map_err(|e| {
                    format!(
                        "malformed VCF position '{}' in {}: {}",
                        cols[1],
                        hc_file.display(),
                        e
                    )
                })?;
                let is_del = cols[3].len() > cols[4].len();
                let shifted = if is_del { pos + 1 } else { pos };
                writeln!(out_f, "{}\t{}\t{}", cols[0], shifted, shifted)
                    .map_err(|e| format!("write {}: {}", out_var.display(), e))?;
            }
        }
    }

    log_message("VAR file preparation complete");
    Ok(())
}

fn run_bam_readcount(paths: &Paths) -> AppResult<()> {
    log_message("=== STAGE 8: bam-readcount ===");
    check_file_exists(Path::new(FILE_PAIRS_LIST))?;

    let pairs = read_pairs()?;
    let mut children = Vec::new();

    for (entry1, entry2) in pairs {
        let samplet = get_sample_name(&entry1);
        let samplen = get_sample_name(&entry2);
        let tumor_bam = get_bam_path(paths, &entry1);
        let normal_bam = get_bam_path(paths, &entry2);
        let prefix = format!("{}_{}", samplen, samplet);

        check_bam_index(&tumor_bam)?;
        check_bam_index(&normal_bam)?;

        let var_file = paths.snp_var_dir.join(format!("{}.snp.Somatic.hc.var", prefix));
        if var_file.is_file() {
            log_message(&format!("bam-readcount: Somatic SNP — tumour {}", samplet));
            let args = vec![
                "-q".to_string(),
                BRC_MAP_QUAL.to_string(),
                "-b".to_string(),
                BRC_BASE_QUAL.to_string(),
                "-f".to_string(),
                GENOMEIDX1.to_string(),
                "-l".to_string(),
                var_file.to_string_lossy().to_string(),
                tumor_bam.to_string_lossy().to_string(),
            ];
            let out = paths.readcount_dir.join(format!("{}.snp.Somatic.hc.readcount", prefix));
            let child = spawn_command("bam-readcount", &args, Some(&out))?;
            push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
        }

        let indel_var = paths.indel_var_dir.join(format!("{}.indel.Somatic.hc.var", prefix));
        if indel_var.is_file() {
            log_message(&format!("bam-readcount: Somatic INDEL — tumour {}", samplet));
            let args = vec![
                "-q".to_string(),
                BRC_MAP_QUAL.to_string(),
                "-b".to_string(),
                BRC_BASE_QUAL.to_string(),
                "-f".to_string(),
                GENOMEIDX1.to_string(),
                "-l".to_string(),
                indel_var.to_string_lossy().to_string(),
                tumor_bam.to_string_lossy().to_string(),
            ];
            let out = paths.readcount_dir.join(format!("{}.indel.Somatic.hc.readcount", prefix));
            let child = spawn_command("bam-readcount", &args, Some(&out))?;
            push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
        }

        let germ_var = paths.snp_var_dir.join(format!("{}.snp.Germline.hc.var", prefix));
        if germ_var.is_file() {
            log_message(&format!("bam-readcount: Germline SNP — normal {}", samplen));
            let args = vec![
                "-q".to_string(),
                BRC_MAP_QUAL.to_string(),
                "-b".to_string(),
                BRC_BASE_QUAL.to_string(),
                "-f".to_string(),
                GENOMEIDX1.to_string(),
                "-l".to_string(),
                germ_var.to_string_lossy().to_string(),
                normal_bam.to_string_lossy().to_string(),
            ];
            let out = paths.readcount_dir.join(format!("{}.snp.Germline.hc.readcount", prefix));
            let child = spawn_command("bam-readcount", &args, Some(&out))?;
            push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
        }

        let loh_var = paths.snp_var_dir.join(format!("{}.snp.LOH.hc.var", prefix));
        if loh_var.is_file() {
            log_message(&format!("bam-readcount: LOH SNP — tumour {}", samplet));
            let args = vec![
                "-q".to_string(),
                BRC_MAP_QUAL.to_string(),
                "-b".to_string(),
                BRC_BASE_QUAL.to_string(),
                "-f".to_string(),
                GENOMEIDX1.to_string(),
                "-l".to_string(),
                loh_var.to_string_lossy().to_string(),
                tumor_bam.to_string_lossy().to_string(),
            ];
            let out = paths.readcount_dir.join(format!("{}.snp.LOH.hc.readcount", prefix));
            let child = spawn_command("bam-readcount", &args, Some(&out))?;
            push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
        }
    }

    wait_all(&mut children)?;
    log_message("bam-readcount complete");
    Ok(())
}

fn run_fpfilter(paths: &Paths) -> AppResult<()> {
    log_message("=== STAGE 9: False positive filtering (fpfilter.pl) ===");
    check_file_exists(Path::new(&format!("{}/fpfilter.pl", SCRIPTSDIR)))?;

    let mut children = Vec::new();
    for rc_file in list_files_matching(&paths.readcount_dir, "*.readcount")? {
        let base = rc_file
            .file_name()
            .and_then(OsStr::to_str)
            .unwrap_or_default()
            .trim_end_matches(".readcount")
            .to_string();
        let vcf = paths.somatic_dir.join(format!("{}.vcf", base));
        if !vcf.is_file() {
            log_message(&format!("WARNING: No matching VCF for {} — skipping fpfilter", base));
            continue;
        }

        let filtered = paths.filtered_dir.join(format!("{}.fpfilter.vcf", base));
        log_message(&format!("fpfilter: {}", base));
        let args = vec![
            format!("{}/fpfilter.pl", SCRIPTSDIR),
            "--vcf-var-file".to_string(),
            vcf.to_string_lossy().to_string(),
            "--readcount-file".to_string(),
            rc_file.to_string_lossy().to_string(),
            "--output-file".to_string(),
            filtered.to_string_lossy().to_string(),
            "--min-var-count".to_string(),
            "3".to_string(),
            "--min-var-freq".to_string(),
            MIN_VAR_FREQ.to_string(),
            "--min-ref-basequal".to_string(),
            BRC_BASE_QUAL.to_string(),
            "--min-var-basequal".to_string(),
            BRC_BASE_QUAL.to_string(),
        ];
        let child = spawn_command("perl", &args, None)?;
        push_child(&mut children, child, MAX_PARALLEL_JOBS)?;
    }

    wait_all(&mut children)?;
    log_message("False positive filtering complete");
    Ok(())
}

fn count_matches(dir: &Path, pattern: &str) -> AppResult<usize> {
    Ok(list_files_matching(dir, pattern)?.len())
}

fn generate_summary(paths: &Paths) -> AppResult<()> {
    log_message("=== STAGE 10: Generating summary report ===");

    let mut out = String::new();
    out.push_str("========================================================================\n");
    out.push_str("VarScan2 Pipeline Summary\n");
    out.push_str(&format!("Generated: {}\n", now_string()));
    out.push_str("========================================================================\n\n");
    out.push_str("CONFIGURATION\n");
    out.push_str(&format!("  Reference genome       : {}\n", GENOMEIDX1));
    out.push_str(&format!("  BAM directory          : {}\n", paths.bam_dir.display()));
    out.push_str(&format!("  BAM suffix             : {}\n", BAM_SUFFIX));
    out.push_str(&format!("  Pairs file             : {}\n", FILE_PAIRS_LIST));
    out.push_str(&format!("  Pairs suffix (stripped): {}\n", PAIRS_SUFFIX));
    out.push_str(&format!(
        "  Target regions (WES)   : {}\n",
        if TARGET_BED.is_empty() {
            "not set (full-genome pileup)"
        } else {
            TARGET_BED
        }
    ));
    out.push_str(&format!(
        "  VCF sample list        : {}\n\n",
        if VCF_SAMPLE_LIST.is_empty() {
            "not set (generic NORMAL/TUMOR headers)"
        } else {
            VCF_SAMPLE_LIST
        }
    ));

    out.push_str("RESULTS\n");
    out.push_str(&format!(
        "  Somatic SNP .hc   : {} files\n",
        count_matches(&paths.somatic_dir, "*.snp.Somatic.hc.vcf")?
    ));
    out.push_str(&format!(
        "  Somatic INDEL .hc : {} files\n",
        count_matches(&paths.somatic_dir, "*.indel.Somatic.hc.vcf")?
    ));
    out.push_str(&format!(
        "  Germline .hc      : {} files\n",
        count_matches(&paths.somatic_dir, "*.Germline.hc.vcf")?
    ));
    out.push_str(&format!(
        "  LOH .hc           : {} files\n",
        count_matches(&paths.somatic_dir, "*.LOH.hc.vcf")?
    ));
    out.push_str(&format!(
        "  Raw CNV           : {} files\n",
        count_matches(&paths.copy_number_dir, "*.copynumber")?
    ));
    out.push_str(&format!(
        "  Called CNV        : {} files\n",
        count_matches(&paths.copy_number_dir, "*.copynumber.called")?
    ));
    out.push_str(&format!(
        "  Homdel            : {} files\n",
        count_matches(&paths.copy_number_dir, "*.copynumber.homdel")?
    ));
    out.push_str(&format!(
        "  FP-filtered VCFs  : {} files\n",
        count_matches(&paths.filtered_dir, "*.fpfilter.vcf")?
    ));

    out.push_str("========================================================================\n");

    fs::write(&paths.summary_file, &out)
        .map_err(|e| format!("write {}: {}", paths.summary_file.display(), e))?;
    log_message(&format!("Summary written to: {}", paths.summary_file.display()));
    println!("{}", out);
    Ok(())
}

fn run_stage<F>(
    stage: u8,
    desc: &str,
    args: &Args,
    paths: &Paths,
    mut f: F,
) -> AppResult<()>
where
    F: FnMut() -> AppResult<()>,
{
    if stage < args.from_stage || stage > args.to_stage {
        log_message(&format!(
            "--- Stage {} skipped (outside range {}–{})",
            stage, args.from_stage, args.to_stage
        ));
        return Ok(());
    }

    if args.resume && resume_stage_match(paths, stage)? {
        log_message(&format!("--- Stage {} skipped — SHA256 resume match: {}", stage, desc));
        return Ok(());
    }

    if args.dry_run {
        log_message(&format!("--- Stage {} [dry-run] would run: {}", stage, desc));
        return Ok(());
    }

    f()?;
    write_stage_marker(paths, stage)?;
    Ok(())
}

fn parse_args() -> AppResult<Args> {
    let mut from_stage: u8 = 1;
    let mut to_stage: u8 = 10;
    let mut resume = false;
    let mut dry_run = false;

    let mut it = env::args().skip(1).peekable();
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "-h" | "--help" => {
                println!(
                    r#"Usage: varscan2_pipeline [--help] [--stage N] [--from N --to N] [--resume] [--dry-run]

VarScan2 paired tumour/normal somatic variant and CNV analysis pipeline.
Runs 10 sequential stages from BAM flagstats through false-positive filtering.

PIPELINE STAGES
  Stage 1   samtools flagstat        Alignment QC; read counts for CNV data-ratio
  Stage 2   samtools mpileup         Paired pileup (normal first, tumour second)
  Stage 3   VarScan somatic          Raw SNP and INDEL calls (.snp.vcf, .indel.vcf)
  Stage 4   VarScan processSomatic   Classify: Somatic / Germline / LOH (.hc.vcf)
  Stage 5   VarScan copynumber       Per-window log2 depth ratios (.copynumber)
  Stage 6   VarScan copyCaller       Segmented CNV calls (.called, .homdel)
  Stage 7   prepare .var files       Convert .hc.vcf to position lists (.var)
  Stage 8   bam-readcount            Per-position allele read support
  Stage 9   fpfilter.pl              False positive removal (.fpfilter.vcf)
  Stage 10  Summary report           varscan_pipeline_summary.txt

QUICK START
  1. Edit GENOMEIDX1 at the top of varscan2_pipeline.rs to point to your GRCh38
     reference FASTA, then rebuild: cargo build --release
  2. Populate sample_pairs.csv (tumour col 1, normal col 2):
       50T_CRC_final.bam,50N_CRC_final.bam
  3. Place VarScan.v2.3.9.jar in software/ and fpfilter.pl in scripts/
     (see software/README.md and scripts/README.md)
  4. Run: ./target/release/varscan2_pipeline

REQUIRED CONFIGURATION (compile-time constants in varscan2_pipeline.rs)
  GENOMEIDX1        Indexed reference FASTA (must match BAM chromosome naming)
  SOFTWAREDIR       Directory containing VarScan.v2.3.9.jar  (default: software)
  SCRIPTSDIR        Directory containing fpfilter.pl          (default: scripts)
  BAM_SUFFIX        BAM filename suffix                       (default: _final.bam)
  FILE_PAIRS_LIST   Tumour/normal pairs CSV                   (default: sample_pairs.csv)

OPTIONAL CONFIGURATION
  TARGET_BED        Capture kit BED file (WES only).
                    When set: restricts mpileup to captured regions, reduces
                    file size ~95%, excludes alt/decoy contigs automatically.
                    Leave empty for WGS.
  VCF_SAMPLE_LIST   Plain-text file with sample names (normal first, tumour
                    second, one per line) for correct VCF column labels.
                    Leave empty to use VarScan generic NORMAL/TUMOR headers.

KEY PARAMETERS (defaults shown — edit constants in varscan2_pipeline.rs)
  MIN_COVERAGE=20           Minimum site depth (both samples)
  MIN_COVERAGE_NORMAL=10    Minimum depth in normal
  MIN_COVERAGE_TUMOR=20     Minimum depth in tumour
  MIN_BASE_QUAL=20          Minimum base quality (mpileup + VarScan)
  MIN_VAR_FREQ=0.10         Minimum variant allele frequency to call
  TUMOR_PURITY=1.0          Tumour purity estimate — SET TO ACTUAL VALUE
                            VarScan uses this in its somatic Fisher's test.
                            Leaving at 1.0 for an impure tumour will
                            under-call subclonal variants.
  SOMATIC_P_VALUE=0.05      Significance threshold for somatic classification
  CNV_AMP_THRESHOLD=0.25    log2 ratio threshold for amplification call
  CNV_DEL_THRESHOLD=0.25    log2 ratio threshold for deletion call
  CNV_RECENTER_UP/DOWN=0    Baseline correction for chromosomally unstable
                            tumours — inspect *.copynumber before setting

PREREQUISITES
  Tool              Min version   Notes
  Rust + cargo      1.70          cargo build --release to compile
  Java JRE          8             Required for VarScan.v2.3.9.jar
  samtools          1.13          Must be on PATH
  bam-readcount     0.8           Must be on PATH
  Perl              5.10          Required for fpfilter.pl

  All BAM files must be coordinate-sorted, indexed (.bai present), and
  duplicate-marked (MarkDuplicates). Duplicates are excluded from pileup
  via -F 0x400. Chromosome naming in BAMs must match the reference FASTA
  (Ensembl: 1,2,...,X,Y,MT  vs  UCSC: chr1,chr2,...,chrX,chrY,chrM).

OUTPUT DIRECTORIES
  flagstats/        samtools flagstat outputs
  mpileup/          Paired .mpileup files
  somatic/          VarScan somatic and processSomatic outputs
  copynumber/       VarScan copynumber and copyCaller outputs
  snp-VAR/          SNP position lists (.var) for bam-readcount
  indel-VAR/        Somatic INDEL position lists (.var) for bam-readcount
  readcount/        bam-readcount outputs
  filtered/         fpfilter output VCFs (.fpfilter.vcf)
  .resume_state/    SHA256 stage markers for --resume validation

STAGE DEPENDENCIES (do not skip intermediate stages in a partial run)
  Stage 2 → requires Stage 1 flagstats (data-ratio for CNV)
  Stage 3 → requires Stage 2 mpileup
  Stage 4 → requires Stage 3 VCF output
  Stages 7–9 → require Stage 4 .hc.vcf output
  Stage 8 → requires Stage 7 .var files
  Stage 9 → requires Stage 8 .readcount files

OPTIONS
  --stage N       Run only stage N (equivalent to --from N --to N)
  --from N        Start execution at stage N (default: 1)
  --to N          Stop execution after stage N (default: 10)
  --resume        Skip a stage only when both are true:
                  1) output files are present, and
                  2) a stored SHA256 marker matches current inputs/params.
  --dry-run       Print which stages would execute without running them.
                  Prerequisite file checks are also skipped in dry-run mode.
  -h, --help      Show this help message and exit

EXAMPLES
  Run all stages:
    ./target/release/varscan2_pipeline

  Run a single stage:
    ./target/release/varscan2_pipeline --stage 5

  Run a range of stages:
    ./target/release/varscan2_pipeline --from 3 --to 6

  Resume after a failed or partial run:
    ./target/release/varscan2_pipeline --resume

  Resume from a specific stage:
    ./target/release/varscan2_pipeline --from 5 --resume

  Preview what would run without executing:
    ./target/release/varscan2_pipeline --dry-run

POST-PIPELINE STEPS
  1. Inspect *.copynumber files; adjust CNV_RECENTER_UP/DOWN if log2
     baseline is offset from zero (common in CIN tumours).
  2. Run circular binary segmentation (CBS) on *.copynumber.called.
  3. Annotate filtered VCFs with VEP or ANNOVAR.
  4. Validate high-impact somatic mutations.

See README.md for full parameter justifications and troubleshooting."#
                );
                std::process::exit(0);
            }
            "--stage" => {
                let v = it.next().ok_or_else(|| "--stage requires a value (1-10)".to_string())?;
                let n = v.parse::<u8>().map_err(|_| "--stage requires an integer".to_string())?;
                from_stage = n;
                to_stage = n;
            }
            "--from" => {
                let v = it.next().ok_or_else(|| "--from requires a value (1-10)".to_string())?;
                from_stage = v.parse::<u8>().map_err(|_| "--from requires an integer".to_string())?;
            }
            "--to" => {
                let v = it.next().ok_or_else(|| "--to requires a value (1-10)".to_string())?;
                to_stage = v.parse::<u8>().map_err(|_| "--to requires an integer".to_string())?;
            }
            "--resume" => resume = true,
            "--dry-run" => dry_run = true,
            _ => return Err(format!("Unknown argument: {}", arg)),
        }
    }

    if !(1..=10).contains(&from_stage) || !(1..=10).contains(&to_stage) {
        return Err("Stage numbers must be integers between 1 and 10".to_string());
    }
    if from_stage > to_stage {
        return Err(format!("--from ({}) must not exceed --to ({})", from_stage, to_stage));
    }

    Ok(Args {
        from_stage,
        to_stage,
        resume,
        dry_run,
    })
}

fn build_paths() -> AppResult<Paths> {
    let cwd = env::current_dir().map_err(|e| format!("current_dir: {}", e))?;
    Ok(Paths {        bam_dir: cwd.clone(),
        resume_state_dir: cwd.join(".resume_state"),
        flagstat_dir: cwd.join("flagstats"),
        mpileup_dir: cwd.join("mpileup"),
        somatic_dir: cwd.join("somatic"),
        copy_number_dir: cwd.join("copynumber"),
        snp_var_dir: cwd.join("snp-VAR"),
        indel_var_dir: cwd.join("indel-VAR"),
        readcount_dir: cwd.join("readcount"),
        filtered_dir: cwd.join("filtered"),
        summary_file: cwd.join("varscan_pipeline_summary.txt"),
    })
}

fn run() -> AppResult<()> {
    println!("========================================================================");
    println!("VarScan2 Somatic Variant and CNV Analysis Pipeline (Rust)");
    println!("========================================================================");

    let args = parse_args()?;
    let paths = build_paths()?;

    log_message("Starting VarScan2 Pipeline (Rust)");
    log_message(&format!(
        "Stages: {}–{}  |  Resume: {}  |  Dry-run: {}",
        args.from_stage, args.to_stage, args.resume, args.dry_run
    ));

    if !args.dry_run {
        for cmd in ["samtools", "java", "bam-readcount", "perl"] {
            check_command_exists(cmd)?;
        }

        check_file_exists(Path::new(FILE_PAIRS_LIST))?;
        if GENOMEIDX1.is_empty() {
            return Err(
                "GENOMEIDX1 constant is empty. Edit the constant at the top of \
                 varscan2_pipeline.rs to point to your reference FASTA, then \
                 rebuild with `cargo build --release`."
                    .to_string(),
            );
        }
        check_file_exists(Path::new(GENOMEIDX1))?;
        check_file_exists(Path::new(&format!("{}/VarScan.v2.3.9.jar", SOFTWAREDIR)))?;
        check_file_exists(Path::new(&format!("{}/fpfilter.pl", SCRIPTSDIR)))?;
    }

    setup_directories(&paths)?;

    run_stage(1, "samtools flagstat", &args, &paths, || generate_flagstats(&paths))?;
    run_stage(2, "samtools mpileup", &args, &paths, || generate_mpileup(&paths))?;
    run_stage(3, "VarScan somatic", &args, &paths, || run_varscan_somatic(&paths))?;
    run_stage(4, "VarScan processSomatic", &args, &paths, || process_somatic_variants(&paths))?;

    if args.to_stage >= 7 && !args.dry_run {
        check_stage4_output(&paths)?;
    }

    run_stage(5, "VarScan copynumber", &args, &paths, || run_varscan_copynumber(&paths))?;
    run_stage(6, "VarScan copyCaller", &args, &paths, || run_copy_caller(&paths))?;
    run_stage(7, "prepare .var files", &args, &paths, || prepare_filter_input(&paths))?;
    run_stage(8, "bam-readcount", &args, &paths, || run_bam_readcount(&paths))?;
    run_stage(9, "fpfilter.pl", &args, &paths, || run_fpfilter(&paths))?;
    run_stage(10, "summary", &args, &paths, || generate_summary(&paths))?;

    if !args.dry_run {
        log_message("VarScan2 Pipeline complete");
    }
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        log_message(&format!("ERROR: {}", e));
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write as IoWrite;

    // ── glob_match ──────────────────────────────────────────────────────────

    #[test]
    fn glob_star_matches_everything() {
        assert!(glob_match("*", "anything.bam"));
        assert!(glob_match("*", ""));
    }

    #[test]
    fn glob_no_wildcard_exact() {
        assert!(glob_match("file.bam", "file.bam"));
        assert!(!glob_match("file.bam", "FILE.bam"));
        assert!(!glob_match("file.bam", "file.bam.bai"));
    }

    #[test]
    fn glob_suffix_wildcard() {
        assert!(glob_match("*.flagstats", "sample.flagstats"));
        assert!(glob_match("*.flagstats", "a.b.flagstats"));
        assert!(!glob_match("*.flagstats", "sample.flagstats.gz"));
        assert!(!glob_match("*.flagstats", "sample.vcf"));
    }

    #[test]
    fn glob_prefix_wildcard() {
        assert!(glob_match("prefix.*", "prefix.bam"));
        assert!(glob_match("prefix.*", "prefix.flagstats"));
        assert!(!glob_match("prefix.*", "xprefix.bam"));
    }

    #[test]
    fn glob_middle_wildcard() {
        assert!(glob_match("sample*final.bam", "sample_CRC_final.bam"));
        assert!(glob_match("sample*final.bam", "samplefinal.bam"));
        assert!(!glob_match("sample*final.bam", "sample_CRC_final.bam.bai"));
    }

    #[test]
    fn glob_multi_wildcard() {
        assert!(glob_match("a*b*c", "aXbYc"));
        assert!(glob_match("a*b*c", "abbc"));
        assert!(glob_match("a*b*c", "abc"));
        assert!(!glob_match("a*b*c", "aXbY"));
        assert!(!glob_match("a*b*c", "XbYc"));
    }

    // ── clean_pair_field ────────────────────────────────────────────────────

    #[test]
    fn clean_pair_strips_spaces_and_cr() {
        assert_eq!(clean_pair_field("  sample.bam  "), "sample.bam");
        assert_eq!(clean_pair_field("sample.bam\r"), "sample.bam");
        assert_eq!(clean_pair_field(" sample.bam\r"), "sample.bam");
        assert_eq!(clean_pair_field("sample.bam"), "sample.bam");
    }

    // ── get_sample_name ─────────────────────────────────────────────────────

    #[test]
    fn get_sample_name_strips_suffix() {
        assert_eq!(get_sample_name("50T_CRC_final.bam"), "50T_CRC");
        assert_eq!(get_sample_name("sample_final.bam"), "sample");
    }

    #[test]
    fn get_sample_name_no_suffix_unchanged() {
        assert_eq!(get_sample_name("sample.bam"), "sample.bam");
        assert_eq!(get_sample_name("sample"), "sample");
    }

    // ── parse_primary_mapped_count ──────────────────────────────────────────

    fn write_tempfile(content: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "varscan_test_{}.flagstats",
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .subsec_nanos()
        ));
        let mut f = std::fs::File::create(&path).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        path
    }

    #[test]
    fn parse_flagstats_samtools_1_13_format() {
        let content = "\
37500000 + 0 in total (QC-passed reads + QC-failed reads)\n\
37500000 + 0 primary\n\
0 + 0 secondary\n\
0 + 0 supplementary\n\
1234567 + 0 duplicates\n\
1234567 + 0 primary duplicates\n\
36800000 + 0 mapped (98.13% : N/A)\n\
36800000 + 0 primary mapped (98.13% : N/A)\n\
";
        let path = write_tempfile(content);
        let result = parse_primary_mapped_count(&path).unwrap();
        assert_eq!(result, Some(36800000u64));
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn parse_flagstats_missing_primary_mapped_line() {
        let content = "\
37500000 + 0 in total (QC-passed reads + QC-failed reads)\n\
37500000 + 0 mapped (98.13% : N/A)\n\
";
        let path = write_tempfile(content);
        let result = parse_primary_mapped_count(&path).unwrap();
        assert_eq!(result, None);
        std::fs::remove_file(path).ok();
    }

    // ── sha256_of_string ────────────────────────────────────────────────────

    #[test]
    fn sha256_of_empty_string() {
        let h = sha256_of_string("").unwrap();
        assert_eq!(
            h,
            "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
        );
    }

    #[test]
    fn sha256_of_known_string() {
        let h = sha256_of_string("hello").unwrap();
        assert_eq!(
            h,
            "2cf24dba5fb0a30e26e83b2ac5b9e29e1b161e5c1fa7425e73043362938b9824"
        );
    }
}
