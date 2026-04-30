use sha2::{Digest, Sha256};
use serde::Deserialize;
use std::env;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::time::UNIX_EPOCH;

// ── Config structs ────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct PathsConfig {
    pub reference:       String,
    pub bam_dir:         String,
    pub bam_suffix:      String,
    pub pairs_file:      String,
    pub pairs_suffix:    String,
    pub target_bed:      String,
    pub vcf_sample_list: String,
    pub software_dir:    String,
    pub scripts_dir:     String,
}

impl Default for PathsConfig {
    fn default() -> Self {
        Self {
            reference:       String::new(),
            bam_dir:         ".".to_string(),
            bam_suffix:      "_final.bam".to_string(),
            pairs_file:      "sample_pairs.csv".to_string(),
            pairs_suffix:    "_final.bam".to_string(),
            target_bed:      String::new(),
            vcf_sample_list: String::new(),
            software_dir:    "software".to_string(),
            scripts_dir:     "scripts".to_string(),
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct SomaticConfig {
    pub min_coverage:        i32,
    pub min_coverage_normal: i32,
    pub min_coverage_tumor:  i32,
    pub min_base_qual:       i32,
    pub min_var_freq:        f64,
    pub min_freq_for_hom:    f64,
    pub normal_purity:       f64,
    pub tumor_purity:        f64,
    pub p_value:             f64,
    pub somatic_p_value:     f64,
    pub strand_filter:       i32,
}

impl Default for SomaticConfig {
    fn default() -> Self {
        Self {
            min_coverage:        20,
            min_coverage_normal: 10,
            min_coverage_tumor:  20,
            min_base_qual:       20,
            min_var_freq:        0.10,
            min_freq_for_hom:    0.75,
            normal_purity:       1.0,
            tumor_purity:        1.0,
            p_value:             0.99,
            somatic_p_value:     0.05,
            strand_filter:       1,
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct ProcessSomaticConfig {
    pub min_tumor_freq:  f64,
    pub max_normal_freq: f64,
    pub p_value:         f64,
}

impl Default for ProcessSomaticConfig {
    fn default() -> Self {
        Self {
            min_tumor_freq:  0.10,
            max_normal_freq: 0.05,
            p_value:         0.05,
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct CnvConfig {
    pub min_coverage:     i32,
    pub p_value:          f64,
    pub min_segment_size: i32,
    pub max_segment_size: i32,
    pub amp_threshold:    f64,
    pub del_threshold:    f64,
    pub recenter_up:      f64,
    pub recenter_down:    f64,
}

impl Default for CnvConfig {
    fn default() -> Self {
        Self {
            min_coverage:     20,
            p_value:          0.01,
            min_segment_size: 10,
            max_segment_size: 100,
            amp_threshold:    0.25,
            del_threshold:    0.25,
            recenter_up:      0.0,
            recenter_down:    0.0,
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct ReadcountConfig {
    pub map_qual:          i32,
    pub base_qual:         i32,
    pub max_parallel_jobs: usize,
}

impl Default for ReadcountConfig {
    fn default() -> Self {
        Self {
            map_qual:          10,
            base_qual:         15,
            max_parallel_jobs: 30,
        }
    }
}

#[derive(Debug, Clone, Default, Deserialize)]
#[serde(default)]
pub struct Config {
    pub paths:           PathsConfig,
    pub somatic:         SomaticConfig,
    pub process_somatic: ProcessSomaticConfig,
    pub cnv:             CnvConfig,
    pub readcount:       ReadcountConfig,
}


type AppResult<T> = Result<T, String>;

fn load_config(path: Option<&str>) -> AppResult<Config> {
    let toml_path = path.unwrap_or("config.toml");
    if !std::path::Path::new(toml_path).is_file() {
        if path.is_some() {
            return Err(format!("Config file not found: {}", toml_path));
        }
        return Ok(Config::default());
    }
    let content = fs::read_to_string(toml_path)
        .map_err(|e| format!("read config {}: {}", toml_path, e))?;
    toml::from_str::<Config>(&content)
        .map_err(|e| format!("parse config {}: {}", toml_path, e))
}

fn apply_env_overrides(cfg: &mut Config) -> AppResult<()> {
    macro_rules! env_str {
        ($var:expr, $field:expr) => {
            if let Ok(v) = std::env::var($var) { $field = v; }
        };
    }
    macro_rules! env_int {
        ($var:expr, $field:expr) => {
            if let Ok(v) = std::env::var($var) {
                $field = v.parse().map_err(|_| {
                    format!("{} must be an integer, got: {}", $var, v)
                })?;
            }
        };
    }
    macro_rules! env_float {
        ($var:expr, $field:expr) => {
            if let Ok(v) = std::env::var($var) {
                $field = v.parse().map_err(|_| {
                    format!("{} must be a float, got: {}", $var, v)
                })?;
            }
        };
    }

    env_str!  ("VARSCAN_REFERENCE",           cfg.paths.reference);
    env_str!  ("VARSCAN_BAM_DIR",             cfg.paths.bam_dir);
    env_str!  ("VARSCAN_BAM_SUFFIX",          cfg.paths.bam_suffix);
    env_str!  ("VARSCAN_PAIRS_FILE",          cfg.paths.pairs_file);
    env_str!  ("VARSCAN_PAIRS_SUFFIX",        cfg.paths.pairs_suffix);
    env_str!  ("VARSCAN_TARGET_BED",          cfg.paths.target_bed);
    env_str!  ("VARSCAN_VCF_SAMPLE_LIST",     cfg.paths.vcf_sample_list);
    env_str!  ("VARSCAN_SOFTWARE_DIR",        cfg.paths.software_dir);
    env_str!  ("VARSCAN_SCRIPTS_DIR",         cfg.paths.scripts_dir);
    env_int!  ("VARSCAN_MIN_COVERAGE",        cfg.somatic.min_coverage);
    env_int!  ("VARSCAN_MIN_COVERAGE_NORMAL", cfg.somatic.min_coverage_normal);
    env_int!  ("VARSCAN_MIN_COVERAGE_TUMOR",  cfg.somatic.min_coverage_tumor);
    env_int!  ("VARSCAN_MIN_BASE_QUAL",       cfg.somatic.min_base_qual);
    env_float!("VARSCAN_MIN_VAR_FREQ",        cfg.somatic.min_var_freq);
    env_float!("VARSCAN_MIN_FREQ_FOR_HOM",    cfg.somatic.min_freq_for_hom);
    env_float!("VARSCAN_NORMAL_PURITY",       cfg.somatic.normal_purity);
    env_float!("VARSCAN_TUMOR_PURITY",        cfg.somatic.tumor_purity);
    env_float!("VARSCAN_P_VALUE",             cfg.somatic.p_value);
    env_float!("VARSCAN_SOMATIC_P_VALUE",     cfg.somatic.somatic_p_value);
    env_int!  ("VARSCAN_STRAND_FILTER",       cfg.somatic.strand_filter);
    env_float!("VARSCAN_MIN_TUMOR_FREQ",      cfg.process_somatic.min_tumor_freq);
    env_float!("VARSCAN_MAX_NORMAL_FREQ",     cfg.process_somatic.max_normal_freq);
    env_float!("VARSCAN_PROCESS_P_VALUE",     cfg.process_somatic.p_value);
    env_int!  ("VARSCAN_CNV_MIN_COVERAGE",    cfg.cnv.min_coverage);
    env_float!("VARSCAN_CNV_P_VALUE",         cfg.cnv.p_value);
    env_int!  ("VARSCAN_MIN_SEGMENT_SIZE",    cfg.cnv.min_segment_size);
    env_int!  ("VARSCAN_MAX_SEGMENT_SIZE",    cfg.cnv.max_segment_size);
    env_float!("VARSCAN_CNV_AMP_THRESHOLD",   cfg.cnv.amp_threshold);
    env_float!("VARSCAN_CNV_DEL_THRESHOLD",   cfg.cnv.del_threshold);
    env_float!("VARSCAN_CNV_RECENTER_UP",     cfg.cnv.recenter_up);
    env_float!("VARSCAN_CNV_RECENTER_DOWN",   cfg.cnv.recenter_down);
    env_int!  ("VARSCAN_BRC_MAP_QUAL",        cfg.readcount.map_qual);
    env_int!  ("VARSCAN_BRC_BASE_QUAL",       cfg.readcount.base_qual);

    Ok(())
}

fn config_example_toml() -> String {
    r#"# VarScan2 Pipeline Configuration
# Generated by: varscan2_pipeline --init-config
# Edit 'reference' and optionally 'bam_dir', then run:
#   varscan2_pipeline --validate
#   varscan2_pipeline --resume

[paths]
reference    = ""           # REQUIRED: absolute path to indexed GRCh38 FASTA
bam_dir      = "."          # directory containing BAM files (default: cwd)
bam_suffix   = "_final.bam"
pairs_file   = "sample_pairs.csv"
pairs_suffix = "_final.bam"
target_bed   = ""           # WES: path to capture BED; leave empty for WGS
software_dir = "software"   # contains VarScan.v2.3.9.jar
scripts_dir  = "scripts"    # contains fpfilter.pl

[somatic]
min_coverage        = 20
min_coverage_normal = 10
min_coverage_tumor  = 20
min_base_qual       = 20
min_var_freq        = 0.10
min_freq_for_hom    = 0.75
normal_purity       = 1.0
tumor_purity        = 1.0   # override per-sample in pairs CSV column 3
p_value             = 0.99
somatic_p_value     = 0.05
strand_filter       = 1

[process_somatic]
min_tumor_freq  = 0.10
max_normal_freq = 0.05
p_value         = 0.05

[cnv]
min_coverage     = 20
p_value          = 0.01
min_segment_size = 10
max_segment_size = 100
amp_threshold    = 0.25
del_threshold    = 0.25
recenter_up      = 0.0
recenter_down    = 0.0

[readcount]
map_qual          = 10
base_qual         = 15
max_parallel_jobs = 30
"#.to_string()
}

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

fn validate_config(cfg: &Config) -> bool {
    let mut ok = true;

    macro_rules! check_ok {
        ($cond:expr, $msg:expr) => {
            if $cond { println!("[OK]   {}", $msg); }
            else      { println!("[FAIL] {}", $msg); ok = false; }
        };
    }
    macro_rules! warn {
        ($msg:expr) => { println!("[WARN] {}", $msg); };
    }

    // reference
    if cfg.paths.reference.is_empty() {
        println!("[FAIL] reference: not set — add [paths] reference = \"...\" to config.toml");
        ok = false;
    } else {
        let ref_path = Path::new(&cfg.paths.reference);
        check_ok!(ref_path.is_file(), format!("reference: {} (file present)", cfg.paths.reference));
        if ref_path.is_file() {
            let fai = format!("{}.fai", cfg.paths.reference);
            check_ok!(Path::new(&fai).is_file(), format!("reference index: {}", fai));
        }
    }

    // pairs file
    let pairs_ok = Path::new(&cfg.paths.pairs_file).is_file();
    check_ok!(pairs_ok, format!("pairs file: {}", cfg.paths.pairs_file));
    let pair_count = if pairs_ok {
        match std::fs::File::open(&cfg.paths.pairs_file)
            .map(std::io::BufReader::new)
        {
            Ok(r) => {
                use std::io::BufRead;
                let n = r.lines().filter(|l| {
                    l.as_ref().map(|s| !s.trim().is_empty()).unwrap_or(false)
                }).count();
                println!("[OK]   pairs count: {} pair(s)", n);
                n
            }
            Err(_) => 0,
        }
    } else { 0 };

    // bam files from pairs CSV
    if pairs_ok {
        let cfg_clone = cfg.clone();
        if let Ok(pairs) = read_pairs(&cfg_clone) {
            for (entry1, entry2, _) in &pairs {
                let bam_dir_pb = if cfg.paths.bam_dir.is_empty() || cfg.paths.bam_dir == "." {
                    std::env::current_dir().unwrap_or_default()
                } else {
                    PathBuf::from(&cfg.paths.bam_dir)
                };
                let fake_paths = Paths {
                    bam_dir: bam_dir_pb,
                    resume_state_dir: PathBuf::new(),
                    flagstat_dir: PathBuf::new(),
                    mpileup_dir: PathBuf::new(),
                    somatic_dir: PathBuf::new(),
                    copy_number_dir: PathBuf::new(),
                    snp_var_dir: PathBuf::new(),
                    indel_var_dir: PathBuf::new(),
                    readcount_dir: PathBuf::new(),
                    filtered_dir: PathBuf::new(),
                    summary_file: PathBuf::new(),
                };
                for entry in [entry1, entry2] {
                    let bam = get_bam_path(&fake_paths, entry, cfg);
                    if bam.is_file() {
                        match check_bam_index(&bam) {
                            Ok(_)  => println!("[OK]   {}: indexed", bam.display()),
                            Err(e) => { println!("[FAIL] {}", e); ok = false; }
                        }
                    } else {
                        println!("[FAIL] {}: not found", bam.display());
                        ok = false;
                    }
                }
            }
        }
    }

    // tools on PATH
    for cmd in ["samtools", "java", "bam-readcount", "perl"] {
        check_ok!(check_command_exists(cmd).is_ok(), format!("{} on PATH", cmd));
    }

    // jars / scripts
    let jar = format!("{}/VarScan.v2.3.9.jar", cfg.paths.software_dir);
    check_ok!(Path::new(&jar).is_file(), format!("VarScan.v2.3.9.jar: {}", jar));
    let fp = format!("{}/fpfilter.pl", cfg.paths.scripts_dir);
    check_ok!(Path::new(&fp).is_file(), format!("fpfilter.pl: {}", fp));

    // warnings
    if cfg.paths.target_bed.is_empty() {
        warn!("target_bed not set — full-genome mpileup (WGS mode)");
    }
    if cfg.somatic.tumor_purity == 1.0 && pair_count > 0 {
        warn!("tumor_purity=1.0 (default) for all pairs — set col 3 in pairs CSV if purity is known");
    }

    ok
}

#[derive(Clone)]
struct Args {
    from_stage:  u8,
    to_stage:    u8,
    resume:      bool,
    dry_run:     bool,
    config_path: Option<String>,
    init_config: bool,
    validate:    bool,
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

fn get_sample_name(entry: &str, pairs_suffix: &str) -> String {
    if entry.ends_with(pairs_suffix) {
        entry.trim_end_matches(pairs_suffix).to_string()
    } else {
        entry.to_string()
    }
}

fn get_bam_path(paths: &Paths, entry: &str, cfg: &Config) -> PathBuf {
    let mut p = paths.bam_dir.clone();
    p.push(format!(
        "{}{}",
        get_sample_name(entry, &cfg.paths.pairs_suffix),
        cfg.paths.bam_suffix
    ));
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

fn read_pairs(cfg: &Config) -> AppResult<Vec<(String, String, Option<f64>)>> {
    let file = File::open(&cfg.paths.pairs_file)
        .map_err(|e| format!("open {}: {}", cfg.paths.pairs_file, e))?;
    let reader = BufReader::new(file);
    let mut pairs = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("read {}: {}", cfg.paths.pairs_file, e))?;
        let mut parts = line.splitn(3, ',');
        let p1 = clean_pair_field(parts.next().unwrap_or(""));
        let p2 = clean_pair_field(parts.next().unwrap_or(""));
        if p1.is_empty() || p2.is_empty() {
            log_message(&format!("WARNING: Incomplete pair (empty field) — skipping: '{}','{}'", p1, p2));
            continue;
        }
        let purity = parts.next()
            .map(clean_pair_field)
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<f64>().map_err(|_| {
                format!("invalid tumor_purity '{}' for pair {},{} in {}", s, p1, p2, cfg.paths.pairs_file)
            }))
            .transpose()?;
        pairs.push((p1, p2, purity));
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

fn hash_pairs_cleaned(cfg: &Config) -> AppResult<String> {
    if !Path::new(&cfg.paths.pairs_file).is_file() {
        return Ok("MISSING".to_string());
    }
    let pairs = read_pairs(cfg)?;
    if pairs.is_empty() {
        return Ok("EMPTY".to_string());
    }
    let normalized = pairs
        .iter()
        .map(|(t, n, p)| match p {
            Some(v) => format!("{},{},{}", t, n, v),
            None    => format!("{},{}", t, n),
        })
        .collect::<Vec<_>>()
        .join("\n");
    sha256_of_string(&normalized)
}

fn hash_paired_bam_metadata(paths: &Paths, cfg: &Config) -> AppResult<String> {
    if !Path::new(&cfg.paths.pairs_file).is_file() {
        return Ok("MISSING".to_string());
    }

    let pairs = read_pairs(cfg)?;
    if pairs.is_empty() {
        return Ok("EMPTY".to_string());
    }

    let mut lines = Vec::new();
    for (entry1, entry2, _) in pairs {
        let normal_bam = get_bam_path(paths, &entry2, cfg);
        let tumor_bam  = get_bam_path(paths, &entry1, cfg);
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

fn compute_stage_hash(stage: u8, paths: &Paths, cfg: &Config) -> AppResult<String> {
    let exe = env::current_exe().map_err(|e| format!("current_exe: {}", e))?;
    let mut payload = String::new();
    payload.push_str(&format!("exe_sha256={}\n", sha256_of_file_or_missing(&exe)?));
    payload.push_str(&format!("pairs_sha256={}\n", hash_pairs_cleaned(cfg)?));
    payload.push_str(&format!("bam_sha256={}\n",   hash_paired_bam_metadata(paths, cfg)?));
    payload.push_str(&format!("stage={}\n", stage));

    let jar = format!("{}/VarScan.v2.3.9.jar", cfg.paths.software_dir);

    match stage {
        1 => {
            payload.push_str(&format!("bam_dir={}\n",    paths.bam_dir.display()));
            payload.push_str(&format!("bam_suffix={}\n", cfg.paths.bam_suffix));
        }
        2 => {
            payload.push_str(&format!("reference_sha256={}\n",  file_sha(&cfg.paths.reference)?));
            payload.push_str(&format!("target_bed_sha256={}\n", file_sha(&cfg.paths.target_bed)?));
            payload.push_str(&format!("BRC_MAP_QUAL={}\n",      cfg.readcount.map_qual));
            payload.push_str(&format!("MIN_BASE_QUAL={}\n",     cfg.somatic.min_base_qual));
        }
        3 => {
            payload.push_str(&format!("jar_sha256={}\n",       file_sha(&jar)?));
            payload.push_str(&format!("mpileup_manifest={}\n", hash_manifest(&paths.mpileup_dir, "*.mpileup")?));
            payload.push_str(&format!("vcf_sample_sha={}\n",   file_sha(&cfg.paths.vcf_sample_list)?));
            payload.push_str(&format!("MIN_COVERAGE={}\n",        cfg.somatic.min_coverage));
            payload.push_str(&format!("MIN_COVERAGE_NORMAL={}\n", cfg.somatic.min_coverage_normal));
            payload.push_str(&format!("MIN_COVERAGE_TUMOR={}\n",  cfg.somatic.min_coverage_tumor));
            payload.push_str(&format!("MIN_VAR_FREQ={}\n",        cfg.somatic.min_var_freq));
            payload.push_str(&format!("MIN_FREQ_FOR_HOM={}\n",    cfg.somatic.min_freq_for_hom));
            payload.push_str(&format!("NORMAL_PURITY={}\n",       cfg.somatic.normal_purity));
            payload.push_str(&format!("TUMOR_PURITY={}\n",        cfg.somatic.tumor_purity));
            payload.push_str(&format!("P_VALUE={}\n",             cfg.somatic.p_value));
            payload.push_str(&format!("SOMATIC_P_VALUE={}\n",     cfg.somatic.somatic_p_value));
            payload.push_str(&format!("STRAND_FILTER={}\n",       cfg.somatic.strand_filter));
            payload.push_str(&format!("BRC_MAP_QUAL={}\n",        cfg.readcount.map_qual));
            payload.push_str(&format!("MIN_BASE_QUAL={}\n",       cfg.somatic.min_base_qual));
        }
        4 => {
            payload.push_str(&format!("jar_sha256={}\n",     file_sha(&jar)?));
            payload.push_str(&format!("snp_manifest={}\n",   hash_manifest(&paths.somatic_dir, "*.snp.vcf")?));
            payload.push_str(&format!("indel_manifest={}\n", hash_manifest(&paths.somatic_dir, "*.indel.vcf")?));
            payload.push_str(&format!("MIN_TUMOR_FREQ={}\n",   cfg.process_somatic.min_tumor_freq));
            payload.push_str(&format!("MAX_NORMAL_FREQ={}\n",  cfg.process_somatic.max_normal_freq));
            payload.push_str(&format!("PROCESS_P_VALUE={}\n",  cfg.process_somatic.p_value));
        }
        5 => {
            payload.push_str(&format!("jar_sha256={}\n",          file_sha(&jar)?));
            payload.push_str(&format!("mpileup_manifest={}\n",    hash_manifest(&paths.mpileup_dir, "*.mpileup")?));
            payload.push_str(&format!("flagstats_manifest={}\n",  hash_manifest(&paths.flagstat_dir, "*.flagstats")?));
            payload.push_str(&format!("CNV_MIN_COVERAGE={}\n",    cfg.cnv.min_coverage));
            payload.push_str(&format!("MIN_BASE_QUAL={}\n",       cfg.somatic.min_base_qual));
            payload.push_str(&format!("BRC_MAP_QUAL={}\n",        cfg.readcount.map_qual));
            payload.push_str(&format!("MIN_SEGMENT_SIZE={}\n",    cfg.cnv.min_segment_size));
            payload.push_str(&format!("MAX_SEGMENT_SIZE={}\n",    cfg.cnv.max_segment_size));
            payload.push_str(&format!("CNV_P_VALUE={}\n",         cfg.cnv.p_value));
        }
        6 => {
            payload.push_str(&format!("jar_sha256={}\n",      file_sha(&jar)?));
            payload.push_str(&format!("cnv_manifest={}\n",    hash_manifest(&paths.copy_number_dir, "*.copynumber")?));
            payload.push_str(&format!("CNV_AMP_THRESHOLD={}\n",  cfg.cnv.amp_threshold));
            payload.push_str(&format!("CNV_DEL_THRESHOLD={}\n",  cfg.cnv.del_threshold));
            payload.push_str(&format!("CNV_RECENTER_UP={}\n",    cfg.cnv.recenter_up));
            payload.push_str(&format!("CNV_RECENTER_DOWN={}\n",  cfg.cnv.recenter_down));
        }
        7 => {
            payload.push_str(&format!("hc_manifest={}\n", hash_manifest(&paths.somatic_dir, "*.hc.vcf")?));
        }
        8 => {
            payload.push_str(&format!("reference_sha256={}\n",   file_sha(&cfg.paths.reference)?));
            payload.push_str(&format!("snp_var_manifest={}\n",   hash_manifest(&paths.snp_var_dir, "*.var")?));
            payload.push_str(&format!("indel_var_manifest={}\n", hash_manifest(&paths.indel_var_dir, "*.var")?));
            payload.push_str(&format!("BRC_MAP_QUAL={}\n",  cfg.readcount.map_qual));
            payload.push_str(&format!("BRC_BASE_QUAL={}\n", cfg.readcount.base_qual));
        }
        9 => {
            payload.push_str(&format!("fpfilter_sha256={}\n",    file_sha(&format!("{}/fpfilter.pl", cfg.paths.scripts_dir))?));
            payload.push_str(&format!("readcount_manifest={}\n", hash_manifest(&paths.readcount_dir, "*.readcount")?));
            payload.push_str(&format!("somatic_manifest={}\n",   hash_manifest(&paths.somatic_dir, "*.vcf")?));
            payload.push_str(&format!("MIN_VAR_FREQ={}\n",  cfg.somatic.min_var_freq));
            payload.push_str(&format!("BRC_BASE_QUAL={}\n", cfg.readcount.base_qual));
        }
        10 => {
            payload.push_str(&format!("somatic_manifest={}\n",  hash_manifest(&paths.somatic_dir, "*.vcf")?));
            payload.push_str(&format!("copy_manifest={}\n",     hash_manifest(&paths.copy_number_dir, "*.copynumber*")?));
            payload.push_str(&format!("filtered_manifest={}\n", hash_manifest(&paths.filtered_dir, "*.fpfilter.vcf")?));
        }
        _ => {}
    }

    sha256_of_string(&payload)
}

fn stage_marker_path(paths: &Paths, stage: u8) -> PathBuf {
    paths.resume_state_dir.join(format!("stage_{}.sha256", stage))
}

fn write_stage_marker(paths: &Paths, stage: u8, cfg: &Config) -> AppResult<()> {
    let marker = stage_marker_path(paths, stage);
    let hash = compute_stage_hash(stage, paths, cfg)?;
    fs::write(&marker, format!("{}\n", hash))
        .map_err(|e| format!("write marker {}: {}", marker.display(), e))
}

fn resume_stage_match(paths: &Paths, stage: u8, cfg: &Config) -> AppResult<bool> {
    let marker = stage_marker_path(paths, stage);
    if !marker.is_file() {
        return Ok(false);
    }
    let expected = fs::read_to_string(&marker)
        .map_err(|e| format!("read marker {}: {}", marker.display(), e))?
        .trim()
        .to_string();
    let current = compute_stage_hash(stage, paths, cfg)?;
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

fn generate_flagstats(paths: &Paths, cfg: &Config) -> AppResult<()> {
    log_message("=== STAGE 1: Generating BAM flagstats ===");
    if !paths.bam_dir.is_dir() {
        return Err(format!("BAM directory not found: {}", paths.bam_dir.display()));
    }

    let mut children = Vec::new();
    let mut bams = list_files_matching(&paths.bam_dir, &format!("*{}", cfg.paths.bam_suffix))?;
    bams.sort();
    for bam in bams {
        let sample = bam
            .file_name()
            .and_then(OsStr::to_str)
            .unwrap_or_default()
            .trim_end_matches(cfg.paths.bam_suffix.as_str())
            .to_string();
        log_message(&format!("Flagstats: {}", sample));
        let out = paths.flagstat_dir.join(format!("{}.flagstats", sample));
        let args = vec!["flagstat".to_string(), bam.to_string_lossy().to_string()];
        let child = spawn_command("samtools", &args, Some(&out))?;
        push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
    }
    wait_all(&mut children)?;
    log_message("Flagstats complete");
    Ok(())
}

fn generate_mpileup(paths: &Paths, cfg: &Config) -> AppResult<()> {
    log_message("=== STAGE 2: Generating paired mpileup (normal first, tumour second) ===");
    check_file_exists(Path::new(&cfg.paths.pairs_file))?;
    check_file_exists(Path::new(&cfg.paths.reference))?;
    if !cfg.paths.target_bed.is_empty() {
        check_file_exists(Path::new(&cfg.paths.target_bed))?;
    }

    let pairs = read_pairs(cfg)?;
    let mut children = Vec::new();

    for (entry1, entry2, _) in pairs {
        let samplet = get_sample_name(&entry1, &cfg.paths.pairs_suffix);
        let samplen = get_sample_name(&entry2, &cfg.paths.pairs_suffix);
        let tumor_bam = get_bam_path(paths, &entry1, cfg);
        let normal_bam = get_bam_path(paths, &entry2, cfg);
        check_file_exists(&normal_bam)?;
        check_file_exists(&tumor_bam)?;

        let out = paths.mpileup_dir.join(format!("{}_{}.mpileup", samplen, samplet));
        log_message(&format!("Mpileup: Normal={}  Tumour={}", samplen, samplet));

        // -B: disables BAQ to avoid over-penalising reads at WES capture boundaries;
        //     increases INDEL FP rate slightly — remove for strict WGS variant calling.
        // -q map_qual (10): more conservative than VarScan manual's -q 1; excludes
        //     low-confidence multi-mapped reads (MAPQ 1-9). Intentional.
        // -Q min_base_qual: pre-filter bases before pileup; VarScan --min-base-qual
        //     below is then redundant but left for explicit documentation of intent.
        let mut args = vec![
            "mpileup".to_string(),
            "-B".to_string(),
            "-q".to_string(),
            cfg.readcount.map_qual.to_string(),
            "-Q".to_string(),
            cfg.somatic.min_base_qual.to_string(),
            "-F".to_string(),
            "0x400".to_string(),
        ];
        if !cfg.paths.target_bed.is_empty() {
            args.push("-l".to_string());
            args.push(cfg.paths.target_bed.clone());
        }
        args.push("-f".to_string());
        args.push(cfg.paths.reference.clone());
        args.push(normal_bam.to_string_lossy().to_string());
        args.push(tumor_bam.to_string_lossy().to_string());

        let child = spawn_command("samtools", &args, Some(&out))?;
        push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
    }

    wait_all(&mut children)?;
    log_message("Mpileup complete");
    Ok(())
}

fn run_varscan_somatic(paths: &Paths, cfg: &Config) -> AppResult<()> {
    log_message("=== STAGE 3: VarScan somatic variant calling ===");
    check_file_exists(Path::new(&cfg.paths.pairs_file))?;
    check_file_exists(Path::new(&format!("{}/VarScan.v2.3.9.jar", cfg.paths.software_dir)))?;

    let pairs = read_pairs(cfg)?;
    let mut children = Vec::new();

    for (entry1, entry2, sample_purity) in pairs {
        let samplet = get_sample_name(&entry1, &cfg.paths.pairs_suffix);
        let samplen = get_sample_name(&entry2, &cfg.paths.pairs_suffix);
        let mpileup_file = paths.mpileup_dir.join(format!("{}_{}.mpileup", samplen, samplet));
        if !mpileup_file.is_file() {
            log_message(&format!("WARNING: Mpileup not found: {} — skipping pair", mpileup_file.display()));
            continue;
        }

        let effective_purity = sample_purity.unwrap_or(cfg.somatic.tumor_purity);
        log_message(&format!("VarScan somatic: Normal={}  Tumour={}  tumor_purity={}", samplen, samplet, effective_purity));
        let mut args = vec![
            "-jar".to_string(),
            format!("{}/VarScan.v2.3.9.jar", cfg.paths.software_dir),
            "somatic".to_string(),
            mpileup_file.to_string_lossy().to_string(),
            paths.somatic_dir.join(format!("{}_{}", samplen, samplet)).to_string_lossy().to_string(),
            "--mpileup".to_string(),
            "1".to_string(),
            "--min-coverage".to_string(),
            cfg.somatic.min_coverage.to_string(),
            "--min-coverage-normal".to_string(),
            cfg.somatic.min_coverage_normal.to_string(),
            "--min-coverage-tumor".to_string(),
            cfg.somatic.min_coverage_tumor.to_string(),
            "--min-var-freq".to_string(),
            cfg.somatic.min_var_freq.to_string(),
            "--min-freq-for-hom".to_string(),
            cfg.somatic.min_freq_for_hom.to_string(),
            "--normal-purity".to_string(),
            cfg.somatic.normal_purity.to_string(),
            "--tumor-purity".to_string(),
            effective_purity.to_string(),
            "--p-value".to_string(),
            cfg.somatic.p_value.to_string(),
            "--somatic-p-value".to_string(),
            cfg.somatic.somatic_p_value.to_string(),
            "--strand-filter".to_string(),
            cfg.somatic.strand_filter.to_string(),
            "--min-base-qual".to_string(),
            cfg.somatic.min_base_qual.to_string(),
            "--output-vcf".to_string(),
            "1".to_string(),
        ];

        if !cfg.paths.vcf_sample_list.is_empty() && Path::new(&cfg.paths.vcf_sample_list).is_file() {
            args.push("--vcf-sample-list".to_string());
            args.push(cfg.paths.vcf_sample_list.clone());
        }

        let child = spawn_command("java", &args, None)?;
        push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
    }

    wait_all(&mut children)?;
    log_message("VarScan somatic calling complete");
    Ok(())
}

fn process_somatic_variants(paths: &Paths, cfg: &Config) -> AppResult<()> {
    log_message("=== STAGE 4: processSomatic — classifying into Somatic / Germline / LOH ===");
    let mut children = Vec::new();

    for snp_file in list_files_matching(&paths.somatic_dir, "*.snp.vcf")? {
        log_message(&format!("processSomatic SNP: {}", snp_file.file_name().and_then(OsStr::to_str).unwrap_or("")));
        let args = vec![
            "-jar".to_string(),
            format!("{}/VarScan.v2.3.9.jar", cfg.paths.software_dir),
            "processSomatic".to_string(),
            snp_file.to_string_lossy().to_string(),
            "--min-tumor-freq".to_string(),
            cfg.process_somatic.min_tumor_freq.to_string(),
            "--max-normal-freq".to_string(),
            cfg.process_somatic.max_normal_freq.to_string(),
            "--p-value".to_string(),
            cfg.process_somatic.p_value.to_string(),
        ];
        let child = spawn_command("java", &args, None)?;
        push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
    }

    for indel_file in list_files_matching(&paths.somatic_dir, "*.indel.vcf")? {
        log_message(&format!("processSomatic INDEL: {}", indel_file.file_name().and_then(OsStr::to_str).unwrap_or("")));
        let args = vec![
            "-jar".to_string(),
            format!("{}/VarScan.v2.3.9.jar", cfg.paths.software_dir),
            "processSomatic".to_string(),
            indel_file.to_string_lossy().to_string(),
            "--min-tumor-freq".to_string(),
            cfg.process_somatic.min_tumor_freq.to_string(),
            "--max-normal-freq".to_string(),
            cfg.process_somatic.max_normal_freq.to_string(),
            "--p-value".to_string(),
            cfg.process_somatic.p_value.to_string(),
        ];
        let child = spawn_command("java", &args, None)?;
        push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
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

fn run_varscan_copynumber(paths: &Paths, cfg: &Config) -> AppResult<()> {
    log_message("=== STAGE 5: VarScan copy number analysis ===");
    check_file_exists(Path::new(&cfg.paths.pairs_file))?;

    let pairs = read_pairs(cfg)?;
    let mut children = Vec::new();

    for (entry1, entry2, _) in pairs {
        let samplet = get_sample_name(&entry1, &cfg.paths.pairs_suffix);
        let samplen = get_sample_name(&entry2, &cfg.paths.pairs_suffix);

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
            format!("{}/VarScan.v2.3.9.jar", cfg.paths.software_dir),
            "copynumber".to_string(),
            mpileup_file.to_string_lossy().to_string(),
            paths.copy_number_dir.join(format!("{}_{}", samplen, samplet)).to_string_lossy().to_string(),
            "--mpileup".to_string(),
            "1".to_string(),
            "--min-coverage".to_string(),
            cfg.cnv.min_coverage.to_string(),
            "--min-base-qual".to_string(),
            cfg.somatic.min_base_qual.to_string(),
            "--min-segment-size".to_string(),
            cfg.cnv.min_segment_size.to_string(),
            "--max-segment-size".to_string(),
            cfg.cnv.max_segment_size.to_string(),
            "--p-value".to_string(),
            cfg.cnv.p_value.to_string(),
            "--data-ratio".to_string(),
            format!("{:.6}", dataratio),
        ];
        let child = spawn_command("java", &args, None)?;
        push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
    }

    wait_all(&mut children)?;
    log_message("VarScan copy number complete");
    Ok(())
}

fn run_copy_caller(paths: &Paths, cfg: &Config) -> AppResult<()> {
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
            format!("{}/VarScan.v2.3.9.jar", cfg.paths.software_dir),
            "copyCaller".to_string(),
            cnv_file.to_string_lossy().to_string(),
            "--output-file".to_string(),
            paths.copy_number_dir.join(format!("{}.copynumber.called", base)).to_string_lossy().to_string(),
            "--output-homdel-file".to_string(),
            paths.copy_number_dir.join(format!("{}.copynumber.homdel", base)).to_string_lossy().to_string(),
            "--amp-threshold".to_string(),
            cfg.cnv.amp_threshold.to_string(),
            "--del-threshold".to_string(),
            cfg.cnv.del_threshold.to_string(),
            "--recenter-up".to_string(),
            cfg.cnv.recenter_up.to_string(),
            "--recenter-down".to_string(),
            cfg.cnv.recenter_down.to_string(),
        ];

        let child = spawn_command("java", &args, None)?;
        push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
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

fn run_bam_readcount(paths: &Paths, cfg: &Config) -> AppResult<()> {
    log_message("=== STAGE 8: bam-readcount ===");
    check_file_exists(Path::new(&cfg.paths.pairs_file))?;

    let pairs = read_pairs(cfg)?;
    let mut children = Vec::new();

    for (entry1, entry2, _) in pairs {
        let samplet = get_sample_name(&entry1, &cfg.paths.pairs_suffix);
        let samplen = get_sample_name(&entry2, &cfg.paths.pairs_suffix);
        let tumor_bam = get_bam_path(paths, &entry1, cfg);
        let normal_bam = get_bam_path(paths, &entry2, cfg);
        let prefix = format!("{}_{}", samplen, samplet);

        check_bam_index(&tumor_bam)?;
        check_bam_index(&normal_bam)?;

        let var_file = paths.snp_var_dir.join(format!("{}.snp.Somatic.hc.var", prefix));
        if var_file.is_file() {
            log_message(&format!("bam-readcount: Somatic SNP — tumour {}", samplet));
            let args = vec![
                "-q".to_string(),
                cfg.readcount.map_qual.to_string(),
                "-b".to_string(),
                cfg.readcount.base_qual.to_string(),
                "-f".to_string(),
                cfg.paths.reference.clone(),
                "-l".to_string(),
                var_file.to_string_lossy().to_string(),
                tumor_bam.to_string_lossy().to_string(),
            ];
            let out = paths.readcount_dir.join(format!("{}.snp.Somatic.hc.readcount", prefix));
            let child = spawn_command("bam-readcount", &args, Some(&out))?;
            push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
        }

        let indel_var = paths.indel_var_dir.join(format!("{}.indel.Somatic.hc.var", prefix));
        if indel_var.is_file() {
            log_message(&format!("bam-readcount: Somatic INDEL — tumour {}", samplet));
            let args = vec![
                "-q".to_string(),
                cfg.readcount.map_qual.to_string(),
                "-b".to_string(),
                cfg.readcount.base_qual.to_string(),
                "-f".to_string(),
                cfg.paths.reference.clone(),
                "-l".to_string(),
                indel_var.to_string_lossy().to_string(),
                tumor_bam.to_string_lossy().to_string(),
            ];
            let out = paths.readcount_dir.join(format!("{}.indel.Somatic.hc.readcount", prefix));
            let child = spawn_command("bam-readcount", &args, Some(&out))?;
            push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
        }

        let germ_var = paths.snp_var_dir.join(format!("{}.snp.Germline.hc.var", prefix));
        if germ_var.is_file() {
            log_message(&format!("bam-readcount: Germline SNP — normal {}", samplen));
            let args = vec![
                "-q".to_string(),
                cfg.readcount.map_qual.to_string(),
                "-b".to_string(),
                cfg.readcount.base_qual.to_string(),
                "-f".to_string(),
                cfg.paths.reference.clone(),
                "-l".to_string(),
                germ_var.to_string_lossy().to_string(),
                normal_bam.to_string_lossy().to_string(),
            ];
            let out = paths.readcount_dir.join(format!("{}.snp.Germline.hc.readcount", prefix));
            let child = spawn_command("bam-readcount", &args, Some(&out))?;
            push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
        }

        let loh_var = paths.snp_var_dir.join(format!("{}.snp.LOH.hc.var", prefix));
        if loh_var.is_file() {
            log_message(&format!("bam-readcount: LOH SNP — tumour {}", samplet));
            let args = vec![
                "-q".to_string(),
                cfg.readcount.map_qual.to_string(),
                "-b".to_string(),
                cfg.readcount.base_qual.to_string(),
                "-f".to_string(),
                cfg.paths.reference.clone(),
                "-l".to_string(),
                loh_var.to_string_lossy().to_string(),
                tumor_bam.to_string_lossy().to_string(),
            ];
            let out = paths.readcount_dir.join(format!("{}.snp.LOH.hc.readcount", prefix));
            let child = spawn_command("bam-readcount", &args, Some(&out))?;
            push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
        }
    }

    wait_all(&mut children)?;
    log_message("bam-readcount complete");
    Ok(())
}

fn run_fpfilter(paths: &Paths, cfg: &Config) -> AppResult<()> {
    log_message("=== STAGE 9: False positive filtering (fpfilter.pl) ===");
    check_file_exists(Path::new(&format!("{}/fpfilter.pl", cfg.paths.scripts_dir)))?;

    let mut children = Vec::new();
    for rc_file in list_files_matching(&paths.readcount_dir, "*.readcount")? {
        let base = rc_file
            .file_name()
            .and_then(OsStr::to_str)
            .unwrap_or_default()
            .trim_end_matches(".readcount")
            .to_string();
        // readcount base is e.g. "samplen_samplet.snp.Somatic.hc"
        // matching VCF has .hc.vcf suffix: somatic/samplen_samplet.snp.Somatic.hc.vcf
        let vcf = paths.somatic_dir.join(format!("{}.vcf", base));
        if !vcf.is_file() {
            log_message(&format!("WARNING: No matching VCF for {} — skipping fpfilter", base));
            continue;
        }

        let filtered = paths.filtered_dir.join(format!("{}.fpfilter.vcf", base));
        log_message(&format!("fpfilter: {}", base));
        let args = vec![
            format!("{}/fpfilter.pl", cfg.paths.scripts_dir),
            "--vcf-var-file".to_string(),
            vcf.to_string_lossy().to_string(),
            "--readcount-file".to_string(),
            rc_file.to_string_lossy().to_string(),
            "--output-file".to_string(),
            filtered.to_string_lossy().to_string(),
            "--min-var-count".to_string(),
            "3".to_string(),
            "--min-var-freq".to_string(),
            cfg.somatic.min_var_freq.to_string(),
            "--min-ref-basequal".to_string(),
            cfg.readcount.base_qual.to_string(),
            "--min-var-basequal".to_string(),
            cfg.readcount.base_qual.to_string(),
        ];
        let child = spawn_command("perl", &args, None)?;
        push_child(&mut children, child, cfg.readcount.max_parallel_jobs)?;
    }

    wait_all(&mut children)?;
    log_message("False positive filtering complete");
    Ok(())
}

fn count_matches(dir: &Path, pattern: &str) -> AppResult<usize> {
    Ok(list_files_matching(dir, pattern)?.len())
}

fn generate_summary(paths: &Paths, cfg: &Config) -> AppResult<()> {
    log_message("=== STAGE 10: Generating summary report ===");

    let mut out = String::new();
    out.push_str("========================================================================\n");
    out.push_str("VarScan2 Pipeline Summary\n");
    out.push_str(&format!("Generated: {}\n", now_string()));
    out.push_str("========================================================================\n\n");
    out.push_str("CONFIGURATION\n");
    out.push_str(&format!("  Reference genome       : {}\n", cfg.paths.reference));
    out.push_str(&format!("  BAM directory          : {}\n", paths.bam_dir.display()));
    out.push_str(&format!("  BAM suffix             : {}\n", cfg.paths.bam_suffix));
    out.push_str(&format!("  Pairs file             : {}\n", cfg.paths.pairs_file));
    out.push_str(&format!("  Pairs suffix (stripped): {}\n", cfg.paths.pairs_suffix));
    out.push_str(&format!(
        "  Target regions (WES)   : {}\n",
        if cfg.paths.target_bed.is_empty() {
            "not set (full-genome pileup)".to_string()
        } else {
            cfg.paths.target_bed.clone()
        }
    ));
    out.push_str(&format!(
        "  VCF sample list        : {}\n\n",
        if cfg.paths.vcf_sample_list.is_empty() {
            "not set (generic NORMAL/TUMOR headers)".to_string()
        } else {
            cfg.paths.vcf_sample_list.clone()
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
    cfg: &Config,
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

    if args.resume && resume_stage_match(paths, stage, cfg)? {
        log_message(&format!("--- Stage {} skipped — SHA256 resume match: {}", stage, desc));
        return Ok(());
    }

    if args.dry_run {
        log_message(&format!("--- Stage {} [dry-run] would run: {}", stage, desc));
        return Ok(());
    }

    f()?;
    write_stage_marker(paths, stage, cfg)?;
    Ok(())
}

fn parse_args() -> AppResult<Args> {
    let mut from_stage:  u8 = 1;
    let mut to_stage:    u8 = 10;
    let mut resume       = false;
    let mut dry_run      = false;
    let mut config_path: Option<String> = None;
    let mut init_config  = false;
    let mut validate     = false;

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
       case01_final.bam,control01_final.bam
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
            "--config" => {
                let v = it.next().ok_or_else(|| "--config requires a path".to_string())?;
                config_path = Some(v);
            }
            "--init-config" => init_config = true,
            "--validate"    => validate = true,
            "--resume"  => resume = true,
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
        config_path,
        init_config,
        validate,
    })
}

fn build_paths(cfg: &Config) -> AppResult<Paths> {
    let cwd = env::current_dir().map_err(|e| format!("current_dir: {}", e))?;
    let bam_dir = if cfg.paths.bam_dir.is_empty() || cfg.paths.bam_dir == "." {
        cwd.clone()
    } else {
        PathBuf::from(&cfg.paths.bam_dir)
    };
    Ok(Paths {
        bam_dir,
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

    if args.init_config {
        print!("{}", config_example_toml());
        return Ok(());
    }

    let mut cfg = load_config(args.config_path.as_deref())?;
    apply_env_overrides(&mut cfg)?;

    if args.validate {
        let passed = validate_config(&cfg);
        std::process::exit(if passed { 0 } else { 1 });
    }

    let paths = build_paths(&cfg)?;

    log_message("Starting VarScan2 Pipeline (Rust)");
    log_message(&format!(
        "Stages: {}–{}  |  Resume: {}  |  Dry-run: {}",
        args.from_stage, args.to_stage, args.resume, args.dry_run
    ));

    if !args.dry_run {
        for cmd in ["samtools", "java", "bam-readcount", "perl"] {
            check_command_exists(cmd)?;
        }

        check_file_exists(Path::new(&cfg.paths.pairs_file))?;
        if cfg.paths.reference.is_empty() {
            return Err(
                "reference path is empty. Set it in config.toml ([paths] reference = \"...\") \
                 or via VARSCAN_REFERENCE env var."
                    .to_string(),
            );
        }
        check_file_exists(Path::new(&cfg.paths.reference))?;
        check_file_exists(Path::new(&format!("{}/VarScan.v2.3.9.jar", cfg.paths.software_dir)))?;
        check_file_exists(Path::new(&format!("{}/fpfilter.pl", cfg.paths.scripts_dir)))?;
    }

    setup_directories(&paths)?;

    run_stage(1, "samtools flagstat",    &args, &paths, &cfg, || generate_flagstats(&paths, &cfg))?;
    run_stage(2, "samtools mpileup",     &args, &paths, &cfg, || generate_mpileup(&paths, &cfg))?;
    run_stage(3, "VarScan somatic",      &args, &paths, &cfg, || run_varscan_somatic(&paths, &cfg))?;
    run_stage(4, "VarScan processSomatic", &args, &paths, &cfg, || process_somatic_variants(&paths, &cfg))?;

    if args.to_stage >= 7 && !args.dry_run {
        check_stage4_output(&paths)?;
    }

    run_stage(5, "VarScan copynumber",   &args, &paths, &cfg, || run_varscan_copynumber(&paths, &cfg))?;
    run_stage(6, "VarScan copyCaller",   &args, &paths, &cfg, || run_copy_caller(&paths, &cfg))?;
    run_stage(7, "prepare .var files",   &args, &paths, &cfg, || prepare_filter_input(&paths))?;
    run_stage(8, "bam-readcount",        &args, &paths, &cfg, || run_bam_readcount(&paths, &cfg))?;
    run_stage(9, "fpfilter.pl",          &args, &paths, &cfg, || run_fpfilter(&paths, &cfg))?;
    run_stage(10, "summary",             &args, &paths, &cfg, || generate_summary(&paths, &cfg))?;

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

    // ── config_example_toml ──────────────────────────────────────────────────

    #[test]
    fn init_config_example_produces_valid_toml() {
        let content = config_example_toml();
        let parsed: Config = toml::from_str(&content).unwrap();
        assert_eq!(parsed.somatic.min_coverage, 20);
        assert_eq!(parsed.readcount.max_parallel_jobs, 30);
        assert_eq!(parsed.paths.bam_suffix, "_final.bam");
    }

    // ── apply_env_overrides ──────────────────────────────────────────────────

    #[test]
    fn env_override_reference() {
        std::env::set_var("VARSCAN_REFERENCE", "/env/ref.fa");
        let mut cfg = Config::default();
        apply_env_overrides(&mut cfg).unwrap();
        assert_eq!(cfg.paths.reference, "/env/ref.fa");
        std::env::remove_var("VARSCAN_REFERENCE");
    }

    #[test]
    fn env_override_min_coverage_parses_int() {
        // uses a distinct var from env_override_invalid_int_errors to avoid parallel conflict
        std::env::set_var("VARSCAN_MIN_BASE_QUAL", "25");
        let mut cfg = Config::default();
        apply_env_overrides(&mut cfg).unwrap();
        assert_eq!(cfg.somatic.min_base_qual, 25);
        std::env::remove_var("VARSCAN_MIN_BASE_QUAL");
    }

    #[test]
    fn env_override_invalid_int_errors() {
        std::env::set_var("VARSCAN_MIN_SEGMENT_SIZE", "notanint");
        let mut cfg = Config::default();
        let result = apply_env_overrides(&mut cfg);
        std::env::remove_var("VARSCAN_MIN_SEGMENT_SIZE");
        assert!(result.is_err());
    }

    #[test]
    fn env_override_tumor_purity_parses_float() {
        std::env::set_var("VARSCAN_TUMOR_PURITY", "0.65");
        let mut cfg = Config::default();
        apply_env_overrides(&mut cfg).unwrap();
        assert!((cfg.somatic.tumor_purity - 0.65).abs() < 1e-9);
        std::env::remove_var("VARSCAN_TUMOR_PURITY");
    }

    // ── load_config ─────────────────────────────────────────────────────────

    #[test]
    fn load_config_from_toml_overrides_defaults() {
        let dir = std::env::temp_dir();
        let path = dir.join("varscan_test_config.toml");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "[somatic]\nmin_coverage = 30").unwrap();
        writeln!(f, "[paths]\nreference = \"/data/ref.fa\"").unwrap();
        let cfg = load_config(Some(path.to_str().unwrap())).unwrap();
        assert_eq!(cfg.somatic.min_coverage, 30);
        assert_eq!(cfg.paths.reference, "/data/ref.fa");
        assert_eq!(cfg.somatic.min_base_qual, 20); // non-overridden stays default
        std::fs::remove_file(path).ok();
    }

    #[test]
    fn load_config_no_file_returns_defaults() {
        let cfg = load_config(None).unwrap();
        assert_eq!(cfg.somatic.min_coverage, 20);
    }

    #[test]
    fn load_config_missing_explicit_file_errors() {
        let result = load_config(Some("/tmp/nonexistent_varscan_xyzzy.toml"));
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("not found"));
    }

    // ── Config defaults ─────────────────────────────────────────────────────

    #[test]
    fn config_defaults_match_original_consts() {
        let cfg = Config::default();
        assert_eq!(cfg.somatic.min_coverage, 20);
        assert_eq!(cfg.somatic.min_base_qual, 20);
        assert!((cfg.somatic.min_var_freq - 0.10).abs() < 1e-9);
        assert!((cfg.cnv.p_value - 0.01).abs() < 1e-9);
        assert_eq!(cfg.readcount.map_qual, 10);
        assert_eq!(cfg.readcount.base_qual, 15);
        assert_eq!(cfg.readcount.max_parallel_jobs, 30);
        assert_eq!(cfg.paths.bam_suffix, "_final.bam");
        assert_eq!(cfg.paths.pairs_file, "sample_pairs.csv");
    }

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
        assert_eq!(get_sample_name("case01_final.bam", "_final.bam"), "case01");
        assert_eq!(get_sample_name("sample_final.bam", "_final.bam"), "sample");
    }

    #[test]
    fn get_sample_name_no_suffix_unchanged() {
        assert_eq!(get_sample_name("sample.bam", "_final.bam"), "sample.bam");
        assert_eq!(get_sample_name("sample", "_final.bam"), "sample");
    }

    // ── read_pairs — per-sample purity ──────────────────────────────────────

    fn write_pairs_tempfile(content: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir();
        let path = dir.join(format!(
            "varscan_pairs_test_{}.csv",
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
    fn read_pairs_no_purity_column() {
        let p = write_pairs_tempfile("t1_final.bam,n1_final.bam\nt2_final.bam,n2_final.bam\n");
        let mut cfg = Config::default();
        cfg.paths.pairs_file = p.to_string_lossy().to_string();
        let pairs = read_pairs(&cfg).unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0], ("t1_final.bam".to_string(), "n1_final.bam".to_string(), None));
        assert_eq!(pairs[1], ("t2_final.bam".to_string(), "n2_final.bam".to_string(), None));
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn read_pairs_with_purity_column() {
        let p = write_pairs_tempfile("t1_final.bam,n1_final.bam,0.65\nt2_final.bam,n2_final.bam\n");
        let mut cfg = Config::default();
        cfg.paths.pairs_file = p.to_string_lossy().to_string();
        let pairs = read_pairs(&cfg).unwrap();
        assert_eq!(pairs.len(), 2);
        assert_eq!(pairs[0].2, Some(0.65));
        assert_eq!(pairs[1].2, None);
        std::fs::remove_file(p).ok();
    }

    #[test]
    fn read_pairs_invalid_purity_errors() {
        let p = write_pairs_tempfile("t1_final.bam,n1_final.bam,notanumber\n");
        let mut cfg = Config::default();
        cfg.paths.pairs_file = p.to_string_lossy().to_string();
        let result = read_pairs(&cfg);
        std::fs::remove_file(p).ok();
        assert!(result.is_err());
        assert!(result.unwrap_err().contains("invalid tumor_purity"));
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
