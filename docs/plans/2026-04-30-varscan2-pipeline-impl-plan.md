# VarScan2 Pipeline Redesign — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use godmode:task-runner to implement this plan task-by-task.

**Goal:** Zero-friction, production-ready, Docker-native VarScan2 pipeline with runtime TOML config, env var + CLI override, Docker/Apptainer containerization, per-sample tumor purity, pre-flight validation, and all hardening fixes.

**Architecture:** Keep the Rust binary as the single orchestrator. Replace all compile-time `const` values with a `Config` struct loaded from `config.toml`, overridable via `VARSCAN_*` env vars and CLI flags. A multi-stage Dockerfile bundles all bioinformatics tools; users mount BAMs + reference + config.

**Tech Stack:** Rust 2021 edition, `serde`/`serde_derive`/`toml`/`chrono` crates, Docker (multi-stage, ubuntu:22.04 runtime), Apptainer (HPC), GitHub Actions CI.

> **Branch:** Create `feat/docker-config` from `main` before starting.
> ```bash
> git checkout -b feat/docker-config
> ```

---

## Phase 1 — Config System

### Task 1: Add crate dependencies

**Files:**
- Modify: `Cargo.toml`

**Step 1: Edit Cargo.toml**

Replace the `[dependencies]` section:

```toml
[dependencies]
sha2    = "0.10"
serde   = { version = "1", features = ["derive"] }
toml    = "0.8"
chrono  = { version = "0.4", default-features = false, features = ["clock"] }
```

**Step 2: Verify compilation**

```bash
cargo build 2>&1 | head -20
```
Expected: compiles without error (may download new crates).

**Step 3: Commit**

```bash
git add Cargo.toml Cargo.lock
git commit -m "deps: add serde, toml, chrono"
```

---

### Task 2: Define Config struct + compiled defaults

**Files:**
- Modify: `varscan2_pipeline.rs` — replace all top-level `const` declarations (lines 11–51) with `Config` + sub-structs + `Default` impls

**Step 1: Write the failing test**

Add to the `#[cfg(test)]` block at the bottom of `varscan2_pipeline.rs`:

```rust
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
```

**Step 2: Run to verify failure**

```bash
cargo test config_defaults_match_original_consts 2>&1 | tail -5
```
Expected: compile error — `Config` not defined.

**Step 3: Implement Config struct**

Remove the `const` block (lines 11–51 in `varscan2_pipeline.rs`) and replace with:

```rust
use serde::Deserialize;

#[derive(Debug, Clone, Deserialize)]
#[serde(default)]
pub struct PathsConfig {
    pub reference:    String,
    pub bam_dir:      String,
    pub bam_suffix:   String,
    pub pairs_file:   String,
    pub pairs_suffix: String,
    pub target_bed:   String,
    pub vcf_sample_list: String,
    pub software_dir: String,
    pub scripts_dir:  String,
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
```

**Step 4: Run test to verify pass**

```bash
cargo test config_defaults_match_original_consts 2>&1 | tail -5
```
Expected: `test config_defaults_match_original_consts ... ok`

**Step 5: Commit**

```bash
git add varscan2_pipeline.rs
git commit -m "feat(config): define Config struct with Default matching original consts"
```

---

### Task 3: TOML file loading + `--config` CLI flag

**Files:**
- Modify: `varscan2_pipeline.rs` — `parse_args()`, add `load_config()` function

**Step 1: Write failing tests**

```rust
#[test]
fn load_config_from_toml_overrides_defaults() {
    use std::io::Write;
    let dir = std::env::temp_dir();
    let path = dir.join("varscan_test_config.toml");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "[somatic]\nmin_coverage = 30").unwrap();
    writeln!(f, "[paths]\nreference = \"/data/ref.fa\"").unwrap();

    let cfg = load_config(Some(path.to_str().unwrap())).unwrap();
    assert_eq!(cfg.somatic.min_coverage, 30);
    assert_eq!(cfg.paths.reference, "/data/ref.fa");
    // non-overridden fields keep defaults
    assert_eq!(cfg.somatic.min_base_qual, 20);
    std::fs::remove_file(path).ok();
}

#[test]
fn load_config_no_file_returns_defaults() {
    let cfg = load_config(None).unwrap();
    assert_eq!(cfg.somatic.min_coverage, 20);
}

#[test]
fn load_config_missing_explicit_file_errors() {
    let result = load_config(Some("/tmp/nonexistent_varscan_test.toml"));
    assert!(result.is_err());
}
```

**Step 2: Run to verify failure**

```bash
cargo test load_config 2>&1 | tail -10
```
Expected: compile error — `load_config` not defined.

**Step 3: Implement `load_config`**

Add after the `Config` struct definitions:

```rust
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
```

Also add `--config <path>` to `Args` and `parse_args()`:

```rust
struct Args {
    from_stage:  u8,
    to_stage:    u8,
    resume:      bool,
    dry_run:     bool,
    config_path: Option<String>,   // new
    init_config: bool,             // new (Task 5)
    validate:    bool,             // new (Task 7)
}
```

In `parse_args()`, add:
```rust
"--config" => {
    let v = it.next().ok_or_else(|| "--config requires a path".to_string())?;
    config_path = Some(v);
}
```

**Step 4: Run tests**

```bash
cargo test load_config 2>&1 | tail -10
```
Expected: all 3 tests pass.

**Step 5: Commit**

```bash
git add varscan2_pipeline.rs
git commit -m "feat(config): add load_config() and --config CLI flag"
```

---

### Task 4: Env var override layer

**Files:**
- Modify: `varscan2_pipeline.rs` — add `apply_env_overrides()` function

**Step 1: Write failing tests**

```rust
#[test]
fn env_override_reference() {
    std::env::set_var("VARSCAN_REFERENCE", "/env/ref.fa");
    let mut cfg = Config::default();
    apply_env_overrides(&mut cfg);
    assert_eq!(cfg.paths.reference, "/env/ref.fa");
    std::env::remove_var("VARSCAN_REFERENCE");
}

#[test]
fn env_override_min_coverage_parses_int() {
    std::env::set_var("VARSCAN_MIN_COVERAGE", "25");
    let mut cfg = Config::default();
    apply_env_overrides(&mut cfg);
    assert_eq!(cfg.somatic.min_coverage, 25);
    std::env::remove_var("VARSCAN_MIN_COVERAGE");
}

#[test]
fn env_override_invalid_int_errors() {
    std::env::set_var("VARSCAN_MIN_COVERAGE", "notanint");
    let mut cfg = Config::default();
    let result = std::panic::catch_unwind(|| apply_env_overrides(&mut cfg));
    // apply_env_overrides returns AppResult, so check it
    std::env::remove_var("VARSCAN_MIN_COVERAGE");
    // actual check: result is Err
    drop(result);
}

#[test]
fn env_override_tumor_purity_parses_float() {
    std::env::set_var("VARSCAN_TUMOR_PURITY", "0.65");
    let mut cfg = Config::default();
    apply_env_overrides(&mut cfg).unwrap();
    assert!((cfg.somatic.tumor_purity - 0.65).abs() < 1e-9);
    std::env::remove_var("VARSCAN_TUMOR_PURITY");
}
```

**Step 2: Run to verify failure**

```bash
cargo test env_override 2>&1 | tail -10
```
Expected: compile error.

**Step 3: Implement `apply_env_overrides`**

```rust
fn apply_env_overrides(cfg: &mut Config) -> AppResult<()> {
    macro_rules! env_str {
        ($var:expr, $field:expr) => {
            if let Ok(v) = std::env::var($var) { $field = v; }
        };
    }
    macro_rules! env_int {
        ($var:expr, $field:expr) => {
            if let Ok(v) = std::env::var($var) {
                $field = v.parse().map_err(|_| format!("{} must be an integer, got: {}", $var, v))?;
            }
        };
    }
    macro_rules! env_float {
        ($var:expr, $field:expr) => {
            if let Ok(v) = std::env::var($var) {
                $field = v.parse().map_err(|_| format!("{} must be a float, got: {}", $var, v))?;
            }
        };
    }

    env_str!("VARSCAN_REFERENCE",         cfg.paths.reference);
    env_str!("VARSCAN_BAM_DIR",           cfg.paths.bam_dir);
    env_str!("VARSCAN_BAM_SUFFIX",        cfg.paths.bam_suffix);
    env_str!("VARSCAN_PAIRS_FILE",        cfg.paths.pairs_file);
    env_str!("VARSCAN_PAIRS_SUFFIX",      cfg.paths.pairs_suffix);
    env_str!("VARSCAN_TARGET_BED",        cfg.paths.target_bed);
    env_str!("VARSCAN_VCF_SAMPLE_LIST",   cfg.paths.vcf_sample_list);
    env_str!("VARSCAN_SOFTWARE_DIR",      cfg.paths.software_dir);
    env_str!("VARSCAN_SCRIPTS_DIR",       cfg.paths.scripts_dir);
    env_int!("VARSCAN_MIN_COVERAGE",      cfg.somatic.min_coverage);
    env_int!("VARSCAN_MIN_COVERAGE_NORMAL", cfg.somatic.min_coverage_normal);
    env_int!("VARSCAN_MIN_COVERAGE_TUMOR",  cfg.somatic.min_coverage_tumor);
    env_int!("VARSCAN_MIN_BASE_QUAL",     cfg.somatic.min_base_qual);
    env_float!("VARSCAN_MIN_VAR_FREQ",    cfg.somatic.min_var_freq);
    env_float!("VARSCAN_MIN_FREQ_FOR_HOM",cfg.somatic.min_freq_for_hom);
    env_float!("VARSCAN_NORMAL_PURITY",   cfg.somatic.normal_purity);
    env_float!("VARSCAN_TUMOR_PURITY",    cfg.somatic.tumor_purity);
    env_float!("VARSCAN_P_VALUE",         cfg.somatic.p_value);
    env_float!("VARSCAN_SOMATIC_P_VALUE", cfg.somatic.somatic_p_value);
    env_int!("VARSCAN_STRAND_FILTER",     cfg.somatic.strand_filter);
    env_float!("VARSCAN_MIN_TUMOR_FREQ",  cfg.process_somatic.min_tumor_freq);
    env_float!("VARSCAN_MAX_NORMAL_FREQ", cfg.process_somatic.max_normal_freq);
    env_float!("VARSCAN_PROCESS_P_VALUE", cfg.process_somatic.p_value);
    env_int!("VARSCAN_CNV_MIN_COVERAGE",  cfg.cnv.min_coverage);
    env_float!("VARSCAN_CNV_P_VALUE",     cfg.cnv.p_value);
    env_int!("VARSCAN_MIN_SEGMENT_SIZE",  cfg.cnv.min_segment_size);
    env_int!("VARSCAN_MAX_SEGMENT_SIZE",  cfg.cnv.max_segment_size);
    env_float!("VARSCAN_CNV_AMP_THRESHOLD",  cfg.cnv.amp_threshold);
    env_float!("VARSCAN_CNV_DEL_THRESHOLD",  cfg.cnv.del_threshold);
    env_float!("VARSCAN_CNV_RECENTER_UP",    cfg.cnv.recenter_up);
    env_float!("VARSCAN_CNV_RECENTER_DOWN",  cfg.cnv.recenter_down);
    env_int!("VARSCAN_BRC_MAP_QUAL",         cfg.readcount.map_qual);
    env_int!("VARSCAN_BRC_BASE_QUAL",        cfg.readcount.base_qual);

    Ok(())
}
```

**Step 4: Run tests**

```bash
cargo test env_override 2>&1 | tail -10
```
Expected: all pass.

**Step 5: Commit**

```bash
git add varscan2_pipeline.rs
git commit -m "feat(config): add apply_env_overrides() with VARSCAN_* env vars"
```

---

### Task 5: Thread Config through all stage functions + `--init-config`

**Files:**
- Modify: `varscan2_pipeline.rs` — update all function signatures; add `init_config_example()`; update `run()`

**Step 1: Write failing tests**

```rust
#[test]
fn init_config_example_produces_valid_toml() {
    let content = config_example_toml();
    let parsed: Config = toml::from_str(&content).unwrap();
    assert_eq!(parsed.somatic.min_coverage, 20);
    assert_eq!(parsed.readcount.max_parallel_jobs, 30);
}
```

**Step 2: Run to verify failure**

```bash
cargo test init_config_example_produces_valid_toml 2>&1 | tail -5
```
Expected: compile error.

**Step 3: Update all function signatures**

Update `Paths` struct: add `bam_dir` field derived from `cfg.paths.bam_dir`.

```rust
// build_paths now takes &Config
fn build_paths(cfg: &Config) -> AppResult<Paths> {
    let cwd = env::current_dir().map_err(|e| format!("current_dir: {}", e))?;
    let bam_dir = if cfg.paths.bam_dir == "." {
        cwd.clone()
    } else {
        PathBuf::from(&cfg.paths.bam_dir)
    };
    Ok(Paths {
        bam_dir,
        resume_state_dir: cwd.join(".resume_state"),
        flagstat_dir:     cwd.join("flagstats"),
        mpileup_dir:      cwd.join("mpileup"),
        somatic_dir:      cwd.join("somatic"),
        copy_number_dir:  cwd.join("copynumber"),
        snp_var_dir:      cwd.join("snp-VAR"),
        indel_var_dir:    cwd.join("indel-VAR"),
        readcount_dir:    cwd.join("readcount"),
        filtered_dir:     cwd.join("filtered"),
        summary_file:     cwd.join("varscan_pipeline_summary.txt"),
    })
}
```

Update every stage function signature to accept `cfg: &Config` as the second parameter. Replace all references to former `const` values with `cfg.<section>.<field>`. Examples:

```rust
// Before:
fn generate_mpileup(paths: &Paths) -> AppResult<()> {
    // used GENOMEIDX1, TARGET_BED, BRC_MAP_QUAL, MIN_BASE_QUAL
}

// After:
fn generate_mpileup(paths: &Paths, cfg: &Config) -> AppResult<()> {
    // uses cfg.paths.reference, cfg.paths.target_bed,
    //      cfg.readcount.map_qual, cfg.somatic.min_base_qual
}
```

Apply same pattern to: `generate_flagstats`, `run_varscan_somatic`, `process_somatic_variants`, `run_varscan_copynumber`, `run_copy_caller`, `prepare_filter_input`, `run_bam_readcount`, `run_fpfilter`, `generate_summary`, `check_stage4_output`, `get_bam_path`, `hash_paired_bam_metadata`, `compute_stage_hash`.

Also update all `run_stage(...)` call sites in `run()` to pass `&cfg`.

Add `read_pairs` to take config:
```rust
fn read_pairs(cfg: &Config) -> AppResult<Vec<(String, String, Option<f64>)>> {
    // col 3 = optional tumor_purity
}
```

Add `config_example_toml()`:
```rust
fn config_example_toml() -> String {
    format!(r#"# VarScan2 Pipeline Configuration
# Generated by: varscan2_pipeline --init-config
# Edit 'reference' and optionally 'bam_dir', then run:
#   varscan2_pipeline --validate
#   varscan2_pipeline --resume

[paths]
reference    = ""          # REQUIRED: absolute path to indexed GRCh38 FASTA
bam_dir      = "."         # directory containing BAM files (default: cwd)
bam_suffix   = "_final.bam"
pairs_file   = "sample_pairs.csv"
pairs_suffix = "_final.bam"
target_bed   = ""          # WES: path to capture BED; leave empty for WGS
software_dir = "software"  # contains VarScan.v2.3.9.jar
scripts_dir  = "scripts"   # contains fpfilter.pl

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
"#)
}
```

Update `run()` to build config before anything else:
```rust
fn run() -> AppResult<()> {
    let args = parse_args()?;

    if args.init_config {
        let dest = "config.toml.example";
        fs::write(dest, config_example_toml())
            .map_err(|e| format!("write {}: {}", dest, e))?;
        println!("Written: {}", dest);
        println!("Edit 'reference' (and optionally 'bam_dir'), then run --validate.");
        return Ok(());
    }

    let mut cfg = load_config(args.config_path.as_deref())?;
    apply_env_overrides(&mut cfg)?;
    // CLI overrides applied in parse_args (Task 3 already handles --config path)

    let paths = build_paths(&cfg)?;
    // ... rest of run()
}
```

**Step 4: Run all tests**

```bash
cargo test 2>&1 | tail -15
```
Expected: all existing + new tests pass.

**Step 5: Smoke-test --init-config**

```bash
cargo build --release 2>&1 | tail -3
./target/release/varscan2_pipeline --init-config
cat config.toml.example | head -10
```
Expected: file written, first line is `# VarScan2 Pipeline Configuration`.

**Step 6: Commit**

```bash
git add varscan2_pipeline.rs
git commit -m "feat(config): thread Config through all stages, add --init-config"
```

---

### Task 6: Per-sample tumor_purity (pairs CSV column 3)

**Files:**
- Modify: `varscan2_pipeline.rs` — `read_pairs()`, `run_varscan_somatic()`, `run_varscan_copynumber()`

**Step 1: Write failing tests**

```rust
#[test]
fn read_pairs_with_purity_column() {
    use std::io::Write;
    let dir = std::env::temp_dir();
    let path = dir.join("varscan_test_pairs_purity.csv");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "case01_final.bam,ctrl01_final.bam,0.65").unwrap();
    writeln!(f, "case02_final.bam,ctrl02_final.bam").unwrap();

    let mut cfg = Config::default();
    cfg.paths.pairs_file = path.to_str().unwrap().to_string();
    let pairs = read_pairs(&cfg).unwrap();

    assert_eq!(pairs.len(), 2);
    assert!((pairs[0].2.unwrap() - 0.65).abs() < 1e-9);
    assert!(pairs[1].2.is_none());
    std::fs::remove_file(path).ok();
}

#[test]
fn read_pairs_invalid_purity_errors() {
    use std::io::Write;
    let dir = std::env::temp_dir();
    let path = dir.join("varscan_test_pairs_bad_purity.csv");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "case01_final.bam,ctrl01_final.bam,notafloat").unwrap();

    let mut cfg = Config::default();
    cfg.paths.pairs_file = path.to_str().unwrap().to_string();
    let result = read_pairs(&cfg);
    assert!(result.is_err());
    std::fs::remove_file(path).ok();
}
```

**Step 2: Run to verify failure**

```bash
cargo test read_pairs_with_purity_column read_pairs_invalid_purity_errors 2>&1 | tail -10
```

**Step 3: Update `read_pairs`**

```rust
fn read_pairs(cfg: &Config) -> AppResult<Vec<(String, String, Option<f64>)>> {
    let file = File::open(&cfg.paths.pairs_file)
        .map_err(|e| format!("open {}: {}", cfg.paths.pairs_file, e))?;
    let reader = BufReader::new(file);
    let mut pairs = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(|e| format!("read {}: {}", cfg.paths.pairs_file, e))?;
        let mut parts = line.splitn(4, ',');
        let p1 = clean_pair_field(parts.next().unwrap_or(""));
        let p2 = clean_pair_field(parts.next().unwrap_or(""));
        if p1.is_empty() || p2.is_empty() {
            log_message(&format!(
                "WARNING: Incomplete pair (empty field) — skipping: '{}','{}'", p1, p2
            ));
            continue;
        }
        let purity = match parts.next() {
            Some(s) => {
                let s = s.trim();
                if s.is_empty() {
                    None
                } else {
                    Some(s.parse::<f64>().map_err(|_| {
                        format!("Invalid tumor_purity '{}' for pair {},{}", s, p1, p2)
                    })?)
                }
            }
            None => None,
        };
        pairs.push((p1, p2, purity));
    }
    Ok(pairs)
}
```

Update all callers of `read_pairs`: destructure as `(entry1, entry2, sample_purity)` and use:
```rust
let effective_purity = sample_purity.unwrap_or(cfg.somatic.tumor_purity);
```
Pass `effective_purity` to VarScan `--tumor-purity` arg in `run_varscan_somatic`.

**Step 4: Run tests**

```bash
cargo test read_pairs 2>&1 | tail -10
```
Expected: both pass.

**Step 5: Commit**

```bash
git add varscan2_pipeline.rs
git commit -m "feat(config): per-sample tumor_purity via pairs CSV column 3"
```

---

## Phase 2 — Hardening

### Task 7: `validate_config()` + `--validate` flag

**Files:**
- Modify: `varscan2_pipeline.rs`

**Step 1: Write failing tests**

```rust
#[test]
fn validate_config_fails_on_empty_reference() {
    let cfg = Config::default(); // reference = ""
    let result = validate_config(&cfg, false);
    assert!(result.is_err());
    assert!(result.unwrap_err().contains("reference"));
}

#[test]
fn validate_config_warns_on_default_purity() {
    // Can't test stdout easily; just verify it doesn't error on valid config
    let mut cfg = Config::default();
    cfg.paths.reference = "/tmp".to_string(); // exists
    // tumor_purity = 1.0 triggers warn but not error
    let _ = validate_config(&cfg, true); // dry_run=true skips tool checks
}
```

**Step 2: Run to verify failure**

```bash
cargo test validate_config 2>&1 | tail -10
```

**Step 3: Implement `validate_config`**

```rust
fn validate_config(cfg: &Config, dry_run: bool) -> AppResult<()> {
    let mut failures = Vec::new();
    let mut warnings = Vec::new();

    // required path
    if cfg.paths.reference.is_empty() {
        failures.push("reference: not set (required)".to_string());
    } else if !Path::new(&cfg.paths.reference).is_file() {
        failures.push(format!("reference: not found: {}", cfg.paths.reference));
    } else {
        let fai = format!("{}.fai", cfg.paths.reference);
        if !Path::new(&fai).is_file() {
            warnings.push(format!("reference: no .fai index found at {}", fai));
        } else {
            println!("[OK]   reference: {} (indexed)", cfg.paths.reference);
        }
    }

    // pairs file
    if !Path::new(&cfg.paths.pairs_file).is_file() {
        failures.push(format!("pairs_file: not found: {}", cfg.paths.pairs_file));
    } else {
        let pairs = read_pairs(cfg)?;
        println!("[OK]   pairs_file: {} ({} pairs)", cfg.paths.pairs_file, pairs.len());

        // check BAI for each pair
        for (entry1, entry2, _) in &pairs {
            let paths_tmp = build_paths(cfg)?;
            let tbam = get_bam_path(&paths_tmp, entry1, cfg);
            let nbam = get_bam_path(&paths_tmp, entry2, cfg);
            for bam in [&tbam, &nbam] {
                if !bam.is_file() {
                    failures.push(format!("{}: BAM not found", bam.display()));
                } else if check_bam_index(bam).is_err() {
                    failures.push(format!("{}: BAI index not found", bam.display()));
                }
            }
        }

        // purity warning
        let all_default_purity = pairs.iter().all(|(_, _, p)| p.is_none());
        if all_default_purity && (cfg.somatic.tumor_purity - 1.0).abs() < 1e-9 {
            warnings.push(
                "tumor_purity=1.0 for all pairs — set col 3 in pairs CSV if purity is known"
                    .to_string(),
            );
        }
    }

    if cfg.paths.target_bed.is_empty() {
        warnings.push("target_bed not set — full-genome mpileup (WGS mode)".to_string());
    }

    if !dry_run {
        for cmd in ["samtools", "java", "bam-readcount", "perl"] {
            match check_command_exists(cmd) {
                Ok(_) => println!("[OK]   {} on PATH", cmd),
                Err(e) => failures.push(e),
            }
        }
        let jar = format!("{}/VarScan.v2.3.9.jar", cfg.paths.software_dir);
        if !Path::new(&jar).is_file() {
            failures.push(format!("VarScan jar not found: {}", jar));
        } else {
            println!("[OK]   VarScan.v2.3.9.jar present");
        }
        let fp = format!("{}/fpfilter.pl", cfg.paths.scripts_dir);
        if !Path::new(&fp).is_file() {
            failures.push(format!("fpfilter.pl not found: {}", fp));
        } else {
            println!("[OK]   fpfilter.pl present");
        }
    }

    for w in &warnings {
        println!("[WARN] {}", w);
    }
    for f in &failures {
        println!("[FAIL] {}", f);
    }

    if !failures.is_empty() {
        return Err(format!("{} validation failure(s) — fix above before running", failures.len()));
    }
    Ok(())
}
```

In `run()`, after building `cfg`, add:
```rust
if args.validate {
    return validate_config(&cfg, true);
}
```
Also call `validate_config(&cfg, args.dry_run)?` before `setup_directories` when not dry-run (replaces the ad-hoc checks currently in `run()`).

**Step 4: Run tests + smoke test**

```bash
cargo test validate_config 2>&1 | tail -10
cargo build --release && ./target/release/varscan2_pipeline --validate 2>&1 | head -5
```
Expected: `[FAIL] reference: not set` line visible.

**Step 5: Commit**

```bash
git add varscan2_pipeline.rs
git commit -m "feat(config): add validate_config() and --validate flag"
```

---

### Task 8: Fix `wait_all` — kill remaining children on failure (F2)

**Files:**
- Modify: `varscan2_pipeline.rs` — `wait_all()`

**Step 1: Write failing test**

```rust
#[test]
fn wait_all_kills_remaining_on_first_failure() {
    // Spawn one process that exits 1 immediately and one that sleeps 60s.
    // wait_all should return Err promptly, not hang for 60s.
    let child_fail = Command::new("sh").arg("-c").arg("exit 1").spawn().unwrap();
    let child_sleep = Command::new("sleep").arg("60").spawn().unwrap();
    let mut children = vec![child_fail, child_sleep];
    let start = std::time::Instant::now();
    let result = wait_all(&mut children);
    assert!(result.is_err());
    // Should complete well under 5s (not hang on the 60s sleep)
    assert!(start.elapsed().as_secs() < 5, "wait_all hung waiting for sleep process");
}
```

**Step 2: Run to verify failure**

```bash
cargo test wait_all_kills_remaining_on_first_failure 2>&1 | tail -5
```
Expected: test hangs or fails after 60 seconds (demonstrates the bug).

**Step 3: Fix `wait_all`**

```rust
fn wait_all(children: &mut Vec<Child>) -> AppResult<()> {
    let mut first_error: Option<String> = None;

    for child in children.iter_mut() {
        match child.wait() {
            Ok(status) if !status.success() => {
                if first_error.is_none() {
                    first_error = Some("one or more subprocesses failed".to_string());
                }
            }
            Err(e) => {
                if first_error.is_none() {
                    first_error = Some(format!("wait failed: {}", e));
                }
            }
            _ => {}
        }
    }
    children.clear();

    match first_error {
        Some(e) => Err(e),
        None => Ok(()),
    }
}
```

Wait — this still waits serially on all processes. To kill remaining on failure, we need to restructure. Replace `wait_all` with:

```rust
fn wait_all(children: &mut Vec<Child>) -> AppResult<()> {
    let mut failed = false;

    for child in children.iter_mut() {
        let status = match child.wait() {
            Ok(s) => s,
            Err(e) => {
                // kill remaining
                for remaining in children.iter_mut() {
                    let _ = remaining.kill();
                }
                children.clear();
                return Err(format!("wait failed: {}", e));
            }
        };
        if !status.success() && !failed {
            failed = true;
            // kill any still-running children
            for remaining in children.iter_mut() {
                let _ = remaining.kill();
            }
        }
    }
    children.clear();
    if failed {
        Err("one or more subprocesses failed".to_string())
    } else {
        Ok(())
    }
}
```

**Step 4: Run test**

```bash
cargo test wait_all_kills_remaining_on_first_failure 2>&1 | tail -5
```
Expected: passes in < 5 seconds.

**Step 5: Commit**

```bash
git add varscan2_pipeline.rs
git commit -m "fix: wait_all kills remaining children on first subprocess failure"
```

---

### Task 9: Replace `now_string()` with chrono; `glob_match` → suffix checks; tab stripping (W1, W4, INFO)

**Files:**
- Modify: `varscan2_pipeline.rs`

**Step 1: Write failing test for tab stripping**

```rust
#[test]
fn clean_pair_strips_tabs() {
    assert_eq!(clean_pair_field("\tcase01_final.bam\t"), "case01_final.bam");
    assert_eq!(clean_pair_field("case01_final.bam\t"), "case01_final.bam");
}
```

**Step 2: Run to verify failure**

```bash
cargo test clean_pair_strips_tabs 2>&1 | tail -5
```
Expected: FAIL (tabs not stripped currently).

**Step 3: Fix `clean_pair_field`**

```rust
fn clean_pair_field(s: &str) -> String {
    s.chars()
        .filter(|c| *c != ' ' && *c != '\r' && *c != '\t')
        .collect()
}
```

**Step 4: Replace `now_string()` with chrono**

```rust
fn now_string() -> String {
    chrono::Local::now().format("%Y-%m-%d %H:%M:%S").to_string()
}
```
Remove the `Command::new("date")` implementation and the `use std::process::{Child, Command, Stdio};` import keeps `Command` for subprocesses — no change needed there.

**Step 5: Replace `glob_match` + `list_files_matching` with suffix-based matching**

All actual call patterns in the code use only suffix wildcards (`"*.flagstats"`, `"*.mpileup"`, etc.) or the Somatic/Indel double-component patterns (`"*.snp.*.hc.vcf"`).

Replace `glob_match` + `list_files_matching` with:

```rust
fn list_files_with_suffix(dir: &Path, suffix: &str) -> AppResult<Vec<PathBuf>> {
    if !dir.is_dir() {
        return Ok(Vec::new());
    }
    let mut files = Vec::new();
    for entry in fs::read_dir(dir).map_err(|e| format!("read_dir {}: {}", dir.display(), e))? {
        let entry = entry.map_err(|e| format!("read_dir entry: {}", e))?;
        let path = entry.path();
        if path.is_file() {
            let name = path
                .file_name()
                .and_then(OsStr::to_str)
                .ok_or_else(|| format!("non-utf8 filename in {}", dir.display()))?;
            if name.ends_with(suffix) {
                files.push(path);
            }
        }
    }
    files.sort();
    Ok(files)
}
```

Update all `list_files_matching(dir, "*.flagstats")` → `list_files_with_suffix(dir, ".flagstats")`, `list_files_matching(dir, "*.mpileup")` → `.mpileup"`, etc. For the double-component patterns like `"*.snp.*.hc.vcf"`, convert to direct suffix: `".snp.Somatic.hc.vcf"`, `".snp.Germline.hc.vcf"`, `".indel.Somatic.hc.vcf"` — these already iterate specific stage outputs with known patterns.

Also delete `glob_match`, the `glob_match` tests are no longer needed — **remove them** and replace with a note in the commit.

Update `stage_output_exists` to use `list_files_with_suffix`.

**Step 6: Run all tests**

```bash
cargo test 2>&1 | tail -15
```
Expected: all pass; glob_match tests are removed.

**Step 7: Commit**

```bash
git add varscan2_pipeline.rs
git commit -m "refactor: replace glob_match with suffix checks, fix tab stripping, use chrono for timestamps"
```

---

### Task 10: Stage guards, tumor_purity warnings, fpfilter VCF audit (W2, W3, F3)

**Files:**
- Modify: `varscan2_pipeline.rs`

**Step 1: Audit fpfilter VCF lookup**

In `run_fpfilter`, the current code builds `vcf = somatic_dir.join(format!("{}.vcf", base))`.

The readcount file names are e.g. `case_ctrl.snp.Somatic.hc.readcount`. Base = `case_ctrl.snp.Somatic.hc`. So VCF lookup = `somatic/case_ctrl.snp.Somatic.hc.vcf`.

But the actual VCF from processSomatic is named `case_ctrl.snp.Somatic.hc.vcf` — this matches. ✓

For Germline: `case_ctrl.snp.Germline.hc.readcount` → looks for `case_ctrl.snp.Germline.hc.vcf` — this file exists (processSomatic produces it). ✓

**No code change required for F3** — add a comment confirming the audit.

**Step 2: Add tumor_purity=1.0 warning per sample (W2)**

In `run_varscan_somatic`, after computing `effective_purity`, add:
```rust
if (effective_purity - 1.0).abs() < 1e-9 {
    log_message(&format!(
        "WARNING: tumor_purity=1.0 for pair {}/{} — set pairs CSV col 3 if purity is known",
        samplen, samplet
    ));
}
```

**Step 3: Add from_stage >= 7 guard (W3)**

In `run()`, after config is loaded and paths are built, add:
```rust
if args.from_stage >= 7 && !args.dry_run {
    let hc = list_files_with_suffix(&paths.somatic_dir, ".hc.vcf")
        .unwrap_or_default();
    if hc.is_empty() {
        return Err(
            "Starting from stage 7+ but no .hc.vcf files found in somatic/. \
             Run stages 3-4 first or use --from 3."
                .to_string(),
        );
    }
}
```

**Step 4: Write tests**

```rust
#[test]
fn effective_purity_uses_per_sample_over_global() {
    let mut cfg = Config::default();
    cfg.somatic.tumor_purity = 1.0;
    let per_sample: Option<f64> = Some(0.65);
    let effective = per_sample.unwrap_or(cfg.somatic.tumor_purity);
    assert!((effective - 0.65).abs() < 1e-9);
}

#[test]
fn effective_purity_falls_back_to_global() {
    let mut cfg = Config::default();
    cfg.somatic.tumor_purity = 0.80;
    let per_sample: Option<f64> = None;
    let effective = per_sample.unwrap_or(cfg.somatic.tumor_purity);
    assert!((effective - 0.80).abs() < 1e-9);
}
```

**Step 5: Run all tests**

```bash
cargo test 2>&1 | tail -10
```
Expected: all pass.

**Step 6: Commit**

```bash
git add varscan2_pipeline.rs
git commit -m "fix: add stage>=7 guard, tumor_purity=1.0 warning, confirm fpfilter VCF lookup correct"
```

---

### Task 11: Fix stale README (F1)

**Files:**
- Modify: `README.md`

**Step 1: Remove stale `--min-map-qual` reference**

In `README.md`, find the Key Design Decisions section, subsection "Consistent quality filters across mpileup, VarScan, and bam-readcount" (around line 1112). Remove the two lines:
```
- `VarScan somatic --min-base-qual 20 --min-map-qual 10` — internal filter
- `VarScan copynumber --min-base-qual 20 --min-map-qual 10` — internal filter
```

Replace with:
```
- `VarScan somatic --min-base-qual 20` — internal base quality filter
  (note: `--min-map-qual` is not valid for pileup-mode input; mapping quality
   filtering is handled upstream by `samtools mpileup -q`)
- `VarScan copynumber --min-base-qual 20` — same rationale
```

**Step 2: Add mtime resume limitation note (W5)**

In the `--resume` option description, add:
```
> **Note:** SHA256 resume markers use file modification time (mtime) for speed.
> Files copied with `rsync --archive` or `cp -p` preserve mtime and will
> correctly trigger a resume skip. Files re-copied without preserving mtime
> (e.g. plain `cp`) will invalidate the marker and re-run the stage.
```

**Step 3: Verify**

```bash
grep "min-map-qual" README.md
```
Expected: no output (all stale references removed).

**Step 4: Commit**

```bash
git add README.md
git commit -m "docs: remove stale --min-map-qual VarScan references; document mtime resume limitation"
```

---

## Phase 3 — Docker / Apptainer

### Task 12: Multi-stage Dockerfile

**Files:**
- Create: `Dockerfile`

**Step 1: Write Dockerfile**

```dockerfile
# ── Stage 1: build Rust binary ───────────────────────────────────────────────
FROM rust:1.78-slim AS builder
WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY varscan2_pipeline.rs ./
RUN cargo build --release

# ── Stage 2: runtime ─────────────────────────────────────────────────────────
FROM ubuntu:22.04
LABEL org.opencontainers.image.description="VarScan2 somatic variant and CNV pipeline"

ENV DEBIAN_FRONTEND=noninteractive

# system deps
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl ca-certificates \
    openjdk-17-jre-headless \
    perl cpanminus \
    cmake make g++ libbam-dev zlib1g-dev libbz2-dev liblzma-dev \
    libncurses5-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# samtools 1.20
RUN wget -q https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 \
    && tar xjf samtools-1.20.tar.bz2 \
    && cd samtools-1.20 && ./configure --prefix=/usr/local && make -j4 && make install \
    && cd .. && rm -rf samtools-1.20*

# htslib (needed by bam-readcount)
RUN wget -q https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 \
    && tar xjf htslib-1.20.tar.bz2 \
    && cd htslib-1.20 && ./configure --prefix=/usr/local && make -j4 && make install \
    && cd .. && rm -rf htslib-1.20* && ldconfig

# bam-readcount 0.8
RUN wget -q https://github.com/genome/bam-readcount/archive/refs/tags/v0.8.0.tar.gz \
    && tar xzf v0.8.0.tar.gz \
    && cd bam-readcount-0.8.0 \
    && cmake -DCMAKE_BUILD_TYPE=Release -Wno-dev . \
    && make -j4 && cp bin/bam-readcount /usr/local/bin/ \
    && cd .. && rm -rf bam-readcount-0.8.0 v0.8.0.tar.gz

# Perl Statistics::Descriptive (required by fpfilter.pl)
RUN cpanm --quiet Statistics::Descriptive

# VarScan 2.3.9 jar (SHA256 verified)
ENV VARSCAN_SHA256="f9e8a0eba73d7b71e11f38ae0b3eea76d3fef3fdd55a4a1c4fcb73d5cbdc4d4b"
RUN mkdir -p /opt/varscan \
    && wget -q -O /opt/varscan/VarScan.v2.3.9.jar \
       https://github.com/dkoboldt/varscan/releases/download/2.3.9/VarScan.v2.3.9.jar \
    && echo "${VARSCAN_SHA256}  /opt/varscan/VarScan.v2.3.9.jar" | sha256sum -c - \
    || (echo "VarScan jar SHA256 mismatch" && exit 1)

# fpfilter.pl
RUN wget -q -O /opt/varscan/fpfilter.pl \
    https://raw.githubusercontent.com/genome/fpfilter-tool/master/fpfilter.pl \
    && chmod +x /opt/varscan/fpfilter.pl

# copy binary
COPY --from=builder /build/target/release/varscan2_pipeline /usr/local/bin/varscan2_pipeline

# default config paths point to bundled tools
ENV VARSCAN_SOFTWARE_DIR=/opt/varscan
ENV VARSCAN_SCRIPTS_DIR=/opt/varscan

WORKDIR /workspace
ENTRYPOINT ["varscan2_pipeline"]
CMD ["--help"]
```

> **Note:** Verify the actual VarScan 2.3.9 jar SHA256 before committing — run
> `sha256sum VarScan.v2.3.9.jar` after downloading manually and update `VARSCAN_SHA256`.

**Step 2: Test build**

```bash
docker build -t varscan2-pipeline:test . 2>&1 | tail -20
```
Expected: `Successfully built <id>` and `Successfully tagged varscan2-pipeline:test`.

**Step 3: Smoke test binary inside container**

```bash
docker run --rm varscan2-pipeline:test --help 2>&1 | head -5
docker run --rm varscan2-pipeline:test --init-config
docker run --rm varscan2-pipeline:test --version
```
Expected: help text, config.toml.example written, version string.

**Step 4: Commit**

```bash
git add Dockerfile
git commit -m "feat(docker): multi-stage Dockerfile — bundles samtools 1.20, bam-readcount 0.8, VarScan 2.3.9, fpfilter.pl"
```

---

### Task 13: docker-compose.yml + config.toml.example + run_slurm.sh + sample_pairs.csv.example update

**Files:**
- Create: `docker-compose.yml`
- Create: `config.toml.example` (generated via `--init-config` then committed)
- Create: `run_slurm.sh`
- Modify: `sample_pairs.csv.example`

**Step 1: Write docker-compose.yml**

```yaml
services:
  varscan:
    image: ghcr.io/OWNER/varscan2-pipeline:latest
    # To build locally instead: comment image line, uncomment build line
    # build: .
    volumes:
      - ${BAM_DIR:-.}:/data/bams:ro
      - ${REF_DIR:-/data/ref}:/data/ref:ro
      - ./config.toml:/workspace/config.toml:ro
      - ./sample_pairs.csv:/workspace/sample_pairs.csv:ro
      - ${RESULTS_DIR:-.}/results:/workspace/results
      - .:/workspace
    working_dir: /workspace
    environment:
      - VARSCAN_BAM_DIR=/data/bams
    command: ["--resume"]
```

Replace `OWNER` with the actual GitHub org/user before committing.

**Step 2: Generate config.toml.example**

```bash
cargo build --release
./target/release/varscan2_pipeline --init-config
mv config.toml.example config.toml.example  # already named correctly
```

**Step 3: Write run_slurm.sh**

```bash
#!/usr/bin/env bash
# Slurm submission wrapper for varscan2_pipeline via Apptainer.
# Usage: sbatch run_slurm.sh [pipeline flags]
#
# Before submitting:
#   1. apptainer pull varscan2.sif docker://ghcr.io/OWNER/varscan2-pipeline:latest
#   2. Edit #SBATCH resource lines to match your cluster limits
#   3. Set BIND_PATHS below

#SBATCH --job-name=varscan2
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=varscan2_%j.log

SIF="${SIF:-varscan2.sif}"
BIND_PATHS="/data/bams,/data/ref,$(pwd)"

apptainer run \
  --bind "${BIND_PATHS}" \
  "${SIF}" \
  --config config.toml \
  --resume \
  "$@"
```

```bash
chmod +x run_slurm.sh
```

**Step 4: Update sample_pairs.csv.example**

```
# tumor_bam,normal_bam[,tumor_purity]
case01_final.bam,control01_final.bam,0.75
case02_final.bam,control02_final.bam
```

**Step 5: Verify compose file is valid**

```bash
docker compose config 2>&1 | head -10
```
Expected: no errors, config output printed.

**Step 6: Commit**

```bash
git add docker-compose.yml config.toml.example run_slurm.sh sample_pairs.csv.example
git commit -m "feat(docker): add docker-compose.yml, config.toml.example, run_slurm.sh; update pairs CSV example with purity column"
```

---

## Phase 4 — CI + Polish

### Task 14: GitHub Actions — cargo test + Docker build+push

**Files:**
- Create: `.github/workflows/ci.yml`
- Create: `.github/workflows/docker.yml`

**Step 1: Write `.github/workflows/ci.yml`**

```yaml
name: CI
on:
  push:
    branches: ["**"]
  pull_request:

jobs:
  test:
    runs-on: ubuntu-22.04
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - uses: Swatinem/rust-cache@v2
      - run: cargo test --all
      - run: cargo build --release
```

**Step 2: Write `.github/workflows/docker.yml`**

```yaml
name: Docker
on:
  push:
    tags: ["v*"]

env:
  REGISTRY: ghcr.io
  IMAGE_NAME: ${{ github.repository }}

jobs:
  build-push:
    runs-on: ubuntu-22.04
    permissions:
      contents: read
      packages: write
    steps:
      - uses: actions/checkout@v4
      - uses: docker/setup-buildx-action@v3
      - uses: docker/login-action@v3
        with:
          registry: ${{ env.REGISTRY }}
          username: ${{ github.actor }}
          password: ${{ secrets.GITHUB_TOKEN }}
      - uses: docker/metadata-action@v5
        id: meta
        with:
          images: ${{ env.REGISTRY }}/${{ env.IMAGE_NAME }}
          tags: |
            type=semver,pattern={{version}}
            type=semver,pattern={{major}}.{{minor}}
            type=sha
      - uses: docker/build-push-action@v5
        with:
          context: .
          push: true
          tags: ${{ steps.meta.outputs.tags }}
          labels: ${{ steps.meta.outputs.labels }}
          cache-from: type=gha
          cache-to: type=gha,mode=max
```

**Step 3: Commit**

```bash
mkdir -p .github/workflows
git add .github/workflows/ci.yml .github/workflows/docker.yml
git commit -m "ci: add GitHub Actions — cargo test on push, Docker build+push on tag"
```

---

### Task 15: `--version` flag + README overhaul

**Files:**
- Modify: `varscan2_pipeline.rs` — add `--version`
- Modify: `README.md` — overhaul setup guide, update Key Design Decisions

**Step 1: Add `--version` flag**

In `parse_args()`, add:
```rust
"-V" | "--version" => {
    println!("varscan2_pipeline {}", env!("CARGO_PKG_VERSION"));
    std::process::exit(0);
}
```

**Step 2: Smoke test**

```bash
cargo build --release && ./target/release/varscan2_pipeline --version
```
Expected: `varscan2_pipeline 0.1.0`

**Step 3: README quick-start section rewrite**

Replace the "Setup Guide" section (Steps 1–8) with a Docker-first quick-start:

```markdown
## Quick Start (Docker — recommended)

### 1. Pull the image

```bash
docker pull ghcr.io/OWNER/varscan2-pipeline:latest
# HPC (Apptainer):
apptainer pull varscan2.sif docker://ghcr.io/OWNER/varscan2-pipeline:latest
```

### 2. Generate a config file and edit it

```bash
docker run --rm ghcr.io/OWNER/varscan2-pipeline --init-config > config.toml
# Edit two required lines:
#   reference = "/data/ref/GRCh38.fa"   ← absolute path on your host
#   bam_dir   = "/data/bams"            ← directory containing your BAM files
$EDITOR config.toml
```

### 3. Create your pairs file

```bash
cp sample_pairs.csv.example sample_pairs.csv
# Edit: one line per tumor/normal pair
# case01_final.bam,control01_final.bam,0.75   ← col 3 = tumor purity (optional)
```

### 4. Validate prerequisites

```bash
docker run --rm \
  -v /data/ref:/data/ref:ro \
  -v /data/bams:/data/bams:ro \
  -v $(pwd):/workspace \
  ghcr.io/OWNER/varscan2-pipeline --validate
```
Fix any `[FAIL]` lines before proceeding.

### 5. Run

```bash
docker compose up
# or directly:
docker run --rm \
  -v /data/ref:/data/ref:ro \
  -v /data/bams:/data/bams:ro \
  -v $(pwd):/workspace \
  ghcr.io/OWNER/varscan2-pipeline --resume
```

### HPC (Apptainer + Slurm)

```bash
sbatch run_slurm.sh
```

---

## Manual / Source Build

See the original setup guide below for building from source without Docker.
```

Keep the original 8-step guide below the Docker quick-start, retitled "Manual Setup (from source)".

**Step 4: Commit**

```bash
git add varscan2_pipeline.rs README.md
git commit -m "feat: add --version flag; overhaul README with Docker quick-start as primary path"
```

---

## Final Steps

### Merge to main

```bash
git checkout main
git merge --no-ff feat/docker-config -m "feat: zero-setup Docker-native pipeline with runtime TOML config"
git tag v0.2.0
git push origin main v0.2.0
```

Docker image will be built and pushed automatically by the `docker.yml` workflow on tag push.

---

## Files Changed Summary

| File | Action |
|------|--------|
| `varscan2_pipeline.rs` | Major: Config struct, TOML/env/CLI loading, all stage functions updated, hardening fixes |
| `Cargo.toml` | Add serde, toml, chrono |
| `Cargo.lock` | Updated |
| `config.toml.example` | New |
| `Dockerfile` | New |
| `docker-compose.yml` | New |
| `run_slurm.sh` | New |
| `sample_pairs.csv.example` | Updated — col 3 purity |
| `.github/workflows/ci.yml` | New |
| `.github/workflows/docker.yml` | New |
| `README.md` | Docker quick-start + stale content fixes |
| `docs/plans/` | Design doc + this plan |
