# VarScan2 Pipeline Redesign — Design Document

**Date:** 2026-04-30
**Approach:** A — Rust orchestrator + TOML config + Docker/Apptainer
**Goal:** Zero-friction, production-ready, Docker-native VarScan2 pipeline

---

## Objectives

1. **Zero-setup:** `git clone` + edit 2 config lines + `docker compose up` → running analysis
2. **Runtime config:** replace compile-time constants with TOML + env var + CLI override
3. **Containerization:** Docker image bundles all tools; Apptainer for HPC
4. **Hardening:** fix all FAIL/WARN items from varscan_review_prompt.md audit
5. **Per-sample purity:** pairs CSV column 3 overrides global tumor_purity
6. **Pre-flight validation:** `--validate` flag catches missing files before wasting compute

---

## Section 1: Config System

### Priority chain (highest to lowest)
```
CLI flag > VARSCAN_* env var > config.toml > compiled defaults
```

### config.toml structure

```toml
[paths]
reference    = "/data/ref/GRCh38.fa"   # required
bam_dir      = "."
bam_suffix   = "_final.bam"
pairs_file   = "sample_pairs.csv"
pairs_suffix = "_final.bam"
target_bed   = ""
software_dir = "software"
scripts_dir  = "scripts"

[somatic]
min_coverage        = 20
min_coverage_normal = 10
min_coverage_tumor  = 20
min_base_qual       = 20
min_var_freq        = 0.10
min_freq_for_hom    = 0.75
normal_purity       = 1.0
tumor_purity        = 1.0
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
```

### Env var mapping
Flat `VARSCAN_` prefix, section-unaware. Examples:
- `VARSCAN_REFERENCE=/data/ref/GRCh38.fa`
- `VARSCAN_MIN_COVERAGE=20`
- `VARSCAN_TUMOR_PURITY=0.65`

### Per-sample tumor_purity
Pairs CSV gains optional column 3:
```
case01_final.bam,ctrl01_final.bam,0.65
case02_final.bam,ctrl02_final.bam       # uses [somatic] tumor_purity
```

### New Cargo.toml dependencies
- `serde` + `serde_derive`
- `toml`

### New flags
| Flag | Behavior |
|------|----------|
| `--config <path>` | Load config from path (default: `./config.toml`) |
| `--init-config` | Write `config.toml.example` to cwd, exit 0 |
| `--validate` | Pre-flight check all paths/tools/BAI indexes, exit without running |
| `--version` | Print semver from Cargo.toml, exit 0 |

---

## Section 2: Docker / Apptainer Containerization

### Multi-stage Dockerfile

**Stage 1 — builder:** `rust:1.78-slim`
- `cargo build --release`
- Output: `/usr/local/bin/varscan2_pipeline`

**Stage 2 — runtime:** `ubuntu:22.04`
- samtools 1.20 (bioconda or compiled)
- bam-readcount 0.8 (compiled from source)
- openjdk-17-jre-headless
- perl + cpanm Statistics::Descriptive
- VarScan.v2.3.9.jar → `/opt/varscan/VarScan.v2.3.9.jar` (SHA256 verified at build time)
- fpfilter.pl → `/opt/varscan/fpfilter.pl`
- Default config paths: `software_dir=/opt/varscan`, `scripts_dir=/opt/varscan`

**Image target size:** ~600 MB

**Version pinning:** all tool versions pinned; VarScan jar SHA256 verified at build time.

### docker-compose.yml (local)

```yaml
services:
  varscan:
    image: ghcr.io/<owner>/varscan2-pipeline:latest
    volumes:
      - ./bams:/data/bams
      - ./ref:/data/ref
      - ./config.toml:/workspace/config.toml
      - ./sample_pairs.csv:/workspace/sample_pairs.csv
      - ./results:/workspace
    working_dir: /workspace
    command: ["--resume"]
```

### Zero-setup local flow
```bash
git clone <repo> && cd varscan
docker run --rm ghcr.io/.../varscan2-pipeline --init-config > config.toml
# edit: reference path + bam_dir
docker compose run varscan --validate
docker compose run varscan --resume
```

### HPC (Apptainer)
```bash
apptainer pull varscan2.sif docker://ghcr.io/<owner>/varscan2-pipeline:latest
apptainer run --bind /data/bams,/data/ref,$(pwd) varscan2.sif \
  --config config.toml --resume
```

### Slurm wrapper (run_slurm.sh — ships in repo)
```bash
#!/usr/bin/env bash
#SBATCH --cpus-per-task=8 --mem=32G --time=48:00:00
apptainer run --bind /data/bams,/data/ref,$(pwd) \
  varscan2.sif --config config.toml --resume "$@"
```

### GitHub Actions CI
- Build + push Docker image on tag push → `ghcr.io/<owner>/varscan2-pipeline:<tag>`
- Run `cargo test` on every push/PR

---

## Section 3: Code Hardening

### FAIL — must fix

| ID | Location | Issue | Fix |
|----|----------|-------|-----|
| F1 | README:1112 | Key Design Decisions lists `VarScan somatic --min-map-qual 10` — stale after last commit | Remove from README |
| F2 | `wait_all` | On child failure, remaining children are not killed — leaves zombie processes | On first non-zero exit, kill remaining children, then return error |
| F3 | `run_fpfilter` | Audit VCF lookup pattern: `base + ".vcf"` — verify Germline/LOH readcount → VCF filename alignment | Confirm or fix pattern |
| F4 | Config validation | `reference` path emptiness check only happens at startup; with runtime config, validation must be a dedicated pass before `setup_directories` | Add `validate_config()` called before any stage runs |

### WARN — should fix

| ID | Issue | Fix |
|----|-------|-----|
| W1 | `now_string()` shells out to `date` — fragile in minimal containers | Replace with `chrono` crate or manual `std::time` formatting |
| W2 | No warning when `tumor_purity=1.0` (default) | Emit `[WARN] tumor_purity=1.0 for <sample> — set col 3 in pairs CSV if purity known` |
| W3 | Stage 4 guard skipped when `from_stage >= 7` — `.hc.vcf` files may not exist | Add: if `from_stage >= 7`, verify `.hc.vcf` files present before proceeding |
| W4 | Custom `glob_match` — covers more than needed; direct suffix checks suffice for all actual patterns | Replace with `str::ends_with` checks; remove glob machinery |
| W5 | `hash_manifest` uses mtime — breaks on `rsync --archive` copy | Document limitation prominently; optionally add `--content-hash` flag to use SHA256 instead of mtime |
| W6 | `bam_dir` hardcoded to `cwd` | Add `[paths] bam_dir` to config (done in Section 1) |

### INFO

- `clean_pair_field` strips spaces but not tabs — pairs files with tab-separated fields silently fail; add tab stripping
- `VCF_SAMPLE_LIST` line count not validated against pairs count — add check

---

## Section 4: New Capabilities & Zero-Setup UX

### `--validate` pre-flight output format
```
[OK]   reference: /data/ref/GRCh38.fa (indexed)
[OK]   pairs file: sample_pairs.csv (3 pairs)
[OK]   samtools 1.20 on PATH
[OK]   java 17 on PATH
[OK]   bam-readcount 0.8 on PATH
[OK]   VarScan.v2.3.9.jar present
[OK]   fpfilter.pl present
[WARN] tumor_purity=1.0 for all pairs — set col 3 in pairs CSV if purity known
[WARN] target_bed not set — full-genome mpileup (WGS mode)
[FAIL] case03_final.bam.bai: not found
```
Exits non-zero if any `[FAIL]` present.

### README overhaul
Replace 8-step manual setup guide with 4 steps:
1. `docker pull` (or `apptainer pull`)
2. `--init-config` → edit `reference` + `bam_dir`
3. `--validate`
4. run

Keep parameter reference and Key Design Decisions sections (fix stale content).

---

## Implementation Phases

| Phase | Scope | Deliverable |
|-------|-------|-------------|
| 1 | Config system | TOML + serde, env var override, CLI `--config`/`--init-config`, per-sample purity, `validate_config()` |
| 2 | Hardening | Fix F1–F4, W1–W6, INFO items; add `--validate` flag |
| 3 | Docker | Dockerfile (multi-stage), docker-compose.yml, run_slurm.sh, config.toml.example |
| 4 | CI + polish | GitHub Actions build+push, `--version`, README overhaul |

---

## Files Changed / Created

| File | Action |
|------|--------|
| `varscan2_pipeline.rs` | Major: config system, hardening fixes, new flags |
| `Cargo.toml` | Add serde, toml, chrono deps |
| `config.toml.example` | New |
| `Dockerfile` | New |
| `docker-compose.yml` | New |
| `run_slurm.sh` | New |
| `.github/workflows/docker.yml` | New |
| `README.md` | Overhaul setup guide; fix stale --min-map-qual reference |
| `sample_pairs.csv.example` | Update to show col 3 tumor_purity example |
