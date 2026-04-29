# scripts/

This directory must contain `fpfilter.pl` before running the pipeline.

## Obtaining fpfilter.pl

`fpfilter.pl` is part of the VarScan2 source distribution and is also maintained
separately at:

```
https://github.com/genome/fpfilter-tool
```

Download:

```bash
wget -O scripts/fpfilter.pl \
  https://raw.githubusercontent.com/genome/fpfilter-tool/master/fpfilter.pl
chmod +x scripts/fpfilter.pl
```

## Requirements

- Perl 5.10 or later must be on `PATH` (`perl --version` to check).
- The following Perl modules are required:
  - `Statistics::Descriptive` (install via `cpanm Statistics::Descriptive`)
