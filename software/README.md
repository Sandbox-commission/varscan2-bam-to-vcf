# software/

This directory must contain `VarScan.v2.3.9.jar` before running the pipeline.

## Obtaining VarScan.v2.3.9.jar

Download from the official GitHub release:

```
https://github.com/dkoboldt/varscan/releases/tag/2.3.9
```

Direct download:

```bash
wget -O software/VarScan.v2.3.9.jar \
  https://github.com/dkoboldt/varscan/releases/download/2.3.9/VarScan.v2.3.9.jar
```

Expected SHA256:

```
e07a7d11d5f2e0f93d72b2293b9c5e879ff62742f6fe5a1bc0c16abfe65dd5e6
```

Verify before running:

```bash
sha256sum software/VarScan.v2.3.9.jar
```

## Requirements

- Java JRE 8 or later must be on `PATH` (`java -version` to check).
