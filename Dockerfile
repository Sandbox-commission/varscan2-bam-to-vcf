# ── Stage 1: Rust builder ─────────────────────────────────────────────────────
FROM rust:1.78-slim AS builder

WORKDIR /build
COPY Cargo.toml Cargo.lock ./
COPY varscan2_pipeline.rs ./
RUN cargo build --release

# ── Stage 2: bam-readcount builder ───────────────────────────────────────────
FROM ubuntu:22.04 AS brc-builder

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake git libncurses5-dev libbz2-dev zlib1g-dev \
    liblzma-dev libcurl4-openssl-dev libssl-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 --branch 0.8.0 \
    https://github.com/genome/bam-readcount.git /src/bam-readcount \
    && cmake -S /src/bam-readcount -B /src/brc-build \
       -DCMAKE_BUILD_TYPE=Release \
    && cmake --build /src/brc-build --parallel $(nproc) \
    && install -m 755 /src/brc-build/bin/bam-readcount /usr/local/bin/bam-readcount

# ── Stage 3: runtime image ───────────────────────────────────────────────────
FROM ubuntu:22.04

ARG VARSCAN_SHA256=e67c75b69cb22ac56618d63e414a0d0c10787f35f50e972636f9cfad0e46a298

LABEL org.opencontainers.image.title="varscan2-pipeline" \
      org.opencontainers.image.description="VarScan2 somatic variant and CNV pipeline" \
      org.opencontainers.image.source="https://github.com/Sandbox-commission/varscan"

RUN apt-get update && apt-get install -y --no-install-recommends \
    samtools \
    openjdk-17-jre-headless \
    perl \
    cpanminus \
    wget \
    ca-certificates \
    libncurses5 \
    libbz2-1.0 \
    zlib1g \
    liblzma5 \
    libcurl4 \
    && cpanm --quiet Statistics::Descriptive \
    && rm -rf /var/lib/apt/lists/* /root/.cpanm

# Copy pipeline binary and bam-readcount
COPY --from=builder  /build/target/release/varscan2_pipeline /usr/local/bin/varscan2_pipeline
COPY --from=brc-builder /usr/local/bin/bam-readcount          /usr/local/bin/bam-readcount

# Download and verify VarScan jar
RUN mkdir -p /opt/varscan \
    && wget -q -O /opt/varscan/VarScan.v2.3.9.jar \
       https://github.com/dkoboldt/varscan/releases/download/2.3.9/VarScan.v2.3.9.jar \
    && echo "${VARSCAN_SHA256}  /opt/varscan/VarScan.v2.3.9.jar" | sha256sum -c -

# Download fpfilter.pl
RUN wget -q -O /opt/varscan/fpfilter.pl \
    https://raw.githubusercontent.com/genome/fpfilter-tool/master/fpfilter.pl \
    && chmod 755 /opt/varscan/fpfilter.pl

WORKDIR /workspace

# Default config points to bundled tools; user mounts BAMs + ref + config
ENV VARSCAN_SOFTWARE_DIR=/opt/varscan \
    VARSCAN_SCRIPTS_DIR=/opt/varscan

ENTRYPOINT ["/usr/local/bin/varscan2_pipeline"]
CMD ["--help"]
