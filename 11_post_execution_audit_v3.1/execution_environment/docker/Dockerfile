# Reference Linux container for BRD v3.1 external tools (Compleasm v0.2.10 + Broccoli v2.0.1)
# Route B execution status on host: NOT_EXECUTED_DOCKER_ABSENT (docker binary unavailable on Windows 11 host)
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
WORKDIR /opt/brd

RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates wget bzip2 git build-essential zlib1g-dev \
    python3 python3-pip hmmer fasttree ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

# Pinned Biopython & pandas
RUN pip3 install --no-cache-dir --break-system-packages biopython==1.86 pandas==2.3.3

# Pinned Compleasm v0.2.10 x64-linux release (includes miniprot & hmmsearch Linux x64 binaries)
RUN wget -q https://github.com/huangnengCSU/compleasm/releases/download/v0.2.10/compleasm-0.2.10_x64-linux.tar.bz2 \
    && tar -xjf compleasm-0.2.10_x64-linux.tar.bz2 \
    && rm compleasm-0.2.10_x64-linux.tar.bz2
ENV PATH="/opt/brd/compleasm_kit:${PATH}"

# Pinned Broccoli v2.0.1 (Rust) + DIAMOND 2.1+
RUN wget -q https://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz \
    && tar -xzf diamond-linux64.tar.gz -C /usr/local/bin && rm diamond-linux64.tar.gz \
    && wget -q https://sh.rustup.rs -O /tmp/rustup-init.sh \
    && sh /tmp/rustup-init.sh -y --profile minimal \
    && export PATH="/root/.cargo/bin:${PATH}" \
    && git clone --depth 1 --branch v2.0.1 https://github.com/rderelle/Broccoli.git /opt/brd/Broccoli \
    && cd /opt/brd/Broccoli && cargo build --release \
    && cp target/release/broccoli /usr/local/bin/broccoli \
    && rm -rf /root/.cargo /root/.rustup

WORKDIR /workspace
CMD ["/bin/bash"]
