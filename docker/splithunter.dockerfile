FROM rust:1.82-slim-bookworm AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        pkg-config \
        clang \
        cmake \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        python3 \
        python3-pip \
        python3-dev \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

RUN pip install --no-cache-dir --upgrade pip \
 && pip install --no-cache-dir "maturin>=1.5,<2.0"

WORKDIR /src
COPY . /src

RUN maturin build --release --out /wheels

FROM python:3.11-slim-bookworm AS runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
        libcurl4 \
        libbz2-1.0 \
        liblzma5 \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir /wheels/*.whl boto3 awscli

ENTRYPOINT ["splithunter_run"]
