# somascan-qc : QC + filtering for ANML-normalized SomaScan ADAT files
#   docker build -t somascan-qc:1.1.0 .
#   docker run --rm -v "$PWD":/data somascan-qc:1.1.0 run -i /data/input.adat -o /data/out
FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="somascan-qc" \
      org.opencontainers.image.version="1.1.0" \
      org.opencontainers.image.description="Reproducible QC pipeline for SomaScan ADAT data"

USER root
WORKDIR /opt/somascan-qc
COPY environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

COPY pyproject.toml README.md LICENSE ./
COPY somascan_qc ./somascan_qc
COPY config ./config
RUN /opt/conda/bin/python -m pip install --no-deps --no-cache-dir . && \
    mkdir -p /data && chmod 777 /data

ENV MPLBACKEND=Agg MPLCONFIGDIR=/tmp/mpl PYTHONDONTWRITEBYTECODE=1
WORKDIR /data
ENTRYPOINT ["/opt/conda/bin/somascan-qc"]
# Smoke-test the image:  docker run --rm -v "$PWD":/data somascan-qc:1.1.0 demo -o /data
CMD ["--help"]
