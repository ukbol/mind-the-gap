# Runtime image for mind-the-gap automated runs on GitHub Actions Larger Runners.
# Pre-installs Python deps, VSEARCH, and NCBI EDirect so each scheduled run
# starts in seconds instead of paying ~5 min of apt/conda solves.
FROM mambaorg/micromamba:1.5.8

ARG VSEARCH_VERSION=2.28.1
ARG ENTREZ_DIRECT_VERSION=22.4

USER root

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        git \
        tini \
    && rm -rf /var/lib/apt/lists/*

RUN micromamba install -y -n base -c conda-forge -c bioconda \
        python=3.11 \
        "vsearch=${VSEARCH_VERSION}" \
        "entrez-direct=${ENTREZ_DIRECT_VERSION}" \
    && micromamba clean -ay

ENV PATH=/opt/conda/bin:$PATH

COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

WORKDIR /work

ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["bash"]
