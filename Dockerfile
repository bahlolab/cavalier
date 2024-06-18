FROM mambaorg/micromamba:1.5.8-noble

LABEL \
  author="Jacob Munro" \
  description="Container for Cavalier" \
  maintainer="Bahlo Lab"

# install os deps
USER root
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    procps \
    xvfb \
    xauth \
  && apt-get clean -y \
  && rm -rf /var/lib/apt/lists/*

# install env with micromamba
COPY environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml \ 
  && micromamba clean --all --yes

# create igv.sh
RUN cp /opt/conda/bin/igv /opt/conda/bin/igv.sh \
    && echo "PATH=/opt/conda/bin:${PATH}" >> \
    /opt/conda/lib/R/etc/Renviron

# Install cavalier R package and build caches
COPY . /tmp/cavalier
RUN /opt/conda/bin/R --slave --vanilla -e "\
    devtools::install(pkg = '/tmp/cavalier', force = TRUE, upgrade = 'never'); \
    cavalier::set_cavalier_opt(cache_dir = NULL); \
    cavalier::build_caches() \
    "

ENV PATH="/opt/conda/bin:${PATH}" \
    TZ=Etc/UTC \
    R_HOME=/opt/conda/lib/R/ \
    R_ENVIRON=/opt/conda/lib/R/etc/Renviron \
    R_LIBS_USER=/opt/conda/lib/R/site-library