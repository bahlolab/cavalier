FROM rocker/verse:4.5.2

# use posit package manager for fast binary installs
ENV RSPM="https://packagemanager.posit.co/cran/__linux__/jammy/latest"

# 1. Install BiocManager and GenomicRanges explicitly
RUN R -q -e "install.packages('BiocManager'); BiocManager::install('GenomicRanges')"

# 2. Copy package source
COPY . /tmp/cavalier

# 3. Install cavalier
RUN R -q -e "remotes::install_local('/tmp/cavalier', dependencies = TRUE, upgrade = 'never', build_vignettes = FALSE)"