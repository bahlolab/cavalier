FROM rocker/verse:4.5.2

# use posit package manager for fast binary installs
ENV RSPM="https://packagemanager.posit.co/cran/__linux__/jammy/latest"

# Install BiocManager and GenomicRanges explicitly
RUN R -q -e "install.packages('BiocManager'); BiocManager::install('GenomicRanges')"

# Copy package source
COPY . /tmp/cavalier

# Install cavalier
RUN R -q -e "remotes::install_local('/tmp/cavalier', dependencies = TRUE, upgrade = 'never', build_vignettes = FALSE)"

# Initialise cavalier cache
RUN mkdir /cavalier_cache && 
  R -q -e "cavalier::set_cavalier_opt(cache_dir = '/cavalier_cache'); cavalier::build_caches(PanelApp = FALSE)"