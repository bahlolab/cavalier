FROM rocker/verse:4.5.2

# use posit package manager
ENV RSPM="https://packagemanager.posit.co/cran/__linux__/jammy/latest"

# Install CRAN packages (binary via PPM where available)
# RUN install2.r --error --skipinstalled \
#     multidplyr \
#     broom \
#     janitor \
#     patchwork \
#     DT \
#     ggrepel \
#     pheatmap

COPY . /tmp/cavalier

RUN R -q -e "remotes::install_local('/tmp/cavalier', dependencies = TRUE, upgrade = 'never', build_vignettes = FALSE)"