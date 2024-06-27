# Changelog: cavalier
<!--- https://keepachangelog.com/en/1.0.0/ --->

## [24.06.0]
### Fixed
- HPO API updated to current version

### Added
- Gene4Epilepsy gene lists: e.g. `get_gene_list('G4E:ALL')`
- HGNC gene lists: e.g.  `get_gene_list('HGNC:protein-coding')`
- Caching of all files the rely on web resources
- Option `read_only_cache_dir` to allow shared cavalier cache directory: e.g. `set_cavalier_opt(read_only_cache_dir = '/PATH/TO/DIR')`
- All URLs to remote resources as options that can be changed in case they move
- Roxygen function documentation to many functions

### Changed
- improve compound heterozygous detection in `inheritance.R`
- moved Dockerfile to this repository from nf-cavalier to streamline building docker container

### Removed
- RVIS scores as outdated and certificate has expired on download page