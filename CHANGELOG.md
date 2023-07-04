# Changelog: nf-cavalier
<!--- https://keepachangelog.com/en/1.0.0/ --->

## [Unreleased]
### Changed
- Cache dir now uses package version in case of changes in formatting of cache files between versions
- PanelApp lists stored in cache as well if a specific list version is requested
- Add support for Structural variants
- Better messaging about failed HTTP requests
- Removed reliance on HPO API for disease names by using hpo.annotation build artifacts instead. HPO API seems to be missing some OMIM diseases.