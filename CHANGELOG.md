# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Fixed
  - Missing the ERCC spike-in for GLDS-273

## [0.3.1] - 2021-06-16
### Fixed
  - Added in "spike-in protocol" as a protocol type that indicates ERCC spike-in

## [0.3.0] - 2021-06-16
### Changed
  - Proto run sheet is no longer deleted by default.
  - Final run sheet now uses the 'remove' directive to drop unneeded columns.

### Removed
  - Read length extraction, prior knowledge of raw read length is no longer required for processing.

### Fixed
  - Duplicate sample_name in final run sheet removed.

## [0.2.1] - 2021-04-07
### Added
  - Added additional readlength parsing keys.
    - Reason: GLDS-373 metadata format has read length key'ed differently.

## [0.2.0] - 2021-04-02
### Added
  - Added organism and readlength parsing to runsheet.
    - Tested on GLDS-48, 83, 104, 192, 194
    - Note: GLDS-192 fails but this seems to be due to missing read length metadata.


## [0.1.0] - 2021-04-02
  - First pre-release tracked
