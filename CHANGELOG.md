# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
  - Microarray Runsheet Staging from GLDS
### Changed
  - Converted RNASeq staging from isa object model to pandas dataframe backend.  No anticipated changes on user-end

## [0.3.0] - 2021-05-18
### Added
  - Microarray Proto Runsheet Staging from GLDS
  - CLI arg to specify staging an incomplete runsheet when missing metadata columns are present

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
