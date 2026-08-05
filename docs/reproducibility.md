# Reproducibility notes

## What is reproducible from this repository

The repository records analysis code, script ordering within several modules, statistical thresholds embedded in those scripts, and a synthetic format test. It does not include controlled inputs, complete production configuration, reference assets, or an author-verified environment lock.

## Required run record

For each production run, record:

- Git commit and any local, reviewed patches.
- Operating system, scheduler, CPU architecture, thread count, and memory.
- Python, R, package, and external-tool versions.
- Reference genome build, annotation/database releases, and SHA256 checksums.
- Input manifests using de-identified IDs and separate controlled linkage storage.
- Command lines, environment variables, random seeds, and exit statuses.
- Output checksums and a non-sensitive quality-control summary.

## Reference builds

The code references T2T-CHM13v2.0/GCF_009914755.1 in several assembly, genotyping, and annotation paths, while mapping and some downstream analyses use GRCh38/hg38 labels. These builds are not interchangeable. Confirm every input/output coordinate system and any liftover step.

## Randomness

Some analysis and plotting scripts may set seeds locally; others call tools whose parallel execution can affect reproducibility. Do not add or change seeds without confirming the study method. Capture each explicit seed and tool-specific deterministic setting in the run record.

## Resource effects

Thread count, memory, scheduler behavior, and library versions can change performance and, for some tools, output ordering or numerical details. Retained defaults are source-derived when identifiable, but site-specific resource directives were removed or parameterized where they exposed private infrastructure.

## Validation boundary

Static parsing does not establish scientific reproducibility. The release requires author-reviewed end-to-end testing on authorized data in an approved environment before publication.
