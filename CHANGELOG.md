# Changelog

## 0.40 - 2026-07-21

- Fixed final-contig variant flushing, sequence-replacement coordinates,
  overlapping replacements, bounds checks, multi-allelic/MNP handling, and
  two-VCF reconciliation/filtering.
- Added explicit metagenomic VCF `FILTER` policies and stronger bacterial
  defaults, including minor-allele support and anomalously high-depth filters.
- Made repeated reference/VCF mismatches fatal while allowing a small,
  configurable number of reported mismatches.
- Hardened VCF, GFF, FASTA-header, depth-interval, and command-line parsing.
- Added bacterial genetic-code/phase-aware translation and restricted GFF
  translation to appropriate, functional non-pseudo feature types.
- Corrected coverage/`N` reporting, per-gene accounting, depth overflow, and
  mutation-matrix ownership; depth masks now use a documented multi-source
  union while call-depth means remain source-specific.
- Matched FASTA/GFF/VCF sequence identifiers safely, protected input/output
  paths from accidental truncation, and added end-to-end regression coverage.
- Made `gzstream` optional so an uncompressed Linux build works from a clean
  checkout; updated `-h` to match the implemented interface.
- Removed the generic `-minFS` threshold because caller-specific `FS` semantics
  are not safe to interpret uniformly.
