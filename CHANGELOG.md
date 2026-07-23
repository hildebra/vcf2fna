# Changelog

## 0.42 - 2026-07-23

- Accepted valid subnormal scientific-notation values in VCF numeric fields,
  including extremely small `FS` and `PV4`-style p-values.
- Limited repeated warnings to five examples per category, followed by one
  suppression notice, without weakening mismatch counting or abort limits.
- Disabled relative high-depth filtering by default. Metagenomic coverage is
  intrinsically uneven, so passing major ALT alleles now replace de novo
  assembly bases unless an explicit high-depth scale or another evidence filter
  rejects them.

## 0.41 - 2026-07-22

- Re-enabled automatic gzip input/output for Linux builds using a bundled
  zlib-backed C++ stream and Linux-aware Makefile/test compilation.
- Retained explicit opt-in support on other operating systems and an explicit
  `VCF2FNA_DISABLE_GZIP` escape hatch for Linux builds without zlib.

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
