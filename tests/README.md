# Regression tests

Run the self-contained end-to-end suite from Linux or WSL:

```sh
python3 tests/run_regression_tests.py
```

The runner builds the current source with `g++ -std=c++17` (and `-lz` on Linux) in a temporary
directory, generates all FASTA/VCF/depth/GFF3 fixtures at runtime, and removes
them on completion. Override the compiler with `CXX` or `--cxx`. To test an
existing Linux executable without rebuilding:

```sh
python3 tests/run_regression_tests.py --binary ./vcf2fna
```

The cases cover EOF flushing, VCF-anchored SNP/indel placement, FILTER modes,
multiallelic allele-depth selection, the REF mismatch limit, FASTA/GFF3 seqid
matching, coding-feature selection and multipart/phase/strand handling,
bacterial translation tables 4/11/25, reference `N` versus depth masking, and
one/two-caller consensus agreement and conflict behavior. Additional edge cases
cover CRLF VCFs, symbolic ALT placeholders, duplicate decomposed calls, filtered
deletions, SNP/indel cross-caller conflicts, and indels on CDS boundaries.
They also exercise IDV/IMF-only indel evidence, unsupported translation-table
handling by output type, strict numeric GFF parsing, ambiguous ALT accounting,
and SNP header-coordinate shifts after upstream indels.
The suite additionally protects cross-caller REF validation, independent depth
summaries for two inputs, 64-bit gene-depth accumulation, `pseudo=true` CDS
exclusion, and rejection of overlapping multipart CDS segments.
Final edge coverage includes exact technical-FILTER exceptions, inherited
pseudogene status, MNP application and gene-boundary behavior, union masking
for overlapping replacements/MNP-plus-SNP conflicts, invalid-indel proximity
isolation, and `-minCallDepth 0` with sparse or wholly omitted intervals. It
also covers multipart delins coordinate merging, Sequence Ontology CDS/ORF
handling, strict missing-QUAL behavior, and output-path collision protection.
Linux runs additionally verify transparent gzip input and output.
The closing safety cases cover a delins crossing multipart CDS segments,
Sequence Ontology CDS-versus-ORF phase rules, configurable handling of missing
VCF QUAL, and preflight rejection of duplicate or input-colliding outputs.
