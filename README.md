# VCF2fna

VCF2fna builds bacterial consensus contigs and annotated gene sequences from a
metagenomic reference, one or two VCF call sets, matching depth tracks, and a
GFF annotation. Its filtering model is intentionally haploid/community-aware:
VCF genotype (`GT`) is not used to decide the consensus allele.

Run `vcf2fna -h` for the complete command-line reference. A minimal run is:

```text
vcf2fna -ref reference.fa -inVCF calls.vcf -depthF depth.bed \
  -gff genes.gff -oCtg consensus.fa
```

One or two comma-separated VCFs may be supplied. `-depthF` and
`-minCallDepth` must have the same cardinality as `-inVCF`; minimum depth
defaults to 5 per input. Optional `-seqPlatform` labels must likewise have one
value per VCF, but are informational: platform names alone are too coarse to
justify silently changing thresholds, so all actual filter settings remain
explicit command-line parameters.

Depth files use bedGraph/BED coordinates: `contig start end depth`, with
zero-based, half-open intervals. Within each file, intervals must be grouped by
contig, sorted, and non-overlapping; omitted bases have depth zero. A base is
resolvable when any supplied depth source meets its matching
`-minCallDepth`. Variant filters nevertheless compare each VCF call with the
mean depth from its own matching source, so one technology cannot inflate the
other's expected depth. A minimum depth of 0 disables coverage masking for that
source, including at intervals omitted from a sparse depth file.

## VCF FILTER policy

`-vcfFilterPolicy technical` is the default. It rejects caller FILTER codes
unless they clearly identify a diploid, germline/somatic, or population-model
criterion that is not valid evidence against a pooled bacterial call. Unknown
named failures are rejected conservatively. `all` rejects every
non-`PASS`/non-`.` FILTER code. `ignore` restores the legacy behavior of
ignoring the FILTER column.
Independent depth, allele-support, quality, strand/mapping-bias, and indel
proximity checks still apply under every policy.

The default credible-minor-allele thresholds are four supporting reads and 5%
frequency (`-minAltReads 4 -minAltFreq 0.05`); indels and other complex
replacements retain a conservative 10% frequency floor because their alignment
error background is higher. Missing VCF `QUAL` fails whenever `-minCallQual`
is nonzero. Sites far above their contig's
mean depth are also filtered by default (`-maxDepthFilterScale 3`); set that
option to 0 to disable the high-depth check. Reference/VCF allele mismatches are
reported, and processing aborts after more than ten by default
(`-maxRefMismatches`).

Sequence-resolved multi-nucleotide/block substitutions use the same
conservative replacement path as indels and are included in the indel/complex
statistics. Overlapping replacements at different VCF anchors are not combined
into a guessed haplotype; their complete reference union is masked instead.

Protein output supports NCBI translation tables 1, 4, 11, and 25. If a GFF
feature contains `transl_except`, amino-acid output stops with an explicit
error rather than silently mistranslating a selenocysteine or pyrrolysine;
contig and nucleotide-only output are unaffected.

GFF3 selection retains functional CDS/ORF features (including their Sequence
Ontology forms) and functional gene fallbacks, while excluding pseudogenes and
their descendants. GFF/VCF sequence IDs are matched to the first
whitespace-delimited token of each FASTA header. Full contig headers are
preserved in contig output; gene FASTA uses canonical `seqid_ordinal` headers.
All output paths are checked before processing and must be distinct from one
another and from every input path.

## Building

The source requires a C++17 compiler and has no mandatory compression library:

```text
g++ -std=c++17 -O2 -I. vcf2fna.cpp options.cpp fasta.cpp vcf.cpp -o vcf2fna
```

This default build reads and writes uncompressed files. Transparent `.gz`
support is optional because `gzstream` is not vendored in this repository. To
enable it, provide `gzstream.h` plus its implementation/library and zlib, define
`VCF2FNA_USE_GZSTREAM`, and add the corresponding include/linker flags. A build
without that macro fails with an explicit message only when a `.gz` path is
requested.
