#!/usr/bin/env python3
"""End-to-end regression tests for VCF2fna.

The tests generate their own FASTA, VCF, depth, and GFF3 inputs.  By default
the current source tree is compiled once with a C++17 compiler; --binary can
be used to exercise an already-built executable instead.
"""

from __future__ import annotations

import argparse
import gzip
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
import unittest


VCF_HEADER = (
    "##fileformat=VCFv4.3\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def vcf_call(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    ad: tuple[int, ...] = (1, 19),
    filt: str = "PASS",
    qual: int | str = 100,
) -> str:
    """Return one pooled-sample VCF record with allele-depth evidence."""
    depth = sum(ad)
    depths = ",".join(str(value) for value in ad)
    return (
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filt}\tDP={depth}"
        f"\tAD:DP\t{depths}:{depth}"
    )


def vcf_info_call(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    info: str,
    filt: str = "PASS",
    qual: int = 100,
) -> str:
    """Return a VCF record whose allele evidence is carried only in INFO."""
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filt}\t{info}\t.\t."


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Read a FASTA into first-token id -> (complete header, sequence)."""
    records: dict[str, tuple[str, str]] = {}
    header: str | None = None
    chunks: list[str] = []

    def store() -> None:
        if header is None:
            return
        sequence_id = header.split()[0]
        if sequence_id in records:
            raise AssertionError(f"duplicate FASTA id {sequence_id!r} in {path}")
        records[sequence_id] = (header, "".join(chunks))

    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                store()
                header = line[1:]
                chunks = []
            else:
                if header is None:
                    raise AssertionError(f"sequence before FASTA header in {path}")
                chunks.append(line)
    store()
    return records


class VCF2fnaRegressionTests(unittest.TestCase):
    repo = Path(__file__).resolve().parents[1]
    binary_override: Path | None = None
    cxx = os.environ.get("CXX", "g++")
    extra_cxxflags = os.environ.get("CXXFLAGS", "")

    @classmethod
    def setUpClass(cls) -> None:
        cls._workspace = tempfile.TemporaryDirectory(prefix="vcf2fna-regression-")
        cls.workspace = Path(cls._workspace.name)
        if cls.binary_override is not None:
            cls.binary = cls.binary_override.resolve()
            if not cls.binary.is_file():
                raise RuntimeError(f"test binary does not exist: {cls.binary}")
            return

        cls.binary = cls.workspace / "vcf2fna"
        compiler = shlex.split(cls.cxx)
        if not compiler:
            raise RuntimeError("CXX resolved to an empty command")
        command = [
            *compiler,
            "-std=c++17",
            "-O0",
            "-g",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            *shlex.split(cls.extra_cxxflags),
            "-I",
            str(cls.repo),
            str(cls.repo / "vcf2fna.cpp"),
            str(cls.repo / "options.cpp"),
            str(cls.repo / "fasta.cpp"),
            str(cls.repo / "vcf.cpp"),
            *(["-lz"] if sys.platform.startswith("linux") else []),
            "-o",
            str(cls.binary),
        ]
        try:
            build = subprocess.run(
                command,
                cwd=cls.repo,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=120,
                check=False,
            )
        except FileNotFoundError as error:
            raise RuntimeError(
                f"C++ compiler not found ({compiler[0]!r}); set CXX or pass --binary"
            ) from error
        if build.returncode != 0:
            raise RuntimeError(
                "C++17 build failed\ncommand: "
                + shlex.join(command)
                + "\nstdout:\n"
                + build.stdout
                + "\nstderr:\n"
                + build.stderr
            )

    @classmethod
    def tearDownClass(cls) -> None:
        cls._workspace.cleanup()

    def make_fixture(
        self,
        name: str,
        reference: list[tuple[str, str]],
        calls1: list[str] | None = None,
        *,
        gff_rows: list[str] | None = None,
        depth1: list[str] | None = None,
        calls2: list[str] | None = None,
        depth2: list[str] | None = None,
    ) -> dict[str, object]:
        case = self.workspace / name
        case.mkdir()
        ref_path = case / "reference.fa"
        ref_path.write_text(
            "".join(f">{header}\n{sequence}\n" for header, sequence in reference),
            encoding="utf-8",
        )
        ids_and_lengths = [(header.split()[0], len(sequence)) for header, sequence in reference]
        default_depth = [f"{seqid}\t0\t{length}\t20" for seqid, length in ids_and_lengths]

        gff_path = case / "features.gff3"
        gff_path.write_text(
            "##gff-version 3\n" + "".join(f"{row}\n" for row in (gff_rows or [])),
            encoding="utf-8",
        )
        vcf1_path = case / "calls-1.vcf"
        vcf1_path.write_text(
            VCF_HEADER + "".join(f"{row}\n" for row in (calls1 or [])),
            encoding="utf-8",
        )
        depth1_path = case / "depth-1.bed"
        depth1_path.write_text(
            "".join(f"{row}\n" for row in (depth1 or default_depth)),
            encoding="utf-8",
        )

        result: dict[str, object] = {
            "dir": case,
            "ref": ref_path,
            "gff": gff_path,
            "vcfs": [vcf1_path],
            "depths": [depth1_path],
        }
        if calls2 is not None:
            vcf2_path = case / "calls-2.vcf"
            vcf2_path.write_text(
                VCF_HEADER + "".join(f"{row}\n" for row in calls2),
                encoding="utf-8",
            )
            depth2_path = case / "depth-2.bed"
            depth2_path.write_text(
                "".join(f"{row}\n" for row in (depth2 or default_depth)),
                encoding="utf-8",
            )
            result["vcfs"].append(vcf2_path)  # type: ignore[union-attr]
            result["depths"].append(depth2_path)  # type: ignore[union-attr]
        return result

    def run_fixture(
        self,
        fixture: dict[str, object],
        label: str,
        *,
        policy: str = "ignore",
        genes: bool = False,
        amino_acids: bool = True,
        min_call_depth: int = 1,
        min_call_qual: int = 0,
        indel_range: int = 0,
        max_ref_mismatches: int = 10,
        gzip_outputs: bool = False,
        expect_success: bool = True,
    ) -> tuple[subprocess.CompletedProcess[str], dict[str, Path]]:
        case = fixture["dir"]
        assert isinstance(case, Path)
        vcfs = fixture["vcfs"]
        depths = fixture["depths"]
        assert isinstance(vcfs, list) and isinstance(depths, list)
        gzip_suffix = ".gz" if gzip_outputs else ""
        outputs = {"contig": case / f"consensus-{label}.fa{gzip_suffix}"}
        command = [
            str(self.binary),
            "-ref",
            str(fixture["ref"]),
            "-inVCF",
            ",".join(str(path) for path in vcfs),
            "-depthF",
            ",".join(str(path) for path in depths),
            "-gff",
            str(fixture["gff"]),
            "-oCtg",
            str(outputs["contig"]),
            "-minCallDepth",
            ",".join(str(min_call_depth) for _ in depths),
            "-minAltReads",
            "1",
            "-minAltFreq",
            "0",
            "-minCallQual",
            str(min_call_qual),
            "-depthFilterScale",
            "0",
            "-maxDepthFilterScale",
            "0",
            "-indelRange",
            str(indel_range),
            "-maxRefMismatches",
            str(max_ref_mismatches),
            "-vcfFilterPolicy",
            policy,
        ]
        if genes:
            outputs["gene_nt"] = case / f"genes-{label}.fna{gzip_suffix}"
            command.extend(
                [
                    "-oGeneNT",
                    str(outputs["gene_nt"]),
                ]
            )
            if amino_acids:
                outputs["gene_aa"] = case / f"genes-{label}.faa{gzip_suffix}"
                command.extend(["-oGeneAA", str(outputs["gene_aa"])])
        completed = subprocess.run(
            command,
            cwd=case,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=30,
            check=False,
        )
        diagnostic = (
            f"command: {shlex.join(command)}\n"
            f"exit: {completed.returncode}\nstdout:\n{completed.stdout}\n"
            f"stderr:\n{completed.stderr}"
        )
        if expect_success:
            self.assertEqual(completed.returncode, 0, diagnostic)
            for output in outputs.values():
                self.assertTrue(output.is_file(), diagnostic)
        else:
            self.assertNotEqual(completed.returncode, 0, diagnostic)
        return completed, outputs

    @unittest.skipUnless(sys.platform.startswith("linux"), "Linux auto-gzip build")
    def test_linux_build_reads_and_writes_gzip_automatically(self) -> None:
        fixture = self.make_fixture(
            "linux-gzip",
            [("ctg", "AAAA")],
            [vcf_call("ctg", 1, "A", "G")],
        )

        def compress(path: Path) -> Path:
            compressed = path.with_name(path.name + ".gz")
            with path.open("rb") as source, gzip.open(compressed, "wb") as target:
                target.write(source.read())
            return compressed

        ref = fixture["ref"]
        gff = fixture["gff"]
        vcfs = fixture["vcfs"]
        depths = fixture["depths"]
        assert isinstance(ref, Path) and isinstance(gff, Path)
        assert isinstance(vcfs, list) and isinstance(depths, list)
        fixture["ref"] = compress(ref)
        fixture["gff"] = compress(gff)
        fixture["vcfs"] = [compress(path) for path in vcfs]
        fixture["depths"] = [compress(path) for path in depths]

        _, outputs = self.run_fixture(fixture, "gzip", gzip_outputs=True)
        with gzip.open(outputs["contig"], "rt", encoding="utf-8") as handle:
            rendered = handle.read()
        self.assertIn("\nGAAA\n", rendered)

    def test_final_contig_is_flushed_at_eof(self) -> None:
        fixture = self.make_fixture(
            "final-contig",
            [("first", "AAAAA"), ("last", "CCCCC")],
            [
                vcf_call("first", 1, "A", "G"),
                vcf_call("last", 5, "C", "T"),
            ],
        )
        _, outputs = self.run_fixture(fixture, "run")
        consensus = read_fasta(outputs["contig"])
        self.assertEqual(consensus["first"][1], "GAAAA")
        self.assertEqual(consensus["last"][1], "CCCCT")

    def test_snp_insertions_and_deletions_use_vcf_anchor_coordinates(self) -> None:
        fixture = self.make_fixture(
            "variant-placement",
            [("ctg", "AACCGGTT")],
            [
                vcf_call("ctg", 2, "A", "T"),
                vcf_call("ctg", 4, "C", "CAA"),
                vcf_call("ctg", 6, "GT", "G"),
            ],
        )
        _, outputs = self.run_fixture(fixture, "run")
        # Replacements are anchored at POS and applied right-to-left:
        # AACCGGTT -> ATCCAAGGT.
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "ATCCAAGGT")

    def test_mnp_application_identity_rejection_and_gene_boundary_stability(self) -> None:
        fixture = self.make_fixture(
            "mnp-boundary",
            [("ctg", "CCATGAAA")],
            [vcf_call("ctg", 2, "CA", "GT")],
            gff_rows=[
                "ctg\ttest\tCDS\t3\t8\t.\t+\t0\tID=cds;transl_table=11",
            ],
        )
        _, outputs = self.run_fixture(
            fixture, "applied", genes=True, amino_acids=False
        )
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "CGTTGAAA")
        # The equal-length MNP overlaps the CDS start, but cannot move it.
        self.assertEqual(read_fasta(outputs["gene_nt"])["ctg_1"][1], "TTGAAA")

        identical = self.make_fixture(
            "mnp-identical-alt",
            [("ctg", "AACCGG")],
            [vcf_call("ctg", 3, "CC", "CC")],
        )
        failed, _ = self.run_fixture(identical, "rejected", expect_success=False)
        self.assertIn(
            "selected ALT allele is identical to REF",
            failed.stdout + failed.stderr,
        )

    def test_overlapping_replacements_at_different_anchors_mask_union(self) -> None:
        fixture = self.make_fixture(
            "overlapping-replacements",
            [("ctg", "AACCGGTT")],
            [
                vcf_call("ctg", 2, "ACC", "A"),
                vcf_call("ctg", 4, "CGG", "TTA"),
            ],
        )
        completed, outputs = self.run_fixture(fixture, "run")
        header, sequence = read_fasta(outputs["contig"])["ctg"]
        self.assertEqual(sequence, "ANNNNNTT")
        self.assertRegex(header, r"(?:^| )CONFL=1(?: |$)")
        self.assertIn("Conflicts resolved: 1", completed.stdout)
        self.assertIn("Warning: overlapping sequence replacements", completed.stderr)

    def test_delins_spanning_multipart_cds_emits_replacement_once(self) -> None:
        fixture = self.make_fixture(
            "multipart-delins",
            [("ctg", "AAAAAACCCCGGGGGG")],
            [
                # REF covers the end of CDS segment 1, the complete inter-segment
                # gap, and the start of segment 2. The inserted TT belongs in the
                # reconstructed CDS exactly once.
                vcf_call("ctg", 5, "AACCCCGG", "ATT"),
            ],
            gff_rows=[
                "ctg\ttest\tCDS\t1\t6\t.\t+\t0\tID=part1;Parent=tx;transl_table=11",
                "ctg\ttest\tCDS\t11\t16\t.\t+\t0\tID=part2;Parent=tx;transl_table=11",
            ],
        )
        _, outputs = self.run_fixture(
            fixture, "run", genes=True, amino_acids=False
        )
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "AAAAATTGGGG")
        genes = read_fasta(outputs["gene_nt"])
        self.assertEqual({key: value[1] for key, value in genes.items()}, {"ctg_1": "AAAAATTGGGG"})

    def test_mnp_and_interior_snp_mask_their_reference_union(self) -> None:
        fixture = self.make_fixture(
            "mnp-interior-snp",
            [("ctg", "AACCGG")],
            [
                vcf_call("ctg", 2, "ACC", "TTA"),
                vcf_call("ctg", 3, "C", "G"),
            ],
        )
        completed, outputs = self.run_fixture(fixture, "run")
        header, sequence = read_fasta(outputs["contig"])["ctg"]
        self.assertEqual(sequence, "ANNNGG")
        self.assertRegex(header, r"(?:^| )CONFL=1(?: |$)")
        self.assertRegex(header, r"(?:^| )REPL=0(?: |$)")
        self.assertIn("Found 1 SNPs and 1 INDELS", completed.stdout)

    def test_filter_column_policies(self) -> None:
        fixture = self.make_fixture(
            "filter-policy",
            [("ctg", "AAAAA")],
            [
                vcf_call("ctg", 1, "A", "G", filt="ExcessHet"),
                vcf_call("ctg", 2, "A", "T", filt="LowQual"),
                vcf_call("ctg", 3, "A", "C", filt="."),
            ],
        )
        expected = {
            "technical": "GNCAA",  # caller-model label ignored; technical failure rejected
            "all": "NNCAA",        # every named failure rejected
            "ignore": "GTCAA",     # FILTER column intentionally ignored
        }
        for policy, sequence in expected.items():
            with self.subTest(policy=policy):
                _, outputs = self.run_fixture(fixture, policy, policy=policy)
                self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], sequence)

    def test_technical_filter_uses_exact_model_only_ids(self) -> None:
        fixture = self.make_fixture(
            "technical-filter-exact-ids",
            [("ctg", "AAAA")],
            [
                vcf_call("ctg", 1, "A", "G", filt="SomaticLowQual"),
                vcf_call("ctg", 2, "A", "T", filt="DiploidDepth"),
                vcf_call("ctg", 3, "A", "C", filt="Somatic"),
                vcf_call("ctg", 4, "A", "G", filt="ExcessHet"),
            ],
        )
        completed, outputs = self.run_fixture(
            fixture, "run", policy="technical"
        )
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "NNCG")
        self.assertIn("VCF_FILTER: 2", completed.stdout)
        self.assertIn("VCF_FILTER_ignored: 2", completed.stdout)

    def test_missing_qual_obeys_minimum_quality_setting(self) -> None:
        fixture = self.make_fixture(
            "missing-qual",
            [("ctg", "AAA")],
            [vcf_call("ctg", 2, "A", "G", qual=".")],
        )

        filtered, filtered_outputs = self.run_fixture(
            fixture, "required", min_call_qual=20
        )
        self.assertEqual(read_fasta(filtered_outputs["contig"])["ctg"][1], "ANA")
        self.assertIn("Filtering stats: QUAL: 1", filtered.stdout)

        accepted, accepted_outputs = self.run_fixture(
            fixture, "optional", min_call_qual=0
        )
        self.assertEqual(read_fasta(accepted_outputs["contig"])["ctg"][1], "AGA")
        self.assertIn("Filtering stats: QUAL: 0", accepted.stdout)

    def test_multiallelic_ad_selects_highest_alt_and_uses_total_depth(self) -> None:
        fixture = self.make_fixture(
            "multiallelic",
            [("ctg", "AAAA")],
            [
                vcf_call("ctg", 1, "A", "C,G,T", ad=(1, 4, 10, 2)),
                # G is again the strongest ALT, but 8/20 is a minor allele.  A
                # denominator using only REF+selected ALT would incorrectly
                # make this a major call (8/13) and emit G rather than N.
                vcf_call("ctg", 2, "A", "C,G,T", ad=(5, 7, 8, 0)),
            ],
        )
        _, outputs = self.run_fixture(fixture, "run")
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "GNAA")

    def test_symbolic_multiallelic_placeholders_are_not_selected_as_bases(self) -> None:
        fixture = self.make_fixture(
            "symbolic-multiallelic",
            [("ctg", "AAA")],
            [
                # Symbolic gVCF/spanning-deletion alleles are not candidates
                # for a sequence consensus. Their evidence still contributes
                # to total depth, so the real G call remains a minor allele.
                vcf_call("ctg", 1, "A", "G,<NON_REF>", ad=(1, 10, 20)),
                vcf_call("ctg", 2, "A", "G,*", ad=(1, 10, 20)),
                vcf_call("ctg", 3, "A", "G"),
            ],
        )
        _, outputs = self.run_fixture(fixture, "run")
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "NNG")

    def test_crlf_vcf_records_are_accepted(self) -> None:
        fixture = self.make_fixture(
            "crlf-vcf",
            [("ctg", "AAAA")],
            [vcf_call("ctg", 2, "A", "G")],
        )
        vcf_path = fixture["vcfs"][0]  # type: ignore[index]
        assert isinstance(vcf_path, Path)
        unix_text = vcf_path.read_text(encoding="utf-8")
        vcf_path.write_bytes(unix_text.replace("\n", "\r\n").encode("utf-8"))
        _, outputs = self.run_fixture(fixture, "run")
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "AGAA")

    def test_filtered_major_deletion_masks_its_reference_span(self) -> None:
        fixture = self.make_fixture(
            "filtered-indel",
            [("ctg", "AACCGG")],
            [vcf_call("ctg", 3, "CCG", "C", filt="LowQual")],
        )
        _, outputs = self.run_fixture(fixture, "run", policy="technical")
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "AANNNG")

    def test_idv_or_imf_can_supply_indel_evidence_without_ad_or_dp4(self) -> None:
        fixture = self.make_fixture(
            "idv-imf-indels",
            [("ctg", "AAAA")],
            [
                # IDV/DP supplies both the alternate count and frequency.
                vcf_info_call("ctg", 1, "A", "AT", "DP=20;IDV=18"),
                # IMF/DP supplies both the frequency and alternate count.
                vcf_info_call("ctg", 4, "A", "AG", "DP=20;IMF=0.9"),
            ],
        )
        completed, outputs = self.run_fixture(fixture, "run")
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "ATAAAG")
        self.assertIn("Passing Filters: 0, 0; 2 entries", completed.stdout)

    def test_duplicate_position_prefers_passing_evidence(self) -> None:
        fixture = self.make_fixture(
            "duplicate-position",
            [("ctg", "AAAA")],
            [
                # The higher-frequency decomposition is technically filtered;
                # the lower-frequency but still-major passing call must survive.
                vcf_call("ctg", 2, "A", "G", ad=(1, 19), filt="LowQual"),
                vcf_call("ctg", 2, "A", "T", ad=(8, 12), filt="PASS"),
            ],
        )
        _, outputs = self.run_fixture(fixture, "run", policy="technical")
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "ATAA")

    def test_reference_mismatches_abort_only_after_limit(self) -> None:
        reference = [("ctg", "A" * 12)]
        mismatch_calls = [vcf_call("ctg", pos, "C", "G") for pos in range(1, 12)]

        tolerated = self.make_fixture("mismatch-tolerated", reference, mismatch_calls[:10])
        completed, outputs = self.run_fixture(
            tolerated, "ten", max_ref_mismatches=10
        )
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "A" * 12)
        self.assertIn("Reference/coordinate mismatches skipped: 10", completed.stdout)

        excessive = self.make_fixture("mismatch-excessive", reference, mismatch_calls)
        failed, _ = self.run_fixture(
            excessive,
            "eleven",
            max_ref_mismatches=10,
            expect_success=False,
        )
        self.assertIn(
            "Reference/VCF mismatch limit exceeded (11 > 10)",
            failed.stdout + failed.stderr,
        )

    def test_gff_ids_feature_selection_phase_reverse_and_multipart_cds(self) -> None:
        sequence = list("C" * 50)

        def place(one_based_start: int, value: str) -> None:
            start = one_based_start - 1
            sequence[start : start + len(value)] = value

        place(1, "AATGAAATAA")       # phase 1 -> ATG AAA TAA
        place(14, "TTAGGGCAT")        # reverse complement -> ATG CCC TAA
        place(26, "ATG")
        place(33, "AAATAA")           # multipart -> ATG AAA TAA

        rows = [
            # A parent gene must not duplicate its CDS child on a contig that
            # has coding features.
            "ctg\ttest\tgene\t1\t10\t.\t+\t.\tID=parent;gene_biotype=protein_coding",
            "ctg\ttest\tCDS\t1\t10\t.\t+\t1\tID=phase_cds;transl_table=11",
            "ctg\ttest\tmRNA\t1\t10\t.\t+\t.\tID=transcript",
            "ctg\ttest\texon\t1\t10\t.\t+\t.\tID=exon",
            "ctg\ttest\tCDS\t14\t22\t.\t-\t0\tID=reverse;transl_table=11",
            "ctg\ttest\tCDS\t26\t28\t.\t+\t0\tParent=multi;transl_table=11",
            "ctg\ttest\tCDS\t33\t38\t.\t+\t0\tParent=multi;transl_table=11",
            "ctg\ttest\trRNA\t40\t45\t.\t+\t.\tID=rrna",
            # On a contig with no explicit CDS, a protein-coding gene remains
            # a useful fallback, while a non-coding gene is excluded.
            "fallback\ttest\tgene\t1\t9\t.\t+\t.\tID=only_gene;gene_biotype=protein_coding;transl_table=11",
            "fallback\ttest\tgene\t1\t9\t.\t+\t.\tID=noncoding;gene_biotype=rRNA",
        ]
        fixture = self.make_fixture(
            "gff",
            [
                ("ctg descriptive FASTA header", "".join(sequence)),
                ("fallback another description", "ATGAAATAA"),
            ],
            [],
            gff_rows=rows,
        )
        _, outputs = self.run_fixture(fixture, "run", genes=True)
        contigs = read_fasta(outputs["contig"])
        nt = read_fasta(outputs["gene_nt"])
        aa = read_fasta(outputs["gene_aa"])

        self.assertTrue(contigs["ctg"][0].startswith("ctg descriptive FASTA header "))
        self.assertEqual(
            {key: value[1] for key, value in nt.items()},
            {
                "ctg_1": "AATGAAATAA",
                "ctg_2": "ATGCCCTAA",
                "ctg_3": "ATGAAATAA",
                "fallback_1": "ATGAAATAA",
            },
        )
        self.assertEqual(
            {key: value[1] for key, value in aa.items()},
            {"ctg_1": "MK*", "ctg_2": "MP*", "ctg_3": "MK*", "fallback_1": "MK*"},
        )

    def test_pseudo_true_cds_is_excluded(self) -> None:
        fixture = self.make_fixture(
            "pseudo-cds",
            [("ctg", "ATGAAATAAATGCCCTAA")],
            [],
            gff_rows=[
                "ctg\ttest\tCDS\t1\t9\t.\t+\t0\tID=pseudo_cds;pseudo=true;transl_table=11",
                "ctg\ttest\tCDS\t10\t18\t.\t+\t0\tID=valid_cds;transl_table=11",
            ],
        )
        completed, outputs = self.run_fixture(fixture, "run", genes=True)
        nucleotide_genes = read_fasta(outputs["gene_nt"])
        proteins = read_fasta(outputs["gene_aa"])
        self.assertEqual(
            {key: value[1] for key, value in nucleotide_genes.items()},
            {"ctg_1": "ATGCCCTAA"},
        )
        self.assertEqual(
            {key: value[1] for key, value in proteins.items()},
            {"ctg_1": "MP*"},
        )
        self.assertIn(
            "GFF features retained: 1; skipped non-gene/non-coding rows: 0; "
            "skipped gene parents where coding features were present: 0; "
            "skipped nonfunctional/pseudo features: 1",
            completed.stdout,
        )

    def test_pseudo_status_is_inherited_from_gene_through_mrna_to_cds(self) -> None:
        fixture = self.make_fixture(
            "inherited-pseudo-cds",
            [("ctg", "ATGAAATAAATGCCCTAA")],
            [],
            gff_rows=[
                "ctg\ttest\tgene\t1\t9\t.\t+\t.\tID=pseudo_gene;pseudo=true",
                "ctg\ttest\tmRNA\t1\t9\t.\t+\t.\tID=pseudo_tx;Parent=pseudo_gene",
                # Every translation-specific problem below must be ignored
                # because pseudo status is inherited before CDS validation.
                "ctg\ttest\tCDS\t1\t9\t.\t+\t.\tID=pseudo_child;Parent=pseudo_tx;"
                "transl_table=99;transl_except=(pos:1..3,aa:Sec)",
                "ctg\ttest\tCDS\t10\t18\t.\t+\t0\tID=valid;transl_table=11",
            ],
        )
        completed, outputs = self.run_fixture(fixture, "run", genes=True)
        self.assertEqual(
            {key: value[1] for key, value in read_fasta(outputs["gene_nt"]).items()},
            {"ctg_1": "ATGCCCTAA"},
        )
        self.assertEqual(
            {key: value[1] for key, value in read_fasta(outputs["gene_aa"]).items()},
            {"ctg_1": "MP*"},
        )
        self.assertIn("skipped nonfunctional/pseudo features: 1", completed.stdout)

    def test_overlapping_multipart_cds_segments_are_rejected(self) -> None:
        fixture = self.make_fixture(
            "overlapping-cds",
            [("ctg", "ATGAAATTTTAA")],
            [],
            gff_rows=[
                "ctg\ttest\tCDS\t1\t6\t.\t+\t0\tID=part1;Parent=tx;transl_table=11",
                "ctg\ttest\tCDS\t5\t12\t.\t+\t0\tID=part2;Parent=tx;transl_table=11",
            ],
        )
        failed, _ = self.run_fixture(
            fixture, "run", genes=True, expect_success=False
        )
        self.assertIn(
            "Overlapping or duplicate GFF segments",
            failed.stdout + failed.stderr,
        )

    def test_sequence_ontology_cds_requires_phase_but_orf_does_not(self) -> None:
        missing_phase = self.make_fixture(
            "so-cds-missing-phase",
            [("ctg", "ATGAAATAA")],
            [],
            gff_rows=[
                "ctg\ttest\tSO:0000316\t1\t9\t.\t+\t.\tID=cds;transl_table=11",
            ],
        )
        failed, _ = self.run_fixture(
            missing_phase, "rejected", genes=True, expect_success=False
        )
        self.assertIn(
            "CDS feature is missing its required GFF3 phase",
            failed.stdout + failed.stderr,
        )

        orf = self.make_fixture(
            "so-orf-no-phase",
            [("ctg", "ATGAAATAA")],
            [],
            gff_rows=[
                "ctg\ttest\tSO:0000236\t1\t9\t.\t+\t.\tID=orf;transl_table=11",
            ],
        )
        _, outputs = self.run_fixture(orf, "retained", genes=True)
        self.assertEqual(read_fasta(outputs["gene_nt"])["ctg_1"][1], "ATGAAATAA")
        self.assertEqual(read_fasta(outputs["gene_aa"])["ctg_1"][1], "MK*")

    def test_bacterial_translation_tables_4_11_and_25(self) -> None:
        fixture = self.make_fixture(
            "translation-tables",
            [
                ("table4", "ATATGATAA"),
                ("table11", "GTGTGATAA"),
                ("table25", "GTGTGATAA"),
            ],
            [],
            gff_rows=[
                "table4\ttest\tCDS\t1\t9\t.\t+\t0\tID=t4;transl_table=4",
                "table11\ttest\tCDS\t1\t9\t.\t+\t0\tID=t11;transl_table=11",
                "table25\ttest\tCDS\t1\t9\t.\t+\t0\tID=t25;transl_table=25",
            ],
        )
        _, outputs = self.run_fixture(fixture, "run", genes=True)
        proteins = read_fasta(outputs["gene_aa"])
        self.assertEqual(proteins["table4_1"][1], "MW*")
        self.assertEqual(proteins["table11_1"][1], "M*")
        self.assertEqual(proteins["table25_1"][1], "MG*")

    def test_unsupported_translation_table_is_only_fatal_for_aa_output(self) -> None:
        fixture = self.make_fixture(
            "unsupported-translation-table",
            [("ctg", "ATGAAATAA")],
            [],
            gff_rows=[
                "ctg\ttest\tCDS\t1\t9\t.\t+\t0\tID=cds;transl_table=99",
            ],
        )

        _, nt_outputs = self.run_fixture(
            fixture, "nt-only", genes=True, amino_acids=False
        )
        self.assertEqual(read_fasta(nt_outputs["gene_nt"])["ctg_1"][1], "ATGAAATAA")

        failed, _ = self.run_fixture(
            fixture, "with-aa", genes=True, expect_success=False
        )
        self.assertIn(
            "Unsupported NCBI translation table 99",
            failed.stdout + failed.stderr,
        )

    def test_malformed_numeric_gff_coordinates_are_rejected(self) -> None:
        fixture = self.make_fixture(
            "malformed-gff-coordinate",
            [("ctg", "ATGAAATAA")],
            [],
            gff_rows=[
                "ctg\ttest\tCDS\t1junk\t9\t.\t+\t0\tID=cds;transl_table=11",
            ],
        )
        failed, _ = self.run_fixture(fixture, "run", expect_success=False)
        self.assertIn("Invalid GFF3 coordinates", failed.stdout + failed.stderr)

    def test_reference_ns_and_depth_mask_have_distinct_accounting(self) -> None:
        fixture = self.make_fixture(
            "n-and-depth-mask",
            [("masked", "ANAAA")],
            [
                vcf_call("masked", 1, "A", "G"),
                vcf_call("masked", 4, "A", "T"),
            ],
            depth1=["masked\t0\t3\t20"],
        )
        completed, outputs = self.run_fixture(fixture, "run")
        header, sequence = read_fasta(outputs["contig"])["masked"]
        self.assertEqual(sequence, "GNANN")
        self.assertRegex(header, r"(?:^| )COV=2(?: |$)")
        self.assertRegex(header, r"(?:^| )REPL=1(?: |$)")
        self.assertIn("low Coverage: 1", completed.stdout)
        self.assertIn("Total bp that can be determined: 2", completed.stdout)

    def test_each_depth_source_reports_its_own_kept_positions(self) -> None:
        fixture = self.make_fixture(
            "per-source-depth-summary",
            [("ctg", "AAAAAA")],
            [],
            depth1=["ctg\t0\t6\t20"],
            calls2=[],
            depth2=["ctg\t1\t4\t20"],
        )
        completed, outputs = self.run_fixture(fixture, "run")
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "AAAAAA")
        self.assertIn(
            "Depth filter (1): Number of positions kept: 6 removed: 0",
            completed.stdout,
        )
        self.assertIn(
            "Depth filter (1): Number of positions kept: 3 removed: 3",
            completed.stdout,
        )

    def test_zero_minimum_depth_unmasks_sparse_and_omitted_intervals(self) -> None:
        fixture = self.make_fixture(
            "zero-minimum-sparse-depth",
            [("covered", "AAAA"), ("omitted", "CCCC")],
            [],
            depth1=["covered\t1\t3\t5"],
        )
        completed, outputs = self.run_fixture(
            fixture, "run", min_call_depth=0
        )
        consensus = read_fasta(outputs["contig"])
        self.assertEqual(consensus["covered"][1], "AAAA")
        self.assertEqual(consensus["omitted"][1], "CCCC")
        self.assertRegex(consensus["covered"][0], r"(?:^| )COV=4(?: |$)")
        self.assertRegex(consensus["omitted"][0], r"(?:^| )COV=4(?: |$)")
        self.assertIn(
            "Depth filter (0): Number of positions kept: 8 removed: 0",
            completed.stdout,
        )

    def test_two_base_gene_depth_accumulation_does_not_overflow(self) -> None:
        fixture = self.make_fixture(
            "high-gene-depth",
            [("ctg", "AA")],
            [],
            gff_rows=[
                "ctg\ttest\tgene\t1\t2\t.\t+\t.\tID=deep;gene_biotype=protein_coding",
            ],
            depth1=["ctg\t0\t2\t2147483647"],
        )
        _, outputs = self.run_fixture(
            fixture, "run", genes=True, amino_acids=False
        )
        header, sequence = read_fasta(outputs["gene_nt"])["ctg_1"]
        self.assertEqual(sequence, "AA")
        match = re.search(r"(?:^| )D=([0-9]+(?:\.[0-9]+)?)(?: |$)", header)
        self.assertIsNotNone(match, header)
        assert match is not None
        self.assertGreater(float(match.group(1)), 2_000_000_000.0)

    def test_alt_n_masks_site_and_is_accounted_as_uncertain(self) -> None:
        fixture = self.make_fixture(
            "ambiguous-alt",
            [("ctg", "AAAA")],
            [vcf_call("ctg", 2, "A", "N")],
        )
        completed, outputs = self.run_fixture(fixture, "run")
        header, sequence = read_fasta(outputs["contig"])["ctg"]
        self.assertEqual(sequence, "ANAA")
        self.assertRegex(header, r"(?:^| )COV=3(?: |$)")
        self.assertRegex(header, r"(?:^| )REPL=0(?: |$)")
        self.assertIn("POS= FR=", header)
        self.assertIn("Unsure: 1;0", completed.stdout)
        self.assertIn("ambiguous_ALT: 1", completed.stdout)

    def test_dual_vcf_agreement_conflict_and_second_only_calls(self) -> None:
        first = [
            vcf_call("shared", 1, "A", "G"),
            vcf_call("shared", 2, "A", "C"),
            vcf_call("shared", 4, "A", "C"),
            vcf_call("shared", 5, "A", "C"),
            vcf_call("shared", 6, "A", "G", filt="LowQual"),
        ]
        second = [
            vcf_call("shared", 1, "A", "G"),
            vcf_call("shared", 2, "A", "T"),
            vcf_call("shared", 3, "A", "G"),
            vcf_call("shared", 5, "A", "T", filt="LowQual"),
            vcf_call("shared", 6, "A", "G", filt="LowQual"),
            vcf_call("only2", 2, "A", "T"),
        ]
        fixture = self.make_fixture(
            "dual-vcf",
            [("shared", "AAAAAA"), ("only2", "AAA")],
            first,
            calls2=second,
        )
        completed, outputs = self.run_fixture(
            fixture, "run", policy="technical"
        )
        consensus = read_fasta(outputs["contig"])
        self.assertEqual(consensus["shared"][1], "GNGCCN")
        self.assertEqual(consensus["only2"][1], "ATA")
        self.assertRegex(consensus["shared"][0], r"(?:^| )CONFL=1(?: |$)")
        self.assertIn("Conflicts resolved: 1", completed.stdout)

    def test_mismatched_secondary_ref_does_not_conflict_with_valid_primary(self) -> None:
        fixture = self.make_fixture(
            "secondary-ref-mismatch",
            [("ctg", "AAAA")],
            [vcf_call("ctg", 2, "A", "G")],
            calls2=[vcf_call("ctg", 2, "C", "T")],
        )
        completed, outputs = self.run_fixture(fixture, "run")
        header, sequence = read_fasta(outputs["contig"])["ctg"]
        self.assertEqual(sequence, "AGAA")
        self.assertRegex(header, r"(?:^| )CONFL=0(?: |$)")
        self.assertIn("Reference/coordinate mismatches skipped: 1", completed.stdout)
        self.assertIn("Warning: skipping variant on ctg", completed.stderr)

    def test_mismatched_indel_does_not_proximity_filter_valid_snp(self) -> None:
        fixture = self.make_fixture(
            "mismatched-indel-proximity",
            [("ctg", "AAAAAAAAAA")],
            [
                vcf_call("ctg", 2, "C", "CTT"),
                vcf_call("ctg", 4, "A", "G"),
            ],
        )
        completed, outputs = self.run_fixture(
            fixture, "run", indel_range=5
        )
        self.assertEqual(read_fasta(outputs["contig"])["ctg"][1], "AAAGAAAAAA")
        self.assertIn("Reference/coordinate mismatches skipped: 1", completed.stdout)
        self.assertIn("indelProx: 0", completed.stdout)

    def test_dual_passing_snp_and_indel_at_same_position_conflict(self) -> None:
        fixture = self.make_fixture(
            "dual-snp-indel",
            [("ctg", "AAAA")],
            [vcf_call("ctg", 2, "A", "G")],
            calls2=[vcf_call("ctg", 2, "A", "AT")],
        )
        completed, outputs = self.run_fixture(fixture, "run", policy="technical")
        header, sequence = read_fasta(outputs["contig"])["ctg"]
        self.assertEqual(sequence, "ANAA")
        self.assertRegex(header, r"(?:^| )CONFL=1(?: |$)")
        self.assertIn("Conflicts resolved: 1", completed.stdout)

    def test_indels_at_cds_boundaries_preserve_gene_coordinates(self) -> None:
        fixture = self.make_fixture(
            "cds-boundary-indels",
            [("insert", "ATGAAAC"), ("delete", "CATGAAA")],
            [
                # Inserted TT occurs after the final CDS base and must not be
                # pulled into the gene merely because A is the VCF anchor.
                vcf_call("insert", 6, "A", "ATT"),
                # C is the left anchor outside the CDS; deleted AT are the
                # first two CDS bases, leaving GAAA as the mutated CDS.
                vcf_call("delete", 1, "CAT", "C"),
            ],
            gff_rows=[
                "insert\ttest\tCDS\t1\t6\t.\t+\t0\tID=insert_cds;transl_table=11",
                "delete\ttest\tCDS\t2\t7\t.\t+\t0\tID=delete_cds;transl_table=11",
            ],
        )
        _, outputs = self.run_fixture(fixture, "run", genes=True)
        contigs = read_fasta(outputs["contig"])
        genes = read_fasta(outputs["gene_nt"])
        self.assertEqual(contigs["insert"][1], "ATGAAATTC")
        self.assertEqual(contigs["delete"][1], "CGAAA")
        self.assertEqual(genes["insert_1"][1], "ATGAAA")
        self.assertEqual(genes["delete_1"][1], "GAAA")

    def test_snp_header_positions_follow_upstream_indel_shifts(self) -> None:
        fixture = self.make_fixture(
            "post-indel-snp-positions",
            [("insert", "AAAAAA"), ("delete", "AAAAAAAA")],
            [
                vcf_call("insert", 2, "A", "ATT"),
                vcf_call("insert", 5, "A", "G"),
                vcf_call("delete", 2, "AAA", "A"),
                vcf_call("delete", 7, "A", "G"),
            ],
        )
        _, outputs = self.run_fixture(fixture, "run")
        consensus = read_fasta(outputs["contig"])

        insert_header, insert_sequence = consensus["insert"]
        self.assertEqual(insert_sequence, "AATTAAGA")
        self.assertRegex(insert_header, r"(?:^| )REPL=1(?: |$)")
        self.assertRegex(insert_header, r"(?:^| )POS=6(?: |$)")

        delete_header, delete_sequence = consensus["delete"]
        self.assertEqual(delete_sequence, "AAAAGA")
        self.assertRegex(delete_header, r"(?:^| )REPL=1(?: |$)")
        self.assertRegex(delete_header, r"(?:^| )POS=4(?: |$)")

    def test_output_path_collisions_are_rejected_before_truncation(self) -> None:
        fixture = self.make_fixture(
            "output-path-collisions",
            [("ctg", "AAAA")],
            [],
        )
        case = fixture["dir"]
        vcfs = fixture["vcfs"]
        depths = fixture["depths"]
        assert isinstance(case, Path)
        assert isinstance(vcfs, list) and isinstance(depths, list)
        common = [
            str(self.binary),
            "-ref",
            str(fixture["ref"]),
            "-inVCF",
            ",".join(str(path) for path in vcfs),
            "-depthF",
            ",".join(str(path) for path in depths),
            "-gff",
            str(fixture["gff"]),
        ]

        shared_output = case / "shared-output.fa"
        sentinel = "do-not-truncate\n"
        shared_output.write_text(sentinel, encoding="utf-8")
        duplicate = subprocess.run(
            [
                *common,
                "-oCtg",
                str(shared_output),
                "-oGeneNT",
                str(shared_output),
            ],
            cwd=case,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=30,
            check=False,
        )
        self.assertNotEqual(duplicate.returncode, 0)
        self.assertIn("must use different output paths", duplicate.stderr)
        self.assertEqual(shared_output.read_text(encoding="utf-8"), sentinel)

        ref_path = fixture["ref"]
        assert isinstance(ref_path, Path)
        original_reference = ref_path.read_text(encoding="utf-8")
        overwrite = subprocess.run(
            [*common, "-oCtg", str(ref_path)],
            cwd=case,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=30,
            check=False,
        )
        self.assertNotEqual(overwrite.returncode, 0)
        self.assertIn("must not overwrite an input file", overwrite.stderr)
        self.assertEqual(ref_path.read_text(encoding="utf-8"), original_reference)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--binary",
        type=Path,
        help="use an existing Linux/WSL vcf2fna executable instead of compiling",
    )
    parser.add_argument(
        "--cxx",
        default=os.environ.get("CXX", "g++"),
        help="C++ compiler command used for the test build (default: CXX or g++)",
    )
    parser.add_argument("-q", "--quiet", action="store_true", help="compact test output")
    args = parser.parse_args()

    VCF2fnaRegressionTests.binary_override = args.binary
    VCF2fnaRegressionTests.cxx = args.cxx
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(VCF2fnaRegressionTests)
    result = unittest.TextTestRunner(verbosity=1 if args.quiet else 2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
