"""Tests for the anatomy layer.

These tests are where the engine's claims are checked against the published
record. Every assertion here is either a published fact about a licensed
vaccine or a baseline that this repository now promises to reproduce.

Sources for the published facts:

* BNT162b2 coding sequence: BioNTech/WHO Comirnaty vaccine description
  document, as transcribed in berthubert/bnt162b2 and archived in
  ``data/reference/raw/``.
* BNT162b2 antigen mutations K986P and V987P: TGA, Nonclinical Evaluation
  Report for BNT162b2; MHRA, COVID-19 mRNA Vaccine BNT162b2 Public
  Assessment Report.
* SARS-CoV-2 spike reference protein: UniProt P0DTC2.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mrna.anatomy import (
    align_proteins,
    analyze_cds,
    codon_adaptation_index,
    codon_counts,
    find_orf,
    kozak_context,
    protein_variants,
    relative_adaptiveness,
    split_construct,
)
from mrna.fasta import read_fasta
from mrna.sequence import AA_TO_CODONS, CODON_TABLE, translate

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"
RAW = REFERENCE / "raw"

BNT162B2_CDS_LENGTH = 3822
BNT162B2_PROTEIN_LENGTH = 1273


def load_usage_table() -> dict[str, float]:
    """Load the committed Kazusa human codon usage table."""
    table: dict[str, float] = {}
    for line in (REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) >= 3 and fields[1] != "*":
            table[fields[0]] = float(fields[2])
    return table


def bnt162b2_cds() -> str:
    """Load the BNT162b2 coding sequence from the WHO-document transcription."""
    return read_fasta(REFERENCE / "bnt162b2_cds.fasta")[0][1]


def reference_protein() -> str:
    """Load the UniProt SARS-CoV-2 spike reference protein."""
    return read_fasta(REFERENCE / "sars_cov2_spike_P0DTC2.fasta")[0][1].strip()


def assemblies() -> tuple[str, str]:
    """Load the two read-derived vaccine mRNA assemblies."""
    records = read_fasta(RAW / "naalytics_figure1figure2_032321.fasta")
    return records[0][1], records[1][1]


def test_find_orf_locates_the_published_start_of_the_bnt162b2_orf() -> None:
    """The published ORF for BNT162b2 starts at nucleotide 55.

    The ORF coordinates are published in the Comirnaty assessment literature
    (nt 55 to 3879). Our ORF boundary includes the stop codon, so it closes
    three bases earlier, at nt 3876.
    """
    orf = find_orf(assemblies()[0])
    assert orf.start + 1 == 55
    assert orf.end == 3876
    assert len(orf.cds) == BNT162B2_CDS_LENGTH


def test_signal_peptide_is_48_nt_as_published() -> None:
    """The published signal peptide spans nt 55 to 102, that is 48 nt."""
    cds = bnt162b2_cds()
    assert len(cds[:48]) == 48
    assert translate(cds[:48]) == "MFVFLVLLPLVSSQCV"


def test_bnt162b2_cds_length_matches_the_published_record() -> None:
    """3,822 nt is 1,273 codons plus the stop codon."""
    cds = bnt162b2_cds()
    assert len(cds) == BNT162B2_CDS_LENGTH
    protein = translate(cds)
    assert protein.endswith("*")
    assert len(protein) - 1 == BNT162B2_PROTEIN_LENGTH


def test_bnt162b2_cds_ends_on_a_stop_codon() -> None:
    """A coding sequence for a therapeutic protein must close on a stop."""
    assert bnt162b2_cds()[-3:] == "TGA"


def test_who_viral_spike_translates_exactly_to_uniprot() -> None:
    """The document's viral S CDS and UniProt P0DTC2 must be the same protein.

    This is the check that validates the reference: if the viral sequence in
    the primary document did not reproduce the authoritative UniProt entry,
    every downstream comparison would be worthless.
    """
    viral_cds = read_fasta(RAW / "berthubert_ncov-s.fasta")[0][1]
    assert translate(viral_cds).rstrip("*") == reference_protein()


def test_vaccine_cds_equals_the_read_derived_assembly() -> None:
    """The primary document and the sequencing assembly must agree exactly.

    One is a human transcription of an official document, the other is a
    read-derived assembly. If they agree base for base across all 3,822
    coding nucleotides, neither is an artefact.
    """
    orf = find_orf(assemblies()[0])
    assert orf.cds == bnt162b2_cds()


def test_the_only_protein_changes_are_the_two_prolines() -> None:
    """BNT162b2 encodes K986P and V987P against the wild-type spike.

    Published in the TGA Nonclinical Evaluation Report for BNT162b2 and in
    the MHRA Public Assessment Report. This construct is the "2P" antigen:
    two prolines at the junction between the HR1 helix and the central
    helix, which lock the spike in its prefusion shape.
    """
    variants = protein_variants(reference_protein(), translate(bnt162b2_cds()).rstrip("*"))
    assert [(v.position, v.reference_aa, v.variant_aa) for v in variants] == [
        (986, "K", "P"),
        (987, "V", "P"),
    ]
    assert [v.describe() for v in variants] == ["K986P", "V987P"]


def test_both_vaccine_antigens_encode_the_same_protein() -> None:
    """The two licensed vaccines carry the same antigen, codon-optimised differently."""
    bnt_orf = find_orf(assemblies()[0])
    other_orf = find_orf(assemblies()[1])
    assert translate(bnt_orf.cds) == translate(other_orf.cds)


def test_the_two_codon_optimisations_differ_only_synonymously() -> None:
    """Different codon choices, identical protein: 363 nt across 294 codons."""
    bnt_orf = find_orf(assemblies()[0])
    other_orf = find_orf(assemblies()[1])
    differing_nt = sum(
        1 for left, right in zip(bnt_orf.cds, other_orf.cds, strict=True) if left != right
    )
    codons_a = [bnt_orf.cds[i : i + 3] for i in range(0, len(bnt_orf.cds), 3)]
    codons_b = [other_orf.cds[i : i + 3] for i in range(0, len(other_orf.cds), 3)]
    differing_codons = [
        index
        for index, (left, right) in enumerate(zip(codons_a, codons_b, strict=True))
        if left != right
    ]
    assert differing_nt == 363
    assert len(differing_codons) == 294
    assert all(
        CODON_TABLE[left] == CODON_TABLE[right]
        for left, right in zip(codons_a, codons_b, strict=True)
    )


def test_gc_and_cai_baselines_are_reproducible() -> None:
    """The measured numbers that this repository publishes must not drift."""
    usage = load_usage_table()
    report = analyze_cds("BNT162b2", "test", bnt162b2_cds(), usage)
    assert report.gc_total == pytest.approx(0.5699, abs=5e-5)
    assert report.cai == pytest.approx(0.9515, abs=5e-5)

    other = analyze_cds("mRNA-1273", "test", find_orf(assemblies()[1]).cds, usage)
    assert other.gc_total == pytest.approx(0.6230, abs=5e-5)
    assert other.cai == pytest.approx(0.9786, abs=5e-5)


def test_kozak_context_is_read_from_the_real_junction() -> None:
    """Both vaccines use the textbook GCCACC Kozak context.

    The context is GCCACC followed by ATG, so position -3 is a purine, as the
    consensus wants. Position +4 is T rather than G, which is a deviation
    from the consensus but evidently not a problem for two licensed products.
    """
    utr5, cds, _, _ = split_construct(assemblies()[0])
    assert utr5.endswith("GCCACC")
    context = kozak_context(f"{utr5[-3:]}{cds[:6]}", 3)
    assert context.context == "ACCATGTTC"
    assert context.minus3 == "A"
    assert context.plus4 == "T"
    assert context.matches_consensus is False
    assert context.deviations == ("+4 is not G",)


def test_split_construct_finds_the_four_parts() -> None:
    """A full mRNA must split into UTR, CDS, UTR and tail."""
    utr5, cds, utr3, poly_a = split_construct(assemblies()[0])
    assert len(utr5) == 54
    assert len(cds) == BNT162B2_CDS_LENGTH
    assert len(utr3) == 299
    assert poly_a == ""


def test_relative_adaptiveness_normalises_the_best_codon_to_one() -> None:
    """The best codon of each amino acid must have weight 1.0."""
    weights = relative_adaptiveness(load_usage_table())
    for amino_acid_codons in {
        "F": ["TTT", "TTC"],
        "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    }.values():
        assert max(weights[codon] for codon in amino_acid_codons) == pytest.approx(1.0)


def test_cai_is_one_for_a_perfectly_optimised_coding_sequence() -> None:
    """A sequence built only from preferred codons must score exactly 1.0."""
    usage = load_usage_table()
    weights = relative_adaptiveness(usage)
    best = {
        amino_acid: max(codons, key=lambda codon: (weights[codon], codon))
        for amino_acid, codons in AA_TO_CODONS.items()
    }
    cds = "".join(best[amino_acid] for amino_acid in "MKV")
    assert codon_adaptation_index(cds, usage) == pytest.approx(1.0)


def test_codon_counts_sum_to_the_coding_sequence() -> None:
    """Codon counts must account for every codon in the sequence."""
    cds = bnt162b2_cds()
    assert sum(codon_counts(cds).values()) == len(cds) // 3


def test_align_proteins_handles_an_n_terminal_extension() -> None:
    """A variant with extra N-terminal residues must still align."""
    reference = "MKVLA"
    variant = "GGMKVLA"
    offset, variants = align_proteins(reference, variant)
    assert offset == 2
    assert variants == ()


def test_align_proteins_reports_the_difference_positions() -> None:
    """Aligned proteins must report their differences in reference coordinates."""
    offset, variants = align_proteins("MKVLA", "MKVPA")
    assert offset == 0
    assert [v.describe() for v in variants] == ["L4P"]


def test_protein_variants_refuses_mismatched_lengths() -> None:
    """A strict diff must refuse sequences it cannot align."""
    with pytest.raises(ValueError, match="equal lengths"):
        protein_variants("MKV", "MKVLA")
