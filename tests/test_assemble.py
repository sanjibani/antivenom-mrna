"""Tests for the full-construct assembler.

The headline test in this file rebuilds a licensed mRNA vaccine from its
published parts and compares the result with the primary regulatory
document. If the assembler cannot reproduce a product that a regulator has
approved, it has no business assembling anything else.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mrna.assemble import (
    A110_POLYA,
    CAP1,
    Construct,
    ConstructPart,
    assemble_construct,
    full_construct,
)
from mrna.fasta import read_fasta
from mrna.motifs import DEFAULT_MOTIFS
from mrna.sequence import translate
from mrna.utr import AES_MTRNR1_3UTR, ALPHA_GLOBIN_5UTR

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"

#: The WHO document prints two consecutive stop codons at 3874 to 3879. The
#: published coding-sequence transcription ends at the first of them, because
#: that is where translation stops, so the second one has to be added back
#: before the construct can be reproduced.
SECOND_STOP_CODON = "TGA"

#: The element table published in WHO INN document 11889.
PUBLISHED_LAYOUT: tuple[tuple[str, int, int], ...] = (
    ("cap", 1, 2),
    ("5'-UTR", 3, 54),
    ("sig", 55, 102),
    ("S protein_mut", 103, 3879),
    ("3'-UTR", 3880, 4174),
    ("poly(A)", 4175, 4284),
)


def who_mrna() -> str:
    """Return the full Comirnaty mRNA parsed from the WHO document."""
    return read_fasta(REFERENCE / "bnt162b2_mrna_who.fasta")[0][1]


def published_cds() -> str:
    """Return the published Comirnaty coding sequence."""
    return read_fasta(REFERENCE / "bnt162b2_cds.fasta")[0][1]


def comirnaty() -> Construct:
    """Rebuild Comirnaty from the published parts."""
    return full_construct(
        "BNT162b2 (Comirnaty)",
        utr5=ALPHA_GLOBIN_5UTR,
        coding=published_cds() + SECOND_STOP_CODON,
        utr3=AES_MTRNR1_3UTR,
        poly_a=A110_POLYA,
        signal_peptide_nt=48,
        coding_name="S protein_mut",
        coding_description=(
            "Codon-optimized sequence encoding full-length SARS-CoV-2 spike "
            "(S) glycoprotein containing mutations K986P and V987P to ensure "
            "the S glycoprotein remains in an antigenically optimal "
            "pre-fusion conformation; stop codons: 3874-3879 (underlined)"
        ),
    )


def test_the_cap_is_two_printed_nucleotides() -> None:
    """The cap footprint is the m7G and the first backbone nucleotide."""
    assert CAP1.prefix == "GA"
    assert CAP1.length == 2


def test_the_assembler_rebuilds_comirnaty_base_for_base() -> None:
    """The assembled construct must equal the published sequence exactly."""
    rebuilt = comirnaty()
    document = who_mrna()
    assert len(rebuilt.sequence) == 4284
    assert rebuilt.sequence == document


def test_the_element_table_matches_the_published_table() -> None:
    """Every published element name and coordinate must be reproduced."""
    rebuilt = comirnaty()
    observed = tuple(
        (element.name, element.start, element.end) for element in rebuilt.elements
    )
    assert observed == PUBLISHED_LAYOUT


def test_the_element_table_uses_the_regulatory_format() -> None:
    """The table must have the three columns a regulatory document uses."""
    table = comirnaty().element_table()
    assert table.splitlines()[0] == "| Element | Description | Position |"
    assert "| cap |" in table
    assert "| 1-2 |" in table


def test_the_rebuilt_construct_translates_to_the_published_protein() -> None:
    """The reconstructed coding region must encode the same protein."""
    rebuilt = comirnaty()
    orf = rebuilt.open_reading_frame()
    assert orf.start == 54
    assert orf.end == 3876
    assert translate(orf.cds).rstrip("*") == translate(published_cds()).rstrip("*")


def test_the_two_stop_codons_sit_where_the_document_says() -> None:
    """The construct must print both stop codons of the published design."""
    assert comirnaty().sequence[3873:3879] == "TGATGA"


def test_the_poly_a_tail_is_the_published_one() -> None:
    """30 adenosines, the 10-nt linker, then 70 adenosines."""
    tail = comirnaty().part("poly(A)").sequence
    assert len(tail) == 110
    assert tail[:30] == "A" * 30
    assert tail[30:40] == "GCATATGACT"
    assert tail[40:] == "A" * 70


def test_the_coding_region_carries_no_forbidden_motif() -> None:
    """The coding element must be free of the polyadenylation signal."""
    rebuilt = comirnaty()
    hits = rebuilt.find_motifs(DEFAULT_MOTIFS)
    assert [hit for hit in hits if hit.element == "S protein_mut"] == []


def test_the_untranslated_regions_do_carry_the_polyadenylation_signal() -> None:
    """The licensed product carries AAUAAA twice, both outside the coding region.

    This is the reason motif scanning has to be scoped and reported per
    element. A polyadenylation signal inside a coding region can truncate
    the transcript; the same signal in a 5' UTR has no cleavage site
    upstream of it and is harmless.
    """
    rebuilt = comirnaty()
    hits = rebuilt.find_motifs(DEFAULT_MOTIFS)
    assert {hit.element for hit in hits} == {"5'-UTR", "3'-UTR"}
    assert len(hits) == 2


def test_the_five_prime_utr_has_no_upstream_start_codon() -> None:
    """A 5' UTR with an ATG would let the ribosome start in the wrong place."""
    assert comirnaty().upstream_starts() == ()


def test_an_upstream_start_codon_is_reported_with_its_context() -> None:
    """A synthetic upstream ATG must be found and its context judged."""
    construct = assemble_construct(
        "with an upstream start",
        [
            ConstructPart("cap", "GA"),
            ConstructPart("5'-UTR", "GCCACCATGGGG"),  # strong context at -3 A +4 G
            ConstructPart("antigen", "ATG" + "GCT" * 10),
        ],
    )
    starts = construct.upstream_starts()
    assert len(starts) == 1
    assert starts[0].index == 2 + 6
    assert starts[0].matches_consensus
    assert starts[0].is_strong


def test_a_weak_upstream_start_codon_is_found_but_not_flagged_strong() -> None:
    """A start codon with a poor context is reported but not called strong."""
    construct = assemble_construct(
        "with a weak upstream start",
        [
            ConstructPart("cap", "GA"),
            ConstructPart("5'-UTR", "TTTATGAAA"),
            ConstructPart("antigen", "ATG" + "GCT" * 10),
        ],
    )
    starts = construct.upstream_starts()
    assert len(starts) == 1
    assert not starts[0].matches_consensus


def test_the_assembly_is_deterministic() -> None:
    """The same inputs must produce the same construct every time."""
    assert comirnaty().sequence == comirnaty().sequence
    assert comirnaty().element_table() == comirnaty().element_table()


def test_the_signal_peptide_is_its_own_element() -> None:
    """A secreted antigen needs a signal peptide row of its own."""
    signal = comirnaty().part("sig")
    assert signal.sequence == who_mrna()[54:102]
    assert translate(signal.sequence) == "MFVFLVLLPLVSSQCV"


def test_a_construct_without_a_signal_peptide_has_no_sig_row() -> None:
    """A cytosolic antigen gets no signal peptide element."""
    construct = full_construct(
        "no signal peptide",
        utr5=ALPHA_GLOBIN_5UTR,
        coding="ATG" + "GCT" * 10 + "TGA",
        utr3=AES_MTRNR1_3UTR,
        poly_a=A110_POLYA,
    )
    assert "sig" not in {part.name for part in construct.parts}


def test_an_empty_part_is_rejected() -> None:
    """An empty element is a broken design, not an empty row."""
    with pytest.raises(ValueError, match="empty"):
        assemble_construct("broken", [ConstructPart("cap", "GA"), ConstructPart("x", "")])


def test_two_parts_with_the_same_name_are_rejected() -> None:
    """Duplicate element names would make the table ambiguous."""
    with pytest.raises(ValueError, match="two parts named"):
        assemble_construct("broken", [ConstructPart("cap", "GA"), ConstructPart("cap", "AT")])


def test_a_part_with_a_character_outside_the_dna_alphabet_is_rejected() -> None:
    """A part that is not DNA cannot be assembled."""
    with pytest.raises(ValueError, match="outside"):
        assemble_construct("broken", [ConstructPart("cap", "GAN")])
