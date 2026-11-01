"""Tests for the UTR library.

The library's contract is that every entry carries the source that justifies
it. These tests protect that contract as much as they protect the sequences.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mrna.fasta import read_fasta
from mrna.utr import AES_ELEMENT, AES_MTRNR1_3UTR, ALPHA_GLOBIN_5UTR, MTRNR1_ELEMENT, Utr

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"

ALL_ENTRIES: tuple[Utr, ...] = (
    ALPHA_GLOBIN_5UTR,
    AES_ELEMENT,
    MTRNR1_ELEMENT,
    AES_MTRNR1_3UTR,
)


def who_mrna() -> str:
    """Return the full Comirnaty mRNA parsed from the WHO document."""
    return read_fasta(REFERENCE / "bnt162b2_mrna_who.fasta")[0][1]


def test_the_five_prime_utr_is_the_published_fifty_two_nucleotides() -> None:
    """The 5' UTR must equal the document's 5'-UTR element at positions 3-54."""
    assert ALPHA_GLOBIN_5UTR.length == 52
    assert ALPHA_GLOBIN_5UTR.sequence == who_mrna()[2:54]


def test_the_five_prime_utr_ends_with_the_kozak_sequence() -> None:
    """GCCACC immediately before the start codon is the optimised Kozak."""
    assert ALPHA_GLOBIN_5UTR.sequence.endswith("GCCACC")


def test_the_three_prime_utr_is_the_published_two_hundred_ninety_five() -> None:
    """The 3' UTR must equal the document's 3'-UTR element at 3880-4174."""
    assert AES_MTRNR1_3UTR.length == 295
    assert AES_MTRNR1_3UTR.sequence == who_mrna()[3879:4174]


def test_the_three_prime_utr_decomposes_into_its_two_named_elements() -> None:
    """AES and the 12S rRNA element must account for the whole published UTR."""
    prefix = AES_MTRNR1_3UTR.sequence[:6]
    suffix = AES_MTRNR1_3UTR.sequence[-11:]
    assert AES_MTRNR1_3UTR.sequence == (
        prefix + AES_ELEMENT.sequence + MTRNR1_ELEMENT.sequence + suffix
    )
    assert prefix == "CTCGAG"


def test_the_two_elements_are_in_the_order_the_source_names_them() -> None:
    """AES comes first, then the mitochondrial 12S rRNA element."""
    assert AES_MTRNR1_3UTR.sequence.index(AES_ELEMENT.sequence) < AES_MTRNR1_3UTR.sequence.index(
        MTRNR1_ELEMENT.sequence,
    )


@pytest.mark.parametrize("entry", ALL_ENTRIES, ids=lambda entry: entry.name)
def test_every_entry_carries_a_citation(entry: Utr) -> None:
    """Nothing enters the library without a source."""
    assert entry.citation.strip()
    assert entry.description.strip()


@pytest.mark.parametrize("entry", ALL_ENTRIES, ids=lambda entry: entry.name)
def test_every_entry_is_plain_dna(entry: Utr) -> None:
    """A UTR is DNA in this project, written with T for U."""
    assert set(entry.sequence) <= set("ACGT")


@pytest.mark.parametrize("entry", ALL_ENTRIES, ids=lambda entry: entry.name)
def test_every_entry_is_not_extremely_at_rich_or_extremely_gc_rich(entry: Utr) -> None:
    """A region that is all one base or all one other base is a typo."""
    fraction = entry.gc_content()
    assert 0.2 < fraction < 0.8


def test_the_contains_helper_finds_a_literal_sequence() -> None:
    """The helper is used to check element boundaries."""
    assert ALPHA_GLOBIN_5UTR.contains("GCCACC")
    assert not ALPHA_GLOBIN_5UTR.contains("TTTTTT")
