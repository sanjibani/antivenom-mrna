"""Tests for the primitive sequence layer.

The hard rule for this project is that no number enters the engine without a
source, and no number leaves it without a test. These tests cover the
primitives with hand-computed answers where a hand-computed answer exists,
and with property-based tests where the property is the real guarantee.
"""

from __future__ import annotations

import pytest
from hypothesis import given
from hypothesis import strategies as st

from mrna.sequence import (
    AA_TO_CODONS,
    CODON_TABLE,
    STOP_CODONS,
    back_translate,
    gc_content,
    translate,
)

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def full_usage_table() -> dict[str, float]:
    """Build a complete codon usage table for tests.

    The first codon of every amino acid gets usage 1.0 and the rest 0.5, so
    the table has a clear winner per amino acid and exercises the ranking
    logic without depending on a downloaded reference table.
    """
    table: dict[str, float] = {}
    for codons in AA_TO_CODONS.values():
        for index, codon in enumerate(codons):
            table[codon] = 1.0 if index == 0 else 0.5
    return table


def test_codon_table_covers_every_triplet() -> None:
    """All 64 codons must be present and must be three bases long."""
    assert len(CODON_TABLE) == 64
    assert all(len(codon) == 3 for codon in CODON_TABLE)
    assert set(CODON_TABLE) == {
        first + second + third
        for first in "ACGT"
        for second in "ACGT"
        for third in "ACGT"
    }


def test_codon_table_encodes_20_amino_acids_and_one_stop() -> None:
    """The standard code has 20 sense codons and 3 stop codons."""
    residues = set(CODON_TABLE.values())
    assert residues - {"*"} == set(AMINO_ACIDS)
    assert len(STOP_CODONS) == 3
    assert all(CODON_TABLE[codon] == "*" for codon in STOP_CODONS)


def test_aa_to_codons_groups_are_sorted_and_complete() -> None:
    """The synonym groups must cover every sense codon exactly once."""
    grouped = [codon for codons in AA_TO_CODONS.values() for codon in codons]
    assert len(grouped) == 61
    assert len(set(grouped)) == 61
    assert set(grouped) == {codon for codon, aa in CODON_TABLE.items() if aa != "*"}
    for codons in AA_TO_CODONS.values():
        assert list(codons) == sorted(codons)


def test_translate_known_segment() -> None:
    """A hand-read segment must translate to the hand-written protein."""
    assert translate("ATGGGCCATTGTGATGGCTGG") == "MGHCDGW"


def test_translate_ignores_case_and_whitespace() -> None:
    """Formatting must not change the biology."""
    assert translate("atg ggc cat tgt gat ggc tgg") == "MGHCDGW"


def test_translate_requires_a_multiple_of_three() -> None:
    """A partial codon cannot be translated, so it must be refused."""
    with pytest.raises(ValueError, match="multiple of 3"):
        translate("ATGGGCC")


def test_translate_rejects_non_acgt_characters() -> None:
    """Ambiguity codes must fail loudly by default."""
    with pytest.raises(ValueError, match="non-ACGT"):
        translate("ATGNNNGGC")
    assert translate("ATGNNNGGC", allow_ambiguous=True) == "MXG"


def test_translate_stop_at_stop_truncates_but_keeps_the_marker() -> None:
    """Truncation must stop at the first stop codon and keep the marker."""
    assert translate("ATGTAAATGAAA", stop_at_stop=True) == "M*"
    assert translate("ATGTAAATGAAA") == "M*M K".replace(" ", "")


def test_translate_renders_every_stop_codon() -> None:
    """The three stop codons must all render as the stop marker."""
    assert translate("TAATAGTGA") == "***"


def test_back_translate_round_trip_on_a_known_protein() -> None:
    """Back-translating then translating must return the original protein."""
    protein = "MGHCDGWYTSKNRQL"
    cds = back_translate(protein, full_usage_table())
    assert translate(cds) == protein


def test_back_translate_prefers_the_most_used_codon() -> None:
    """The highest-usage codon of each amino acid must be chosen."""
    table = full_usage_table()
    table.update(
        {"CTG": 0.9, "CTA": 0.05, "CTC": 0.05, "CTT": 0.05, "TTA": 0.05, "TTG": 0.05},
    )
    assert back_translate("MWL", table) == "ATGTGGCTG"


def test_back_translate_breaks_ties_deterministically() -> None:
    """Equal usage must resolve the same way on every run."""
    table = full_usage_table()
    table.update({"CTG": 0.5, "CTC": 0.5, "CTT": 0.5, "TTA": 0.5, "TTG": 0.5, "CTA": 0.5})
    first = back_translate("ML", table)
    assert first == back_translate("ML", table)
    assert first == "ATGCTA"


def test_back_translate_rejects_stop_codons_in_the_protein() -> None:
    """Stop codons belong to the assembler, not to the back-translator."""
    with pytest.raises(ValueError, match="stop codons"):
        back_translate("MG*W", full_usage_table())


def test_back_translate_reports_missing_usage_entries() -> None:
    """A table with a hole must produce a message naming the hole."""
    table = full_usage_table()
    del table["CTG"]
    with pytest.raises(ValueError, match="CTG"):
        back_translate("ML", table)


def test_gc_content_hand_calculated() -> None:
    """GC fractions must match values computed by hand."""
    assert gc_content("ATGC") == 0.5
    assert gc_content("GGGG") == 1.0
    assert gc_content("ATAT") == 0.0
    assert gc_content("ATGCGC") == pytest.approx(4 / 6)
    assert gc_content("AATTGGCC") == 0.5


def test_gc_content_window_values() -> None:
    """The sliding window must return the correct fraction for each window."""
    assert gc_content("ATGCGC", window=3) == pytest.approx([1 / 3, 2 / 3, 1.0, 1.0])


def test_gc_content_window_step_controls_stride() -> None:
    """A larger stride must return fewer windows."""
    every_base = gc_content("ATGCGC", window=3)
    every_other = gc_content("ATGCGC", window=3, step=2)
    assert len(every_base) == 4
    assert len(every_other) == 2
    assert every_other == every_base[::2]


def test_gc_content_rejects_an_oversized_window() -> None:
    """A window longer than the sequence is a bug in the caller."""
    with pytest.raises(ValueError, match="longer than the sequence"):
        gc_content("ATGC", window=10)


def test_gc_content_rejects_a_non_positive_window() -> None:
    """A zero or negative window is meaningless."""
    with pytest.raises(ValueError, match="must be positive"):
        gc_content("ATGC", window=0)


def test_gc_content_counts_ambiguity_in_the_denominator() -> None:
    """An N is neither G nor C but still occupies a position."""
    assert gc_content("GCNN") == 0.5


@given(st.text(alphabet="ACGT", min_size=1, max_size=500))
def test_gc_content_is_always_a_fraction(seq: str) -> None:
    """GC content can never leave the unit interval."""
    value = gc_content(seq)
    assert 0.0 <= value <= 1.0


@given(st.text(alphabet="ACGT", min_size=1, max_size=200), st.integers(min_value=1, max_value=50))
def test_gc_content_window_count_is_correct(seq: str, window: int) -> None:
    """The number of windows must match the sequence length and the stride."""
    size = min(window, len(seq))
    windows = gc_content(seq, window=size)
    assert len(windows) == len(seq) - size + 1


@given(st.text(alphabet=AMINO_ACIDS, min_size=1, max_size=120))
def test_translate_back_translate_round_trip_is_lossless(protein: str) -> None:
    """Any protein must survive a back-translation round trip unchanged."""
    cds = back_translate(protein, full_usage_table())
    assert translate(cds) == protein
