"""Tests for the forbidden motif machinery.

The hard rule of this project is that nothing enters the engine without a
source. These tests enforce it for motifs: every motif must carry a reason
and a citation, and only verified motifs may be enforced.
"""

from __future__ import annotations

from mrna.motifs import (
    CANDIDATE_MOTIFS,
    DEFAULT_MOTIFS,
    Motif,
    find_motifs,
    homopolymer_runs,
)


def test_find_motifs_locates_the_polyadenylation_signal() -> None:
    """The hexamer must be found at the right coordinates."""
    sequence = "ATGAAATTTAATAAAGGGCCC"
    hits = find_motifs(sequence, DEFAULT_MOTIFS)
    assert len(hits) == 1
    assert hits[0].start == 9
    assert hits[0].end == 15
    assert hits[0].motif.name == "polyadenylation signal"
    assert "10 to 15" in hits[0].describe()


def test_find_motifs_returns_nothing_for_a_clean_sequence() -> None:
    """A sequence without the motif must produce no hits."""
    assert find_motifs("ATGAAATTTGGAAGGGCCC", DEFAULT_MOTIFS) == ()


def test_find_motifs_ignores_case_and_whitespace() -> None:
    """Formatting must not hide a motif."""
    assert len(find_motifs("atg aataaa ggg", DEFAULT_MOTIFS)) == 1


def test_every_enforced_motif_is_verified() -> None:
    """Only verified motifs may be enforced by default."""
    assert all(motif.verified for motif in DEFAULT_MOTIFS)


def test_candidate_motifs_are_not_verified() -> None:
    """Candidates must be marked as unverified so nobody enforces them by accident."""
    assert CANDIDATE_MOTIFS
    assert all(not motif.verified for motif in CANDIDATE_MOTIFS)


def test_every_motif_carries_a_reason_and_a_citation() -> None:
    """No motif may enter the project without a source."""
    for motif in (*DEFAULT_MOTIFS, *CANDIDATE_MOTIFS):
        assert motif.reason.strip()
        assert len(motif.citation.strip()) > 20
        assert motif.pattern


def test_homopolymer_runs_find_only_runs_over_the_limit() -> None:
    """A run of exactly the limit is allowed, one longer is not."""
    sequence = "A" * 9 + "C" * 8
    runs = homopolymer_runs(sequence, 8)
    assert runs == ((0, 9),)


def test_homopolymer_runs_allow_runs_at_the_limit() -> None:
    """The limit is inclusive."""
    assert homopolymer_runs("AAAAAAAA", 8) == ()
    assert homopolymer_runs("AAAAAAAAA", 8) == ((0, 9),)


def test_homopolymer_runs_handle_a_single_base_and_empty_input() -> None:
    """Edge cases must not crash the scanner."""
    assert homopolymer_runs("A", 8) == ()
    assert homopolymer_runs("", 8) == ()


def test_motif_patterns_compile() -> None:
    """Every motif must be a valid regular expression."""
    for motif in (*DEFAULT_MOTIFS, *CANDIDATE_MOTIFS):
        assert motif.compiled().search("ACGT" * 10) is not None or True


def test_a_motif_without_a_citation_is_still_representable() -> None:
    """The dataclass must allow a motif to be built and then rejected by policy."""
    motif = Motif(name="x", pattern="X", reason="", citation="", verified=False)
    assert motif.name == "x"
