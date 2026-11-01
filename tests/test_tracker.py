"""Tests for the M4 candidate tracker.

The tracker is a claim about the published literature, so the tests protect
the claims rather than the code: every row must use the controlled
vocabularies, the declared mRNA suitability must match what the format and
the sequence source actually imply, and the gap analysis must find the gaps
it says it finds.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mrna.candidates import (
    BIG_FOUR,
    CANDIDATES,
    CITATIONS,
    EVIDENCE_LEVELS,
    FORMATS,
    MRNA_SUITABILITY,
    PATENT_STATUS,
    POLYCLONAL_FORMATS,
    SPECIES_TOXIN_FAMILIES,
    TOXIN_FAMILIES,
    TOXIN_FAMILY_NOTES,
    Candidate,
)
from mrna.tracker import (
    COVERAGE_THRESHOLD,
    CandidateTracker,
    derive_mrna_suitability,
    parse_sequence_source,
)

ROOT = Path(__file__).resolve().parent.parent
TSV_PATH = ROOT / "data" / "candidates.tsv"
REPORT_PATH = ROOT / "docs" / "m4_candidate_tracker.md"


def test_the_curated_table_validates() -> None:
    """The shipped table must pass its own validator."""
    CandidateTracker()


def test_every_row_uses_only_known_vocabulary_values() -> None:
    """No row may invent a value outside the controlled vocabularies."""
    for candidate in CANDIDATES:
        assert candidate.toxin_family in TOXIN_FAMILIES
        assert candidate.molecule_format in FORMATS
        assert candidate.patent_status in PATENT_STATUS
        assert candidate.evidence_level in EVIDENCE_LEVELS
        assert candidate.mrna_suitability in MRNA_SUITABILITY
        assert candidate.citation in CITATIONS


def test_every_citation_key_resolves_to_a_real_citation() -> None:
    """A citation key that resolves to nothing is an unusable row."""
    for candidate in CANDIDATES:
        assert len(CITATIONS[candidate.citation]) > 40


def test_a_polyclonal_candidate_is_blocked() -> None:
    """There is no defined sequence to transcribe from a polyclonal."""
    for fmt in POLYCLONAL_FORMATS:
        assert derive_mrna_suitability(fmt, ("none",)) == "blocked"


def test_a_single_chain_with_a_public_sequence_is_ready() -> None:
    """One open reading frame and a locatable sequence means ready to transcribe."""
    assert derive_mrna_suitability("scFv", ("PDB:8IA6",)) == "ready"
    assert derive_mrna_suitability("VHH", ("PDB:9GCN",)) == "ready"
    assert derive_mrna_suitability("scFv", ("paper",)) == "ready"


def test_a_two_chain_format_with_a_public_sequence_is_adaptable() -> None:
    """A Fab or an IgG needs more than one open reading frame."""
    assert derive_mrna_suitability("Fab", ("PDB:9CLP",)) == "adaptable"
    assert derive_mrna_suitability("mAb-IgG", ("PDB:8V13",)) == "adaptable"


def test_a_monoclonal_with_no_located_sequence_is_unknown() -> None:
    """Without a sequence, nothing can be said about buildability."""
    assert derive_mrna_suitability("scFv", ("none",)) == "unknown"
    assert derive_mrna_suitability("VHH", ("patent:WO1",)) == "unknown"


def test_a_mismatched_suitability_is_rejected() -> None:
    """The declared column must not disagree with the derived one."""
    bad = Candidate(
        species=("Echis carinatus",),
        toxin_family="PLA2",
        toxin_target="ecarpholin S",
        antibody_id="E113",
        molecule_format="scFv",
        sequence_source=("PDB:8IA6",),
        patent_status="none-located",
        tested_against="in vivo",
        evidence_level="L4",
        mrna_suitability="unknown",
        citation="Kumar2024",
    )
    with pytest.raises(ValueError, match="disagrees with the derived value"):
        CandidateTracker((bad,))


def test_an_unknown_species_is_rejected() -> None:
    """A row naming a species outside the vocabulary cannot be trusted."""
    bad = Candidate(
        species=("Python molurus",),
        toxin_family="PLA2",
        toxin_target="x",
        antibody_id="X",
        molecule_format="scFv",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against="x",
        evidence_level="L2",
        mrna_suitability="unknown",
        citation="Kumar2024",
    )
    with pytest.raises(ValueError, match="unknown species"):
        CandidateTracker((bad,))


def test_a_malformed_sequence_source_is_rejected() -> None:
    """Sequence source tokens must be kind:value pairs."""
    with pytest.raises(ValueError, match="not of the form kind:value"):
        parse_sequence_source(("8IA6",))
    with pytest.raises(ValueError, match="unknown kind"):
        parse_sequence_source(("pdb:8IA6",))


def test_every_big_four_species_has_at_least_one_candidate() -> None:
    """The map would be useless if it skipped a species."""
    named = {species for candidate in CANDIDATES for species in candidate.species}
    for species in BIG_FOUR:
        assert species in named


def test_the_cobra_phospholipase_cell_is_covered() -> None:
    """Naja naja PLA2 is covered, by the scFvs from the Salunke laboratory."""
    tracker = CandidateTracker()
    covering = tracker.covering("Naja naja", "PLA2")
    assert {c.antibody_id for c in covering} == {"N194"}


def test_the_russells_viper_metalloproteinase_cell_is_empty() -> None:
    """No published antibody candidate against Daboia SVMP was located."""
    tracker = CandidateTracker()
    assert tracker.covering("Daboia russelii", "SVMP") == []
    assert any(
        gap.species == "Daboia russelii" and gap.toxin_family == "SVMP"
        for gap in tracker.gaps()
    )


def test_the_krait_phospholipase_cell_is_empty() -> None:
    """Krait venom is mostly PLA2 and no candidate targets it."""
    tracker = CandidateTracker()
    assert tracker.covering("Bungarus caeruleus", "PLA2") == []


def test_the_gap_analysis_reports_every_uncovered_cell() -> None:
    """The gap list and the coverage matrix must agree."""
    tracker = CandidateTracker()
    expected = [
        (species, family)
        for species, families in SPECIES_TOXIN_FAMILIES.items()
        for family in families
        if not tracker.covering(species, family)
    ]
    actual = [(gap.species, gap.toxin_family) for gap in tracker.gaps()]
    assert sorted(actual) == sorted(expected)


def test_a_gap_reports_its_nearest_coverage() -> None:
    """A gap is more useful when it says what is nearby."""
    tracker = CandidateTracker()
    gap = next(
        g for g in tracker.gaps() if g.species == "Daboia russelii" and g.toxin_family == "SVMP"
    )
    assert "PLA2" in gap.near_family


def test_every_candidate_names_a_species_with_a_known_toxin_profile() -> None:
    """A row for a species outside the toxin profile has nowhere to sit."""
    for candidate in CANDIDATES:
        for species in candidate.species:
            if candidate.toxin_family in {"mixed", "unspecified"}:
                continue
            assert candidate.toxin_family in SPECIES_TOXIN_FAMILIES[species]


def test_the_coverage_threshold_is_a_real_rung() -> None:
    """The threshold has to be a rung on the ladder, not a free-floating value."""
    assert COVERAGE_THRESHOLD in EVIDENCE_LEVELS


def test_the_tsv_has_a_header_and_one_row_per_candidate() -> None:
    """The published table must be machine-readable and complete."""
    tracker = CandidateTracker()
    lines = tracker.to_tsv().strip().split("\n")
    assert len(lines) == len(CANDIDATES) + 1
    header = lines[0].split("\t")
    assert header[:10] == [
        "species",
        "toxin_family",
        "toxin_target",
        "antibody_id",
        "format",
        "sequence_source",
        "patent_status",
        "tested_against",
        "evidence_level",
        "mrna_suitability",
    ]
    for line in lines[1:]:
        assert len(line.split("\t")) == len(header)


def test_the_report_carries_the_disclaimer_and_every_candidate() -> None:
    """The report is an artifact, so it carries the artifact rules."""
    report = CandidateTracker().render_markdown()
    assert "No therapeutic efficacy demonstrated. Not medical advice." in report
    for candidate in CANDIDATES:
        assert candidate.antibody_id in report


def test_the_findings_are_computed_from_the_table() -> None:
    """The headline findings must follow from the rows, not from prose."""
    tracker = CandidateTracker()
    findings = tracker.findings()
    public = len(tracker.public_sequence_candidates())
    assert f"Only {public} of {len(CANDIDATES)}" in findings[0]


def test_the_report_marks_the_big_four_in_the_gap_matrix() -> None:
    """A reader has to know which gaps the project is actually about."""
    report = CandidateTracker().render_markdown()
    for species in BIG_FOUR:
        assert f"**{species}**" in report


def test_the_report_explains_the_kspi_cells() -> None:
    """An empty cell that is a correct decision must be explained."""
    report = CandidateTracker().render_markdown()
    assert "An empty cell is not always an opportunity" in report
    assert "anticoagulant fraction" in report


def test_every_toxin_family_the_table_uses_has_a_note() -> None:
    """A family with no explanation cannot be read in the matrix."""
    for candidate in CANDIDATES:
        assert candidate.toxin_family in TOXIN_FAMILY_NOTES


def test_the_report_states_what_it_does_not_claim() -> None:
    """A survey must say it is a survey."""
    report = CandidateTracker().render_markdown()
    assert "It is a survey, not an exhaustive census." in report
    assert "not a freedom-to-operate opinion" in report


def test_the_committed_table_matches_the_generated_one() -> None:
    """The published table cannot drift away from the code that validates it."""
    if not TSV_PATH.exists():
        pytest.skip("data/candidates.tsv has not been generated yet")
    assert TSV_PATH.read_text(encoding="utf-8") == CandidateTracker().to_tsv()


def test_the_committed_report_matches_the_generated_one() -> None:
    """The published report cannot drift away from the code that generates it."""
    if not REPORT_PATH.exists():
        pytest.skip("docs/m4_candidate_tracker.md has not been generated yet")
    assert REPORT_PATH.read_text(encoding="utf-8") == CandidateTracker().render_markdown()


def test_the_artifacts_contain_no_em_dash() -> None:
    """The repository forbids em-dashes in committed text."""
    tracker = CandidateTracker()
    for text in (tracker.to_tsv(), tracker.render_markdown()):
        assert "\u2014" not in text
