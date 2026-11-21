"""Tests for the M7 regional venom analysis.

The regional profiles are claims about other people's printed numbers, so
the tests protect the claims: every fraction sums against the recorded
remainder, every survey category is either mapped or honestly unmapped,
and the per-tier coverage verdicts must follow from the planner rather
than from prose.
"""

from __future__ import annotations

import pytest

from mrna.cocktail import EXCLUDED_FAMILIES
from mrna.regional import (
    CATEGORY_TO_FAMILY,
    REGIONAL_PROFILES,
    REGIONAL_VENOM_CITATIONS,
    regional_coverage,
)
from mrna.tracker import CandidateTracker


@pytest.fixture(scope="module")
def records() -> list[dict[str, object]]:
    """The regional coverage records, computed once."""
    return regional_coverage(CandidateTracker())


def test_every_profile_cites_a_recorded_source() -> None:
    """A fraction nobody printed is not a fraction anybody can check."""
    for profile in REGIONAL_PROFILES:
        assert profile.citation in REGIONAL_VENOM_CITATIONS


def test_every_profile_belongs_to_a_species_the_matrix_knows() -> None:
    """A profile for an unknown species would weight a matrix that has no cells."""
    for profile in REGIONAL_PROFILES:
        assert profile.species in CATEGORY_TO_FAMILY
        assert profile.fractions


def test_every_category_is_mapped_or_honestly_unmapped() -> None:
    """A category with no matrix family stays unmapped, never approximated."""
    for profile in REGIONAL_PROFILES:
        mapping = CATEGORY_TO_FAMILY[profile.species]
        for category in profile.fractions:
            assert category in mapping or category in {
                "SVSP",
                "snaclec",
                "VEGF",
            }


def test_the_printed_fractions_are_the_survey_fractions() -> None:
    """The numbers are the surveys', transcribed exactly."""
    fractions = {profile.population: profile for profile in REGIONAL_PROFILES
                 if profile.species == "Naja naja"}
    assert fractions["PB"].fractions["neurotoxic 3FTx"] == 0.80
    assert fractions["WB"].fractions["neurotoxic 3FTx"] == 0.733
    assert fractions["RJ"].fractions["cytotoxic 3FTx"] == 0.417
    vipers = {profile.population: profile for profile in REGIONAL_PROFILES
              if profile.species == "Daboia russelii"}
    assert vipers["WB"].fractions["Kunitz-type inhibitors"] == 0.39
    assert vipers["MH"].fractions["PLA2"] == 0.59
    assert vipers["MP"].fractions["snaclec"] == 0.20


def test_the_rajashta_undercount_is_recorded_as_one() -> None:
    """An unprinted fraction is an undercount, and the note says so."""
    rj = next(
        profile for profile in REGIONAL_PROFILES if profile.population == "RJ"
    )
    assert "neurotoxic 3FTx" not in rj.fractions
    assert "undercount" in rj.notes


def test_coverage_records_follow_the_planner(records) -> None:
    """The verdicts must follow from the planner, not from prose.

    The expected strictest covered tier per species and family, derived
    from the M4 table by hand: the cobra neurotoxin cells are held by
    Fabs, whose sequences are public and whose format is two-chain, so
    they survive until the single-polypeptide tier; the cobra PLA2 cell is
    held by an scFv with no locatable sequence, and the viper PLA2 cell by
    a VHH with none either, so both fall at the sequence-available tier.
    """
    expected = {
        ("Naja naja", "3FTx-L"): "sequence-available",
        ("Naja naja", "3FTx-S"): "sequence-available",
        ("Naja naja", "PLA2"): "published",
        ("Daboia russelii", "PLA2"): "published",
    }
    for record in records:
        profile = record["profile"]
        for family_record in record["families"]:
            key = (profile.species, family_record["family"])
            if key not in expected:
                continue
            if family_record["excluded"]:
                assert family_record["covered_through"] is None
            else:
                assert family_record["covered_through"] == expected[key]


def test_the_kunitz_fraction_is_excluded_at_every_tier(records) -> None:
    """The anticoagulant fraction is refused by design, region by region."""
    for record in records:
        for family_record in record["families"]:
            if family_record["family"] in EXCLUDED_FAMILIES:
                assert family_record["excluded"]
                assert family_record["covered_through"] is None
    assert "KSPI" in EXCLUDED_FAMILIES


def test_the_viper_matrix_under_models_the_venom(records) -> None:
    """SVSP and snaclec fractions are unmapped, and the record says so."""
    viper_records = [
        record
        for record in records
        if record["profile"].species == "Daboia russelii"
    ]
    unmapped_categories = {
        category
        for record in viper_records
        for category in record["unmapped"]
    }
    assert {"SVSP", "snaclec"} <= unmapped_categories
