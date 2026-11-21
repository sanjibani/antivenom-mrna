"""Regional venom variation as a design input, not a footnote.

Biological background
---------------------
A species-level coverage matrix treats a species' venom as one thing. It is
not. The Indian spectacled cobra of Punjab carries about 80 per cent
neurotoxic three-finger toxins by venom mass; the desert population of
Rajasthan spreads its venom across cytotoxins and phospholipase A2
instead, and the commercial antivenom fails against it completely at five
times the lethal dose. The Russell's viper of the Gangetic Plain carries
about 39 per cent Kunitz-type inhibitors by mass, the Western Ghats
population 59 per cent phospholipase A2. The same species, different
venoms, and a cocktail sized for the average is sized for neither.

This module carries the printed regional composition data of the two
Sunagar-group surveys, both open access, and does one thing with it: it
asks, for each population, how much of the venom's mass falls in matrix
cells the cocktail plan can and cannot cover, and which of the remainder
the matrix does not even model.

Two rules govern the mapping, and both are stated rather than smoothed
over.

First, the surveys' toxin categories are not the matrix's toxin families.
The cobra survey's neurotoxic three-finger toxins include both long-chain
and short-chain alpha-neurotoxins, while the matrix separates 3FTx-L, the
long-chain neurotoxins, from 3FTx-S, the cytotoxins. The mapping used here
puts the survey's neurotoxic fraction against the matrix's 3FTx-L cell and
its cytotoxic fraction against 3FTx-S, which is the closest available
correspondence and an approximation the report names.

Second, a survey category with no matrix family, such as the snake venom
serine proteases of the viper survey, is recorded as unmapped rather than
squeezed into the nearest family. The unmapped fraction is itself a
finding: it is the part of the venom the matrix does not model, and the
part a reader would otherwise assume away.

Sources: Senji Laxme RR, Attarde S, Khochare S, et al. 2021. Biogeographical
venom variation in the Indian spectacled cobra (Naja naja) underscores the
pressing need for pan-India efficacious snakebite therapy. PLoS Negl Trop
Dis 15(2):e0009150. doi:10.1371/journal.pntd.0009150. Senji Laxme RR,
Khochare S, Attarde S, et al. 2021. Biogeographic venom variation in
Russell's viper (Daboia russelii) and the preclinical inefficacy of
antivenom therapy in snakebite hotspots. PLoS Negl Trop Dis 15(3):e0009247.
doi:10.1371/journal.pntd.0009247.
"""

from __future__ import annotations

from dataclasses import dataclass

from mrna.cocktail import (
    EXCLUDED_FAMILIES,
    TIERS,
    plan_cocktail,
)
from mrna.tracker import CandidateTracker

__all__ = [
    "CATEGORY_TO_FAMILY",
    "REGIONAL_PROFILES",
    "REGIONAL_VENOM_CITATIONS",
    "RegionalProfile",
    "regional_coverage",
]


@dataclass(frozen=True)
class RegionalProfile:
    """One population's venom composition, as the survey printed it.

    Attributes:
        species: The matrix species name the population belongs to.
        population: The survey's population code.
        state: The Indian state the venom came from.
        zone: The biogeographic zone the survey assigns it.
        fractions: The toxin-family mass fractions the survey prints, keyed
            by the survey's own category names, exactly as printed. Only
            families the survey quantifies appear here; the remainder of
            the venom is the unmapped tail and is never invented.
        citation: Key into :data:`REGIONAL_VENOM_CITATIONS`.
        notes: What a reader needs in order to not be misled.
    """

    species: str
    population: str
    state: str
    zone: str
    fractions: dict[str, float]
    citation: str
    notes: str = ""

    @property
    def mapped_total(self) -> float:
        """The total venom mass the profile assigns to matrix families."""
        return sum(_mapped_fractions(self).values())

    @property
    def unmapped_total(self) -> float:
        """The venom mass the matrix does not model."""
        return sum(
            fraction
            for category, fraction in self.fractions.items()
            if category not in CATEGORY_TO_FAMILY.get(self.species, {})
        )


#: The survey's toxin categories that correspond to a matrix family, per
#: species. A category absent here is unmapped, and stays unmapped.
CATEGORY_TO_FAMILY_COBRA: dict[str, str] = {
    "neurotoxic 3FTx": "3FTx-L",
    "cytotoxic 3FTx": "3FTx-S",
    "PLA2": "PLA2",
}

CATEGORY_TO_FAMILY_VIPER: dict[str, str] = {
    "PLA2": "PLA2",
    "Kunitz-type inhibitors": "KSPI",
}

CATEGORY_TO_FAMILY: dict[str, dict[str, str]] = {
    "Naja naja": CATEGORY_TO_FAMILY_COBRA,
    "Daboia russelii": CATEGORY_TO_FAMILY_VIPER,
}


def _mapped_fractions(profile: RegionalProfile) -> dict[str, float]:
    """The profile's fractions that map onto matrix families.

    Args:
        profile: The regional profile.

    Returns:
        Mapping of matrix family to venom mass fraction.
    """
    mapping = CATEGORY_TO_FAMILY[profile.species]
    mapped: dict[str, float] = {}
    for category, fraction in profile.fractions.items():
        family = mapping.get(category)
        if family is not None:
            mapped[family] = mapped.get(family, 0.0) + fraction
    return mapped


#: The printed regional composition data of the two surveys.
#:
#: Only mass fractions the surveys print in their results text appear here,
#: one entry per surveyed population with proteomics. The cobra survey
#: quantified three populations by mass spectrometry, the viper survey
#: three; the other populations of both surveys were characterised by
#: chromatography alone and have no printed composition to record.
REGIONAL_PROFILES: tuple[RegionalProfile, ...] = (
    RegionalProfile(
        species="Naja naja",
        population="PB",
        state="Punjab",
        zone="North, semi-arid",
        fractions={
            "neurotoxic 3FTx": 0.80,
            "cytotoxic 3FTx": 0.10,
        },
        citation="Laxme2021Cobra",
        notes=(
            "The most neurotoxic cobra venom of the survey. The PLA2 "
            "fraction is inside the survey's 0.04 to 20 per cent range but "
            "the per-population value is printed only in the supplementary "
            "tables, so it is not recorded here."
        ),
    ),
    RegionalProfile(
        species="Naja naja",
        population="WB",
        state="West Bengal",
        zone="East, Gangetic Plain",
        fractions={
            "neurotoxic 3FTx": 0.733,
            "cytotoxic 3FTx": 0.236,
        },
        citation="Laxme2021Cobra",
        notes="",
    ),
    RegionalProfile(
        species="Naja naja",
        population="RJ",
        state="Rajasthan",
        zone="Desert",
        fractions={
            "cytotoxic 3FTx": 0.417,
            "PLA2": 0.20,
        },
        citation="Laxme2021Cobra",
        notes=(
            "The venom the commercial antivenom failed against completely "
            "at five times the lethal dose. Neurotoxic 3FTx was reported "
            "only as about 30 per cent, which is not a printed fraction "
            "and is therefore not recorded; the mapped total of this "
            "profile is accordingly an undercount, not an overcount."
        ),
    ),
    RegionalProfile(
        species="Daboia russelii",
        population="WB",
        state="West Bengal",
        zone="East, Gangetic Plain",
        fractions={
            "Kunitz-type inhibitors": 0.39,
            "PLA2": 0.35,
            "SVSP": 0.08,
            "snaclec": 0.05,
        },
        citation="Laxme2021Viper",
        notes=(
            "The Kunitz-type fraction is the anticoagulant fraction the "
            "planner excludes on purpose. The serine protease and C-type "
            "lectin fractions have no matrix family and stay unmapped."
        ),
    ),
    RegionalProfile(
        species="Daboia russelii",
        population="MH",
        state="Maharashtra",
        zone="Southwest, Western Ghats",
        fractions={
            "PLA2": 0.59,
            "SVSP": 0.27,
            "Kunitz-type inhibitors": 0.04,
            "snaclec": 0.02,
        },
        citation="Laxme2021Viper",
        notes="",
    ),
    RegionalProfile(
        species="Daboia russelii",
        population="MP",
        state="Madhya Pradesh",
        zone="Central, Deccan Plateau",
        fractions={
            "PLA2": 0.45,
            "snaclec": 0.20,
            "SVSP": 0.13,
            "Kunitz-type inhibitors": 0.07,
            "VEGF": 0.08,
        },
        citation="Laxme2021Viper",
        notes=(
            "The snaclec and VEGF fractions have no matrix family. The "
            "printed snaclec value is 'nearly about 20 per cent', recorded "
            "here as 0.20."
        ),
    ),
)

#: Primary sources for the regional profiles.
REGIONAL_VENOM_CITATIONS: dict[str, str] = {
    "Laxme2021Cobra": (
        "Senji Laxme RR, Attarde S, Khochare S, et al. 2021. Biogeographical "
        "venom variation in the Indian spectacled cobra (Naja naja) "
        "underscores the pressing need for pan-India efficacious snakebite "
        "therapy. PLoS Negl Trop Dis 15(2):e0009150. "
        "doi:10.1371/journal.pntd.0009150."
    ),
    "Laxme2021Viper": (
        "Senji Laxme RR, Khochare S, Attarde S, et al. 2021. Biogeographic "
        "venom variation in Russell's viper (Daboia russelii) and the "
        "preclinical inefficacy of antivenom therapy in snakebite hotspots. "
        "PLoS Negl Trop Dis 15(3):e0009247. "
        "doi:10.1371/journal.pntd.0009247."
    ),
}


def regional_coverage(tracker: CandidateTracker) -> list[dict[str, object]]:
    """Ask, per population, what the plan can and cannot touch, tier by tier.

    The question a species-level matrix cannot ask. A cell is covered or
    uncovered at the species level; a population's venom mass falls into
    cells family by family, and the same cell can be covered at the
    published tier by a two-chain Fab that no single transcript can carry.
    A family's status is therefore computed at every tier, in order, and
    the record carries the strictest tier at which the family is still
    covered.

    Args:
        tracker: The validated candidate table.

    Returns:
        One record per regional profile, in profile order. Each record
        carries the profile, and per mapped matrix family: the venom mass
        fraction behind it, the strictest tier at which the plan covers
        it, whether the planner excludes it on purpose, and the reason it
        goes uncovered where it does. Survey categories with no matrix
        family are recorded as unmapped rather than squeezed into the
        nearest family.
    """
    covered_by_tier: dict[str, dict[str, set[str]]] = {}
    for tier in TIERS:
        plan = plan_cocktail(tracker, tier.name)
        covered: dict[str, set[str]] = {}
        for group in plan.payloads:
            for cell in group.cells:
                covered.setdefault(cell.species, set()).add(cell.toxin_family)
        covered_by_tier[tier.name] = covered

    records: list[dict[str, object]] = []
    for profile in REGIONAL_PROFILES:
        mapping = CATEGORY_TO_FAMILY[profile.species]
        mapped = _mapped_fractions(profile)
        family_records: list[dict[str, object]] = []
        for family, fraction in mapped.items():
            excluded = family in EXCLUDED_FAMILIES
            covered_through = None
            for tier in TIERS:
                if family in covered_by_tier[tier.name].get(
                    profile.species,
                    set(),
                ):
                    covered_through = tier.name
            family_records.append(
                {
                    "family": family,
                    "fraction": fraction,
                    "covered_through": covered_through,
                    "excluded": excluded,
                },
            )
        unmapped = {
            category: fraction
            for category, fraction in profile.fractions.items()
            if category not in mapping
        }
        records.append(
            {
                "profile": profile,
                "families": family_records,
                "unmapped": unmapped,
            },
        )
    return records
