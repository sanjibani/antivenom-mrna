"""The cocktail planner: what a Big Four antibody cocktail would take.

M4 mapped the candidates. M5 asks the question that map was built to answer:
if you had to build one cocktail that covers India's Big Four snakes, which
payloads would it contain, and what would stop you?

The answer is computed rather than asserted. The planner treats the problem
as a set cover: each payload covers a set of cells, where a cell is a
species and a toxin family, and the planner finds the smallest set of
payloads that covers as much as the published record allows. The tiers then
tighten the question, one constraint at a time, and the attrition between
them is the finding.

Three things make the plan honest rather than tidy.

The first is the excluded cells. The Kunitz-type serine protease
inhibitors are the anticoagulant fraction of viper venom, so an antibody
against them would make a viper bite's coagulopathy worse. Those cells are
not targets, and the planner says so rather than reporting them as gaps
that somebody ought to fill.

The second is the payload rule. A cocktail is a set of antibodies, not a set
of rows in a table, so an antibody that neutralises one toxin family in two
species is one payload that covers two cells and appears in the plan once.

The third is the minimum-cover rule. The planner never proposes a payload
that covers a cell another payload in the same plan already covers, because
a redundant payload costs a second construct for nothing.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass
from itertools import combinations

from mrna.candidates import (
    BIG_FOUR,
    CLINICAL_SET,
    EVIDENCE_LEVELS,
    SINGLE_ORF_FORMATS,
    SPECIES_TOXIN_FAMILIES,
    Candidate,
)
from mrna.tracker import COVERAGE_THRESHOLD, CandidateTracker

__all__ = [
    "EXCLUDED_FAMILIES",
    "MAXIMUM_SEARCH_CANDIDATES",
    "TIERS",
    "Cell",
    "CocktailPlan",
    "PayloadGroup",
    "Tier",
    "excluded_cells",
    "plan_all_tiers",
    "plan_cocktail",
    "target_cells",
]

#: Toxin families the planner will not try to cover, and why.
#:
#: The Kunitz-type serine protease inhibitors are the anticoagulant fraction
#: of viper venom: they bind the clotting factors directly, so an antibody
#: against them would protect the clotting cascade in a patient whose
#: problem is that the cascade is being consumed. A planner that treated
#: every empty cell as an opportunity would be pointing the next milestone
#: at a design that would make the patient worse off.
EXCLUDED_FAMILIES: dict[str, str] = {
    "KSPI": (
        "the anticoagulant fraction of viper venom; an antibody against it "
        "would make the coagulopathy worse rather than better"
    ),
}

#: The largest payload pool the exact cover search will accept.
#:
#: The search is exhaustive, so its cost is exponential in the number of
#: payloads that cover at least one cell. Sixteen is a few tens of thousands
#: of subsets, which is instant. Rather than quietly becoming inexact on a
#: bigger table, the planner refuses.
MAXIMUM_SEARCH_CANDIDATES = 16


@dataclass(frozen=True)
class Cell:
    """One cell of the coverage matrix: a species and a toxin family.

    Attributes:
        species: The binomial name of the snake.
        toxin_family: The toxin family.
    """

    species: str
    toxin_family: str

    def describe(self) -> str:
        """Return the cell as a readable label."""
        return f"{self.species} / {self.toxin_family}"


@dataclass(frozen=True)
class Tier:
    """One constraint the planner applies on top of the ones before it.

    Attributes:
        name: The tier name, used in the report and in the plan.
        question: The question the tier asks, in one sentence.
        allows: The predicate a candidate has to satisfy.
    """

    name: str
    question: str
    allows: Callable[[Candidate], bool]


def _has_binding_evidence(candidate: Candidate) -> bool:
    """Whether the candidate has reached the coverage threshold."""
    ladder = list(EVIDENCE_LEVELS)
    return ladder.index(candidate.evidence_level) >= ladder.index(COVERAGE_THRESHOLD)


def _sequence_is_locatable(candidate: Candidate) -> bool:
    """Whether the candidate's sequence can be located in a public database."""
    return candidate.mrna_suitability in {"ready", "adaptable"}


#: The four questions the planner asks, in order.
#:
#: Each tier is the previous tier plus one more constraint, so the attrition
#: down the table is the finding.
TIERS: tuple[Tier, ...] = (
    Tier(
        name="published",
        question="Is there a published candidate with binding evidence at all?",
        allows=_has_binding_evidence,
    ),
    Tier(
        name="sequence-available",
        question="Can its sequence be located in a public database?",
        allows=lambda candidate: _has_binding_evidence(candidate)
        and _sequence_is_locatable(candidate),
    ),
    Tier(
        name="single-polypeptide",
        question="Is it one polypeptide, so that one transcript can carry it?",
        allows=lambda candidate: _has_binding_evidence(candidate)
        and _sequence_is_locatable(candidate)
        and candidate.molecule_format in SINGLE_ORF_FORMATS,
    ),
    Tier(
        name="licence-unblocked",
        question="Has a public patent search found no application covering it?",
        allows=lambda candidate: _has_binding_evidence(candidate)
        and _sequence_is_locatable(candidate)
        and candidate.molecule_format in SINGLE_ORF_FORMATS
        and candidate.patent_status == "none-located",
    ),
)


def target_cells(species: Iterable[str] = CLINICAL_SET) -> tuple[Cell, ...]:
    """Return the cells the planner tries to cover.

    Args:
        species: The species to plan for. Defaults to the Big Four plus the
            monocled cobra.

    Returns:
        The cells, ordered by species and then by toxin family.
    """
    cells: list[Cell] = []
    for name in species:
        for family in SPECIES_TOXIN_FAMILIES[name]:
            if family in EXCLUDED_FAMILIES:
                continue
            cells.append(Cell(name, family))
    return tuple(cells)


def excluded_cells(species: Iterable[str] = CLINICAL_SET) -> tuple[Cell, ...]:
    """Return the cells the planner deliberately does not try to cover.

    Args:
        species: The species to plan for. Defaults to the Big Four plus the
            monocled cobra.

    Returns:
        The cells, ordered by species and then by toxin family.
    """
    cells: list[Cell] = []
    for name in species:
        for family in SPECIES_TOXIN_FAMILIES[name]:
            if family in EXCLUDED_FAMILIES:
                cells.append(Cell(name, family))
    return tuple(cells)


def _candidate_cells(candidate: Candidate) -> frozenset[Cell]:
    """Return the cells a candidate row covers.

    A row covers a cell when it names the species and the toxin family. A
    row whose target toxin family is ``mixed`` or ``unspecified`` covers
    nothing, because the published record does not resolve which toxins it
    acts on, and a planner cannot count what it cannot place.

    Args:
        candidate: The candidate row to place.

    Returns:
        The cells the row covers.
    """
    return frozenset(
        Cell(species, candidate.toxin_family) for species in candidate.species
    )


@dataclass(frozen=True)
class PayloadGroup:
    """One antibody, with every cell its rows cover.

    A cocktail is a set of antibodies, not a set of rows. An antibody that
    neutralises the same toxin family in two species is one payload that
    covers two cells, and it belongs in the plan once.

    Attributes:
        antibody_id: The clone name as the source publishes it.
        cells: The cells the antibody covers, across every row it has.
        rows: The rows of the table that describe it.
    """

    antibody_id: str
    cells: frozenset[Cell]
    rows: tuple[Candidate, ...]

    @property
    def molecule_format(self) -> str:
        """The antibody format of the first row describing this antibody."""
        return self.rows[0].molecule_format

    @property
    def mrna_suitability(self) -> str:
        """The mRNA suitability recorded for this antibody."""
        return self.rows[0].mrna_suitability

    @property
    def patent_status(self) -> str:
        """The patent status recorded for this antibody."""
        return self.rows[0].patent_status

    @property
    def species(self) -> tuple[str, ...]:
        """Every species this antibody is recorded against."""
        names: list[str] = []
        for row in self.rows:
            for species in row.species:
                if species not in names:
                    names.append(species)
        return tuple(names)


def _group_candidates(candidates: list[Candidate]) -> list[PayloadGroup]:
    """Merge the rows of the same antibody into one payload.

    Args:
        candidates: The rows to group, in table order.

    Returns:
        One group per antibody, in the order the antibodies first appear.
    """
    rows_by_id: dict[str, list[Candidate]] = {}
    for candidate in candidates:
        rows_by_id.setdefault(candidate.antibody_id, []).append(candidate)
    groups: list[PayloadGroup] = []
    for antibody_id, rows in rows_by_id.items():
        cells: frozenset[Cell] = frozenset()
        for row in rows:
            cells |= _candidate_cells(row)
        groups.append(
            PayloadGroup(antibody_id=antibody_id, cells=cells, rows=tuple(rows)),
        )
    return groups


def _minimum_cover(
    groups: list[PayloadGroup], cells: tuple[Cell, ...],
) -> tuple[PayloadGroup, ...]:
    """Find the smallest set of payloads covering as many cells as possible.

    The search is exhaustive, so the answer is exact rather than greedy, and
    the tie is broken by table order so the plan is reproducible.

    Args:
        groups: The payloads the tier allows.
        cells: The cells to cover.

    Returns:
        The chosen payloads, in table order.

    Raises:
        ValueError: If the search space is too large to enumerate exactly.
    """
    if not cells:
        return ()
    index = {cell: position for position, cell in enumerate(cells)}
    masks: list[tuple[PayloadGroup, int]] = []
    for group in groups:
        mask = 0
        for cell in group.cells:
            position = index.get(cell)
            if position is not None:
                mask |= 1 << position
        if mask:
            masks.append((group, mask))
    if not masks:
        return ()
    if len(masks) > MAXIMUM_SEARCH_CANDIDATES:
        msg = (
            f"{len(masks)} payloads cover a cell, which is more than the "
            f"exact search accepts ({MAXIMUM_SEARCH_CANDIDATES}); the cover "
            "algorithm has to be replaced rather than made approximate"
        )
        raise ValueError(msg)

    target = (1 << len(cells)) - 1
    for size in range(1, len(masks) + 1):
        for combo in combinations(masks, size):
            union = 0
            for _, mask in combo:
                union |= mask
            if union == target:
                return tuple(group for group, _ in combo)

    best: tuple[PayloadGroup, ...] = ()
    best_count = 0
    for size in range(1, len(masks) + 1):
        for combo in combinations(masks, size):
            union = 0
            for _, mask in combo:
                union |= mask
            count = bin(union).count("1")
            if count > best_count:
                best_count = count
                best = tuple(group for group, _ in combo)
    return best


@dataclass(frozen=True)
class CocktailPlan:
    """The result of planning one tier.

    Attributes:
        tier: The name of the tier.
        question: The question the tier asked.
        payloads: The chosen payloads, in table order.
        covered: The cells the plan covers.
        uncovered: The target cells the plan does not cover.
        excluded: The cells the planner deliberately does not target.
    """

    tier: str
    question: str
    payloads: tuple[PayloadGroup, ...]
    covered: tuple[Cell, ...]
    uncovered: tuple[Cell, ...]
    excluded: tuple[Cell, ...]

    @property
    def payload_ids(self) -> tuple[str, ...]:
        """The names of the chosen payloads."""
        return tuple(group.antibody_id for group in self.payloads)

    def matrix(self) -> dict[Cell, bool]:
        """Return the coverage of every target cell.

        Returns:
            Mapping of cell to whether the plan covers it.
        """
        covered = set(self.covered)
        return {cell: cell in covered for cell in self.covered + self.uncovered}

    def species_covered(self) -> tuple[str, ...]:
        """Return the species the plan covers at least one cell of.

        Returns:
            The binomial names, in the order of the clinical set.
        """
        names = {cell.species for cell in self.covered}
        return tuple(species for species in CLINICAL_SET if species in names)

    def big_four_uncovered(self) -> tuple[Cell, ...]:
        """Return the uncovered cells that belong to a Big Four species.

        Returns:
            The cells, in the order of the target cells.
        """
        return tuple(cell for cell in self.uncovered if cell.species in BIG_FOUR)


def plan_cocktail(
    tracker: CandidateTracker,
    tier: str,
    species: Iterable[str] = CLINICAL_SET,
) -> CocktailPlan:
    """Plan one tier of the cocktail.

    Args:
        tracker: The validated candidate table.
        tier: The name of the tier, one of :data:`TIERS`.
        species: The species to plan for. Defaults to the Big Four plus the
            monocled cobra.

    Returns:
        The plan for that tier.

    Raises:
        ValueError: If the tier name is not known.
    """
    matches = [candidate for candidate in TIERS if candidate.name == tier]
    if not matches:
        known = ", ".join(candidate.name for candidate in TIERS)
        msg = f"unknown cocktail tier {tier!r}; known tiers: {known}"
        raise ValueError(msg)
    chosen_tier = matches[0]

    cells = target_cells(species)
    allowed = [
        candidate for candidate in tracker.candidates if chosen_tier.allows(candidate)
    ]
    chosen = _minimum_cover(_group_candidates(allowed), cells)
    chosen_cells: set[Cell] = set()
    for group in chosen:
        chosen_cells |= group.cells

    covered = tuple(cell for cell in cells if cell in chosen_cells)
    uncovered = tuple(cell for cell in cells if cell not in chosen_cells)
    return CocktailPlan(
        tier=chosen_tier.name,
        question=chosen_tier.question,
        payloads=chosen,
        covered=covered,
        uncovered=uncovered,
        excluded=excluded_cells(species),
    )


def plan_all_tiers(
    tracker: CandidateTracker,
    species: Iterable[str] = CLINICAL_SET,
) -> tuple[CocktailPlan, ...]:
    """Plan every tier, in order.

    Args:
        tracker: The validated candidate table.
        species: The species to plan for. Defaults to the Big Four plus the
            monocled cobra.

    Returns:
        One plan per tier.
    """
    return tuple(plan_cocktail(tracker, tier.name, species) for tier in TIERS)
