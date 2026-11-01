"""Tests for the M5 cocktail planner and the constructs it licenses.

The planner is a claim about what the published record can supply, so the
tests protect the claim rather than the code. Three things have to hold.

First, the cover has to be exact. The planner searches exhaustively rather
than greedily, so the plan it returns must be minimal, and removing any
payload from it must lose coverage. A greedy search would satisfy neither
test as often as it should, so a plan that passes these tests is evidence
the search is doing what it says.

Second, the attrition down the tiers has to be real. Each tier adds a
constraint, so a cell that is covered at a permissive tier and uncovered at a
stricter one is a finding about the literature, not a bug. The most important
of these is that Indian Russell's viper drops out of the plan entirely once
the sequence has to be publicly locatable.

Third, the report must not state anything the data does not support. The
prose in ``docs/m5_big_four_cocktail.md`` is generated, and a generated
number that drifts away from the computed one is the same failure as a
hand-written one that was never checked. The report tests below recompute
the counts and compare them against the text.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

from mrna.assemble import A110_POLYA, full_construct
from mrna.candidates import (
    BIG_FOUR,
    CLINICAL_SET,
    SPECIES_FAMILIES,
    SPECIES_TOXIN_FAMILIES,
)
from mrna.cocktail import (
    EXCLUDED_FAMILIES,
    MAXIMUM_SEARCH_CANDIDATES,
    TIERS,
    Cell,
    PayloadGroup,
    excluded_cells,
    plan_all_tiers,
    plan_cocktail,
    target_cells,
)
from mrna.leader import HSA_SIGNAL_PEPTIDE
from mrna.motifs import DEFAULT_MOTIFS
from mrna.optimize import optimize_cds
from mrna.payload import PAYLOAD_SPECS, load_payloads
from mrna.policy import HUMAN_POLICY
from mrna.reference import load_codon_usage
from mrna.tracker import CandidateTracker
from mrna.utr import AES_MTRNR1_3UTR, ALPHA_GLOBIN_5UTR

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"
FASTA_PATH = ROOT / "data" / "designs" / "big_four_cocktail_mrna.fasta"
ELEMENTS_PATH = ROOT / "data" / "designs" / "big_four_cocktail.elements.tsv"
REPORT_PATH = ROOT / "docs" / "m5_big_four_cocktail.md"
M3_FASTA_PATH = ROOT / "data" / "designs" / "e113_scfv_ecarpholin_s_mrna.fasta"


def _load_script() -> object:
    """Import scripts/cocktail_design.py as a module.

    The script is not part of the package, so it is loaded by path. The
    point of importing it rather than shelling out is that the report and
    the design files can be compared against what the script would write,
    without writing anything.

    Returns:
        The imported module.
    """
    path = ROOT / "scripts" / "cocktail_design.py"
    spec = importlib.util.spec_from_file_location("cocktail_design", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def tracker() -> CandidateTracker:
    """The validated candidate table."""
    return CandidateTracker()


@pytest.fixture(scope="module")
def plans(tracker: CandidateTracker) -> tuple[object, ...]:
    """One plan per tier, in tier order."""
    return plan_all_tiers(tracker)


# --------------------------------------------------------------------- cells


def test_the_target_set_is_the_clinical_set() -> None:
    """The planner plans for the Big Four plus the monocled cobra."""
    expected = {
        Cell(species, family)
        for species in CLINICAL_SET
        for family in SPECIES_TOXIN_FAMILIES[species]
        if family not in EXCLUDED_FAMILIES
    }
    assert set(target_cells()) == expected
    assert len(target_cells()) == 14


def test_the_excluded_cells_are_the_kspi_cells_and_nothing_else() -> None:
    """Only the Kunitz-type serine protease inhibitor cells are excluded."""
    excluded = excluded_cells()
    assert excluded
    assert {cell.toxin_family for cell in excluded} == set(EXCLUDED_FAMILIES)
    for cell in excluded:
        assert cell.species in BIG_FOUR


def test_an_excluded_cell_is_never_a_target_cell() -> None:
    """The two sets must be disjoint, or the exclusion is meaningless."""
    targets = set(target_cells())
    for cell in excluded_cells():
        assert cell not in targets


def test_every_target_cell_belongs_to_a_species_the_planner_knows() -> None:
    """A cell whose species has no known toxin families would be a bug."""
    for cell in target_cells():
        assert cell.species in CLINICAL_SET
        assert SPECIES_FAMILIES[cell.species] in {"viper", "elapid"}


# --------------------------------------------------------------------- tiers


def test_the_tiers_are_ordered_from_permissive_to_strict() -> None:
    """Each tier must add a constraint, never remove one."""
    assert [tier.name for tier in TIERS] == [
        "published",
        "sequence-available",
        "single-polypeptide",
        "licence-unblocked",
    ]


def test_coverage_never_increases_as_the_tiers_get_stricter(
    plans: tuple[object, ...],
) -> None:
    """The attrition down the tiers is the finding, so it must be monotonic."""
    counts = [len(plan.covered) for plan in plans]
    assert counts == sorted(counts, reverse=True)


def test_the_published_tier_covers_more_than_the_strictest(
    plans: tuple[object, ...],
) -> None:
    """The tiers must actually differ, or there is no finding."""
    assert len(plans[0].covered) > len(plans[-1].covered)


def test_the_plan_at_a_strict_tier_is_also_valid_at_every_looser_one(
    plans: tuple[object, ...],
) -> None:
    """A stricter plan may not cover a cell a looser plan leaves uncovered.

    This is the check that a constraint is actually being applied. Without
    it, a tier whose predicate is wrong would still produce a plausible
    looking table.
    """
    for stricter, looser in zip(plans[1:], plans, strict=False):
        for cell in stricter.covered:
            assert cell in looser.covered


def test_the_licence_unblocked_tier_holds_only_one_payload(
    plans: tuple[object, ...],
) -> None:
    """One row of the table has had a patent search performed and found none."""
    final = plans[-1]
    assert final.tier == "licence-unblocked"
    assert final.payload_ids == ("E113",)


def test_the_single_polypeptide_tier_is_the_two_transcribable_payloads(
    plans: tuple[object, ...],
) -> None:
    """The two payloads that can be transcribed as they stand."""
    buildable = plans[2]
    assert buildable.tier == "single-polypeptide"
    assert set(buildable.payload_ids) == {"E113", "TPL1158_01_C09"}


def test_india_renness_viper_drops_out_when_the_sequence_has_to_be_public(
    plans: tuple[object, ...],
) -> None:
    """The single most actionable finding of M5.

    At the permissive tier the viper is covered by a recombinant candidate,
    but that candidate's sequence cannot be located in a public database, so
    the species leaves the plan entirely at the next tier down.
    """
    permissive, restricted = plans[0], plans[1]
    viper_cells = {
        cell for cell in permissive.covered if cell.species == "Daboia russelii"
    }
    assert viper_cells, "the permissive tier should cover the viper somewhere"
    assert not any(
        cell.species == "Daboia russelii" for cell in restricted.covered
    )


def test_the_krait_phospholipase_cell_has_no_candidate_at_any_tier(
    plans: tuple[object, ...],
) -> None:
    """Common krait PLA2 is about two thirds of that venom by mass."""
    krait = Cell("Bungarus caeruleus", "PLA2")
    for plan in plans:
        assert krait not in plan.covered


def test_no_plan_covers_a_cell_the_planner_excluded(
    plans: tuple[object, ...],
) -> None:
    """The KSPI cells must stay empty at every tier."""
    excluded = set(excluded_cells())
    for plan in plans:
        assert not excluded & set(plan.covered)
        assert set(plan.excluded) == excluded


# ------------------------------------------------------------------- the cover


def test_the_cover_is_minimal(tracker: CandidateTracker) -> None:
    """No payload in the plan may be dropped without losing coverage.

    A redundant payload costs a second construct for nothing, so a plan
    that contains one is not a minimum cover even if it happens to cover
    everything.
    """
    for plan in plan_all_tiers(tracker):
        covered = set(plan.covered)
        for index in range(len(plan.payloads)):
            without = set()
            for other_index, group in enumerate(plan.payloads):
                if other_index == index:
                    continue
                without |= group.cells
            assert without != covered, (
                f"{plan.tier}: {plan.payloads[index].antibody_id} is redundant"
            )


def test_the_cover_is_exhaustive_and_exact(tracker: CandidateTracker) -> None:
    """The search must find the smallest cover, not merely a cover.

    The reference implementation below is deliberately unlike the planner's:
    it computes the best coverage achievable at all, then asks what the
    smallest number of payloads achieves it. The planner instead looks for
    an exact cover first and only falls back on the best partial cover.
    Both must land on the same number, so a greedy planner would fail this
    test.
    """
    for tier in TIERS:
        plan = plan_cocktail(tracker, tier.name)
        cells = target_cells()
        index = {cell: position for position, cell in enumerate(cells)}
        pool = []
        for group in _allowed_groups(tracker, tier.name):
            mask = 0
            for cell in group.cells:
                position = index.get(cell)
                if position is not None:
                    mask |= 1 << position
            if mask:
                pool.append(mask)
        assert len(plan.payloads) == _smallest_size_reaching_best_cover(pool)


def _smallest_size_reaching_best_cover(pool: list[int]) -> int:
    """Return the fewest payloads that achieve the best coverage available.

    Args:
        pool: One bitmask per payload, over the bits of the target cells.

    Returns:
        The size of the smallest subset whose union covers as many cells as
        any subset can.
    """
    best_cover = 0
    for size in range(1, len(pool) + 1):
        for combo in _combinations(pool, size):
            union = 0
            for mask in combo:
                union |= mask
            best_cover = max(best_cover, bin(union).count("1"))
    for size in range(1, len(pool) + 1):
        for combo in _combinations(pool, size):
            union = 0
            for mask in combo:
                union |= mask
            if bin(union).count("1") == best_cover:
                return size
    return 0


def _allowed_groups(tracker: CandidateTracker, tier_name: str) -> list:
    """Return the payload groups a tier allows, without using the planner."""
    tier = next(t for t in TIERS if t.name == tier_name)
    allowed = [c for c in tracker.candidates if tier.allows(c)]
    groups: dict[str, PayloadGroup] = {}
    for candidate in allowed:
        group = groups.get(candidate.antibody_id)
        cells = frozenset(
            Cell(species, candidate.toxin_family) for species in candidate.species
        )
        if group is None:
            groups[candidate.antibody_id] = PayloadGroup(
                antibody_id=candidate.antibody_id,
                cells=cells,
                rows=(candidate,),
            )
        else:
            groups[candidate.antibody_id] = PayloadGroup(
                antibody_id=group.antibody_id,
                cells=group.cells | cells,
                rows=(*group.rows, candidate),
            )
    return list(groups.values())


def _combinations(items: list, size: int) -> list:
    """Return every combination of the given size, without importing itertools."""
    if size == 0:
        return [()]
    result = []
    for index in range(len(items) - size + 1):
        head = items[index]
        for tail in _combinations(items[index + 1 :], size - 1):
            result.append((head, *tail))
    return result


def test_an_unknown_tier_is_refused(tracker: CandidateTracker) -> None:
    """A typo in a tier name must fail loudly."""
    with pytest.raises(ValueError, match="unknown cocktail tier"):
        plan_cocktail(tracker, "not-a-tier")


def test_the_search_refuses_to_become_inexact(monkeypatch: pytest.MonkeyPatch) -> None:
    """Past the search limit the planner raises rather than degrading.

    The alternative, quietly switching to a greedy cover, would produce a
    plan that looks the same and is not minimal.
    """
    import mrna.cocktail as cocktail

    monkeypatch.setattr(cocktail, "MAXIMUM_SEARCH_CANDIDATES", 1)
    with pytest.raises(ValueError, match="exact search accepts"):
        plan_cocktail(CandidateTracker(), "published")


def test_the_search_limit_is_documented_as_a_limit() -> None:
    """The constant is a bound, not a tuning knob."""
    assert MAXIMUM_SEARCH_CANDIDATES >= 16


def test_a_cocktail_is_a_set_of_antibodies_not_a_set_of_rows(
    plans: tuple[object, ...],
) -> None:
    """Centi-LNX-D09 covers two species, and appears in the plan once."""
    plan = plans[0]
    assert plan.payload_ids.count("Centi-LNX-D09") == 1
    group = next(
        g for g in plan.payloads if g.antibody_id == "Centi-LNX-D09"
    )
    assert len(group.cells) == 2
    assert {cell.species for cell in group.cells} == {
        "Naja naja",
        "Bungarus caeruleus",
    }


# ------------------------------------------------------------------ payloads


def test_every_payload_loads() -> None:
    """Each payload's committed sequence must load."""
    payloads = load_payloads()
    assert len(payloads) == len(PAYLOAD_SPECS)
    assert {payload.name for payload in payloads} == {
        spec.name for spec in PAYLOAD_SPECS
    }


def test_the_two_single_chain_payloads_are_the_only_ones() -> None:
    """A single polypeptide is what one transcript can carry."""
    single = {p.name for p in load_payloads() if p.is_single_polypeptide}
    assert single == {"E113", "TPL1158_01_C09"}


def test_every_payload_has_the_chains_its_specification_names() -> None:
    """The chain count in the specification must match the loaded payload."""
    for spec in PAYLOAD_SPECS:
        payload = next(p for p in load_payloads() if p.name == spec.name)
        assert len(payload.chains) == len(spec.chain_labels)
        assert payload.chain_lengths() == {
            label: len(sequence) for label, sequence in payload.chains
        }


def test_the_published_payload_lengths_are_the_database_lengths() -> None:
    """The lengths are what the Protein Data Bank holds."""
    expected = {
        "E113": 239,
        "TPL1158_01_C09": 123,
        "Centi-LNX-D09": 451,
        "Centi-SNX-B03": 445,
        "H11": 434,
        "A01": 261,
    }
    for payload in load_payloads():
        assert payload.length == expected[payload.name]


def test_a_multi_chain_payload_does_not_pretend_to_be_one_string() -> None:
    """Concatenating two chains without a linker is not a molecule.

    The planner never builds a construct from one, and the property that
    says so is what stops that from changing by accident.
    """
    for payload in load_payloads():
        if payload.is_single_polypeptide:
            assert len(payload.chains) == 1
        else:
            assert len(payload.chains) == 2


# ----------------------------------------------------------------- constructs


def _usage_table() -> dict[str, float]:
    """Load the human codon usage table."""
    return load_codon_usage(REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv")


def _build(payload) -> tuple[object, object]:
    """Build the construct for one payload, the way the script does."""
    protein = HSA_SIGNAL_PEPTIDE.sequence + payload.sequence
    design = optimize_cds(protein, _usage_table(), HUMAN_POLICY)
    coding = design.cds + "TGA"
    construct = full_construct(
        payload.name,
        utr5=ALPHA_GLOBIN_5UTR,
        coding=coding,
        utr3=AES_MTRNR1_3UTR,
        poly_a=A110_POLYA,
        leader=HSA_SIGNAL_PEPTIDE,
    )
    return construct, design


def test_each_transcribable_payload_translates_back_to_its_published_sequence() -> None:
    """The construct must encode the published protein, not a variant of it."""
    for payload in load_payloads():
        if not payload.is_single_polypeptide:
            continue
        construct, _design = _build(payload)
        protein = construct.protein
        assert protein.startswith(HSA_SIGNAL_PEPTIDE.sequence)
        assert protein[len(HSA_SIGNAL_PEPTIDE.sequence) :] == payload.sequence


def test_the_transcribable_constructs_satisfy_the_human_policy() -> None:
    """GC band, motifs, homopolymers and poly-lysine all have to hold."""
    for payload in load_payloads():
        if not payload.is_single_polypeptide:
            continue
        construct, design = _build(payload)
        assert design.is_clean, f"{payload.name}: {design.block_violations}"
        assert min(design.block_gc) >= HUMAN_POLICY.gc_min
        assert max(design.block_gc) <= HUMAN_POLICY.gc_max
        assert design.cai >= HUMAN_POLICY.cai_floor
        assert len(construct.upstream_starts()) == 0


def test_the_constructs_have_the_expected_shape() -> None:
    """Cap, UTRs, signal peptide, payload, UTR and tail, in that order."""
    for payload in load_payloads():
        if not payload.is_single_polypeptide:
            continue
        construct, _design = _build(payload)
        layout = {element.name: element for element in construct.elements}
        assert layout["cap"].start == 1
        assert layout["cap"].end == 2
        assert layout["5'-UTR"].start == 3
        assert layout["5'-UTR"].end == 54
        assert layout["sig"].start == 55
        assert layout["sig"].end == 108
        assert layout["poly(A)"].length == 110
        assert layout["poly(A)"].end == construct.length
        orf = construct.open_reading_frame()
        assert orf.start == 54
        assert orf.end == layout["3'-UTR"].start - 1
        assert construct.length == 54 + orf.end - orf.start + 295 + 110


def test_the_m3_construct_is_reproduced_byte_for_byte() -> None:
    """M5 must not silently change the construct M3 already published."""
    from mrna.fasta import read_fasta

    m3 = read_fasta(M3_FASTA_PATH)
    m5 = read_fasta(FASTA_PATH)
    m3_sequence = m3[0][1]
    m5_sequences = {header.split(" |")[0]: sequence for header, sequence in m5}
    e113 = m5_sequences[
        "E113 against ecarpholin S, the Ser49 myotoxic phospholipase A2"
    ]
    assert e113 == m3_sequence


def test_the_motif_scan_is_scoped_to_the_construct() -> None:
    """The polyadenylation signal is expected in the UTRs, not the coding region."""
    for payload in load_payloads():
        if not payload.is_single_polypeptide:
            continue
        construct, _design = _build(payload)
        hits = construct.find_motifs(DEFAULT_MOTIFS)
        assert hits
        for hit in hits:
            assert hit.element in {"5'-UTR", "3'-UTR"}


# ------------------------------------------------------------------ artifacts


def test_the_committed_fasta_holds_the_constructs_and_their_proteins() -> None:
    """The committed design file must hold what the script would write.

    The byte-for-byte reproducibility of the file is enforced by CI, which
    regenerates it and diffs the tree. What is checked here is the content:
    every transcribable payload contributes exactly one construct record and
    one protein record, and the construct carries the payload.
    """
    from mrna.fasta import read_fasta

    module = _load_script()
    payloads = {payload.name: payload for payload in load_payloads()}
    buildable = plan_cocktail(CandidateTracker(), "single-polypeptide")
    records = read_fasta(FASTA_PATH)
    headers = [header.split(" |")[0] for header, _sequence in records]

    for name in buildable.payload_ids:
        payload = payloads[name]
        construct, _design = module.build_construct(payload, _usage_table())
        construct_headers = [
            header for header in headers if header.startswith(f"{name} against")
        ]
        assert len(construct_headers) == 1
        protein_headers = [
            header for header in headers if header == f"{name} payload protein"
        ]
        assert len(protein_headers) == 1
        sequences = dict(
            zip(headers, (s for _h, s in records), strict=True),
        )
        assert construct.sequence in sequences.values()
        assert payload.sequence in sequences.values()


def test_the_committed_fasta_holds_nothing_else() -> None:
    """Only the transcribable payloads are designed, and only they."""
    from mrna.fasta import read_fasta

    records = read_fasta(FASTA_PATH)
    buildable = plan_cocktail(CandidateTracker(), "single-polypeptide")
    assert len(records) == 2 * len(buildable.payload_ids)


def test_the_committed_elements_file_matches_the_constructs() -> None:
    """Every row of the elements table must be a real element of a construct."""
    rows = ELEMENTS_PATH.read_text(encoding="utf-8").strip().split("\n")
    header = rows[0].split("\t")
    assert header == ["member", "element", "start", "end", "length", "description"]
    payloads = {payload.name: payload for payload in load_payloads()}
    buildable = plan_cocktail(CandidateTracker(), "single-polypeptide")
    by_member = {name: payloads[name] for name in buildable.payload_ids}
    assert len(rows) - 1 == len(by_member) * 6
    for row in rows[1:]:
        member, _element, start, end, length, _description = row.split("\t")
        assert member in by_member
        assert int(end) - int(start) + 1 == int(length)


def test_the_committed_report_is_what_the_script_would_write() -> None:
    """The published report must be reproducible."""
    module = _load_script()
    generated = module.build_report(CandidateTracker())
    assert REPORT_PATH.read_text(encoding="utf-8") == generated


# ---------------------------------------------------- the report tells the truth


def test_the_report_states_the_uncovered_count_the_data_gives(
    plans: tuple[object, ...],
) -> None:
    """The number of uncovered cells in the report must be the real one.

    This is the test that a hardcoded count would fail. The report used to
    say "four cells" when the planner found five.
    """
    report = REPORT_PATH.read_text(encoding="utf-8")
    uncovered = len(plans[0].uncovered)
    assert f"{uncovered} cells have no candidate at all" in report
    for cell in plans[0].uncovered:
        assert cell.describe() in report


def test_the_report_does_not_contradict_the_patent_record(
    plans: tuple[object, ...],
) -> None:
    """The report must not claim a patent covers a payload it does not cover."""
    report = REPORT_PATH.read_text(encoding="utf-8")
    buildable = plans[2]
    blocked = [
        group.antibody_id
        for group in buildable.payloads
        if group.patent_status != "none-located"
    ]
    if not blocked:
        assert "covering 0 of them" in report
    else:
        assert f"covering {len(blocked)} of them" in report


def test_the_report_does_not_overstate_the_big_four_coverage(
    plans: tuple[object, ...],
) -> None:
    """The report must say which Big Four species the plan reaches."""
    report = REPORT_PATH.read_text(encoding="utf-8")
    buildable = plans[2]
    covered = [
        name for name in buildable.species_covered() if name in BIG_FOUR
    ]
    assert f"It covers {len(covered)} Big Four species out of {len(BIG_FOUR)}" in report


def test_the_report_carries_the_disclaimer() -> None:
    """Every artifact carries the disclaimer.

    The elements table is a tabular view of the same constructs and follows
    the M4 precedent of carrying no prose, so the sequence file and the
    report are the two artifacts that carry it.
    """
    disclaimer = "In-silico designs. No therapeutic efficacy demonstrated."
    report = REPORT_PATH.read_text(encoding="utf-8")
    assert disclaimer in report
    assert "Not medical advice." in report
    fasta = FASTA_PATH.read_text(encoding="utf-8")
    assert disclaimer in fasta
    assert "Not medical advice." in fasta
