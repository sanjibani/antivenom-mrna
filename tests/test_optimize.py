"""Tests for the codon optimiser.

The optimiser has one job that matters: emit a sequence that is correct,
constrained and reproducible. These tests check each of those, and they check
that the optimiser reports a violation rather than hiding one.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

from mrna.anatomy import codon_adaptation_index, find_orf
from mrna.fasta import read_fasta
from mrna.optimize import (
    codon_agreement,
    naive_back_translation,
    optimize_cds,
    repair_sequence,
)
from mrna.policy import HUMAN_POLICY, DesignPolicy
from mrna.sequence import translate

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"
RAW = REFERENCE / "raw"


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


def spike_antigen() -> str:
    """Return the Comirnaty antigen protein."""
    assemblies = read_fasta(RAW / "naalytics_figure1figure2_032321.fasta")
    return translate(find_orf(assemblies[0][1]).cds).rstrip("*")


def comirnaty_cds() -> str:
    """Return the published Comirnaty coding sequence without its stop codon."""
    assemblies = read_fasta(RAW / "naalytics_figure1figure2_032321.fasta")
    return find_orf(assemblies[0][1]).cds[:-3]


def test_output_translates_back_to_the_input_protein() -> None:
    """Correctness comes first."""
    usage = load_usage_table()
    protein = spike_antigen()
    result = optimize_cds(protein, usage, HUMAN_POLICY)
    assert translate(result.cds) == protein


def test_output_is_deterministic() -> None:
    """A design that changes when you re-run it cannot be audited."""
    usage = load_usage_table()
    protein = spike_antigen()
    first = optimize_cds(protein, usage, HUMAN_POLICY)
    second = optimize_cds(protein, usage, HUMAN_POLICY)
    assert first.cds == second.cds


def test_every_gc_block_is_inside_the_band() -> None:
    """The band is the constraint the optimiser exists to satisfy."""
    usage = load_usage_table()
    result = optimize_cds(spike_antigen(), usage, HUMAN_POLICY)
    for index, value in enumerate(result.block_gc):
        assert HUMAN_POLICY.gc_min - 1e-9 <= value <= HUMAN_POLICY.gc_max + 1e-9, index
    assert result.block_violations == ()


def test_the_design_has_no_forbidden_motif() -> None:
    """No motif may survive in a clean design."""
    usage = load_usage_table()
    result = optimize_cds(spike_antigen(), usage, HUMAN_POLICY)
    assert result.motif_hits == ()
    assert result.is_clean


def test_the_design_reaches_the_cai_target() -> None:
    """The CAI target is a target, not a hope."""
    usage = load_usage_table()
    result = optimize_cds(spike_antigen(), usage, HUMAN_POLICY)
    assert result.cai >= HUMAN_POLICY.cai_target


def test_the_optimiser_agrees_with_comirnaty_on_most_codons() -> None:
    """The benchmark: given the protein alone, we rediscover their choices.

    This is the regression baseline for the headline number of M1. If a change
    to the optimiser breaks the agreement, the change made the tool worse,
    not better.
    """
    usage = load_usage_table()
    result = optimize_cds(spike_antigen(), usage, HUMAN_POLICY)
    agreement = codon_agreement(result.cds, comirnaty_cds())
    assert agreement == pytest.approx(0.782, abs=0.01)


def test_the_viral_gene_is_a_much_worse_coding_sequence() -> None:
    """The reason codon optimisation exists, as a test.

    The virus writes this protein for human cells at a CAI of 0.67. If the
    optimiser stopped beating that, something would be wrong.
    """
    usage = load_usage_table()
    viral = read_fasta(RAW / "berthubert_ncov-s.fasta")[0][1][:-3]
    result = optimize_cds(spike_antigen(), usage, HUMAN_POLICY)
    assert codon_adaptation_index(viral, usage) < 0.75
    assert result.cai > 0.95


def test_the_repair_removes_a_forbidden_motif() -> None:
    """The repair pass must actually repair, and must not change the protein.

    Asparagine followed by lysine encodes as AAT AAA, which is the canonical
    polyadenylation signal sitting in the middle of a coding region.
    """
    usage = load_usage_table()
    cds = "ATG" + "AATAAA" + "GCT"
    assert "AATAAA" in cds
    result = repair_sequence(cds, usage, HUMAN_POLICY)
    assert "AATAAA" not in result.cds
    assert result.motif_hits == ()
    assert translate(result.cds) == translate(cds)


def test_a_polylysine_run_is_broken() -> None:
    """Consecutive AAA codons stall ribosomes, so they must not survive."""
    usage = load_usage_table()
    protein = "MKKKKKKKKKM"
    result = optimize_cds(protein, usage, HUMAN_POLICY)
    assert result.polylysine_violations == ()
    assert "AAAAAAAAAA" not in result.cds


def test_an_infeasible_block_is_reported_not_hidden() -> None:
    """A stretch with no codon choice must produce a stated violation."""
    usage = load_usage_table()
    policy = DesignPolicy(gc_min=0.60, gc_max=0.65, gc_block_nt=30)
    result = optimize_cds("M" * 40, usage, policy)
    assert result.block_violations
    assert result.cai_is_optimal is False
    assert not result.is_clean


def test_the_optimality_claim_holds_when_the_band_is_wide() -> None:
    """With a feasible band the design must be reported as optimal."""
    usage = load_usage_table()
    result = optimize_cds(spike_antigen(), usage, HUMAN_POLICY)
    assert result.cai_is_optimal is True


def test_a_tight_target_reports_infeasibility() -> None:
    """Tightening the band beyond what the protein allows must be reported."""
    usage = load_usage_table()
    policy = DesignPolicy(gc_target=0.5699, gc_block_tolerance_nt=1)
    result = optimize_cds(spike_antigen(), usage, policy)
    assert result.block_violations or result.cai_is_optimal is False


def test_naive_back_translation_uses_the_most_used_codon() -> None:
    """The baseline must be what it claims to be."""
    usage = load_usage_table()
    protein = "MKVLA"
    naive = naive_back_translation(protein, usage)
    assert translate(naive) == protein
    assert codon_adaptation_index(naive, usage) == pytest.approx(1.0)


def test_codon_agreement_measures_the_right_thing() -> None:
    """The benchmark helper must be correct."""
    assert codon_agreement("ATGGCTTAA", "ATGGCTTAA") == 1.0
    assert codon_agreement("ATGGCTTAA", "ATGGCTGCT") == pytest.approx(2 / 3)
    with pytest.raises(ValueError, match="equal lengths"):
        codon_agreement("ATGGCTTAA", "ATGGCT")


def test_optimize_rejects_bad_input() -> None:
    """Empty and invalid proteins must be refused."""
    usage = load_usage_table()
    with pytest.raises(ValueError, match="empty protein"):
        optimize_cds("", usage, HUMAN_POLICY)
    with pytest.raises(ValueError, match="non-standard"):
        optimize_cds("MKV!X", usage, HUMAN_POLICY)


def test_policy_block_window_requires_codon_aligned_blocks() -> None:
    """A block that does not align to codons is a bug in the policy."""
    with pytest.raises(ValueError, match="multiple of 3"):
        DesignPolicy(gc_block_nt=31).block_window()


@given(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=120))
def test_random_proteins_round_trip(protein: str) -> None:
    """Any protein must come back unchanged."""
    usage = load_usage_table()
    result = optimize_cds(protein, usage, HUMAN_POLICY)
    assert translate(result.cds) == protein


@given(st.text(alphabet="ACDEFGHIKLMNPQRSTVWY", min_size=1, max_size=120))
def test_random_proteins_never_leave_a_forbidden_motif(protein: str) -> None:
    """No protein may produce a design with a forbidden motif."""
    usage = load_usage_table()
    result = optimize_cds(protein, usage, HUMAN_POLICY)
    assert result.motif_hits == ()
