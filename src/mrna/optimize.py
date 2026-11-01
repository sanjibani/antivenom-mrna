"""The codon optimiser.

The optimiser turns a protein into a coding sequence that is good for human
ribosomes, holds its GC content inside a band, and contains no motif that
would sabotage its manufacture.

How it works
------------
The objective is the Codon Adaptation Index, which is the geometric mean of
the relative adaptiveness of the codons (Sharp and Li, 1987). Maximising the
CAI is the same as maximising the sum of ``log(w)`` over the codons, so the
optimiser is a dynamic program over codons where the state is the GC count
inside the current block and the value is that sum.

Because the GC constraint applies to fixed blocks and blocks do not interact
with each other, each block is solved independently. That makes the result
optimal rather than merely good: *the CAI is maximal subject to the GC
constraint*, and that sentence can go into a design document.

Three repair passes run afterwards, each with a biological reason and each
able to undo itself if it would break the GC band:

* poly-lysine runs of consecutive AAA codons, which stall ribosomes,
* forbidden motifs such as the AAUAAA polyadenylation signal,
* long homopolymer runs, which destabilise templates and make T7 polymerase
  slip.

The optimiser is deterministic. The same inputs produce byte-identical
output, because a design that changes when you re-run it cannot be audited.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from mrna.anatomy import codon_adaptation_index, relative_adaptiveness
from mrna.motifs import MotifHit, find_motifs, homopolymer_runs
from mrna.policy import DesignPolicy
from mrna.sequence import AA_TO_CODONS, translate

__all__ = [
    "OptimizationResult",
    "RepairResult",
    "codon_agreement",
    "naive_back_translation",
    "optimize_cds",
    "repair_sequence",
]


@dataclass(frozen=True)
class OptimizationResult:
    """Everything the optimiser measured about one designed coding sequence.

    Attributes:
        protein: The protein the sequence encodes.
        cds: The designed coding sequence.
        cai: Codon Adaptation Index of the designed sequence.
        gc_total: GC fraction of the whole coding sequence.
        block_gc: GC fraction of every block.
        block_violations: Human-readable descriptions of blocks that could
            not be brought inside the band, if any.
        window_gc_min: Lowest GC fraction in any 50-nucleotide window.
        window_gc_max: Highest GC fraction in any 50-nucleotide window.
        motif_hits: Forbidden motifs still present after repair.
        homopolymer_violations: Homopolymer runs still present after repair.
        polylysine_violations: Poly-lysine runs still present after repair.
        cai_is_optimal: True when every block was solved inside the GC band.
    """

    protein: str
    cds: str
    cai: float
    gc_total: float
    block_gc: tuple[float, ...]
    block_violations: tuple[str, ...]
    window_gc_min: float
    window_gc_max: float
    motif_hits: tuple[MotifHit, ...]
    homopolymer_violations: tuple[tuple[int, int], ...]
    polylysine_violations: tuple[tuple[int, int], ...]
    cai_is_optimal: bool

    @property
    def is_clean(self) -> bool:
        """Whether the design satisfies every constraint."""
        return (
            not self.block_violations
            and not self.motif_hits
            and not self.homopolymer_violations
            and not self.polylysine_violations
        )


def _gc_count(codon: str) -> int:
    """Return the number of G or C bases in a codon."""
    return sum(1 for base in codon if base in "GC")


@dataclass(frozen=True)
class RepairResult:
    """The outcome of running the repair passes over a coding sequence.

    Attributes:
        cds: The repaired coding sequence.
        motif_hits: Forbidden motifs still present.
        homopolymer_violations: Homopolymer runs still present.
        polylysine_violations: Poly-lysine runs still present.
    """

    cds: str
    motif_hits: tuple[MotifHit, ...]
    homopolymer_violations: tuple[tuple[int, int], ...]
    polylysine_violations: tuple[tuple[int, int], ...]

    @property
    def is_clean(self) -> bool:
        """Whether every repair succeeded."""
        return (
            not self.motif_hits
            and not self.homopolymer_violations
            and not self.polylysine_violations
        )


def _candidates(
    protein: str,
    usage_table: Mapping[str, float],
) -> dict[str, tuple[tuple[str, float], ...]]:
    """Build the candidate codons for each amino acid in a protein.

    Args:
        protein: The protein to encode.
        usage_table: Codon usage table for the target host.

    Returns:
        Mapping of amino acid to a tuple of ``(codon, relative adaptiveness)``
        pairs, ordered by descending adaptiveness then alphabetically.

    Raises:
        ValueError: If the usage table does not cover every needed codon.
    """
    weights = relative_adaptiveness(usage_table)
    needed = {codon for amino_acid in set(protein) for codon in AA_TO_CODONS[amino_acid]}
    missing = sorted(needed - set(weights))
    if missing:
        msg = f"the usage table is missing {len(missing)} codons: {', '.join(missing)}"
        raise ValueError(msg)

    candidates: dict[str, tuple[tuple[str, float], ...]] = {}
    for amino_acid in sorted(set(protein)):
        options = [
            (codon, weights[codon])
            for codon in AA_TO_CODONS[amino_acid]
            if weights[codon] > 0
        ]
        options.sort(key=lambda item: (-item[1], item[0]))
        candidates[amino_acid] = tuple(options)
    return candidates


def _window(low: int, high: int, length: int) -> tuple[int, int]:
    """Return the allowed GC count range for a block of a given length."""
    return max(0, low), min(length, high)


def _optimize_block(
    residues: Sequence[str],
    candidates: Mapping[str, Sequence[tuple[str, float]]],
    low: int,
    high: int,
) -> tuple[tuple[str, ...], int, bool]:
    """Solve one GC block optimally.

    Args:
        residues: The amino acids in the block.
        candidates: Candidate codons for each amino acid.
        low: Minimum allowed GC count for the block.
        high: Maximum allowed GC count for the block.

    Returns:
        A tuple of the chosen codons, the GC count they achieve, and whether
        the count was inside the allowed range.
    """
    states: dict[int, tuple[float, tuple[str, ...]]] = {0: (0.0, ())}
    for amino_acid in residues:
        following: dict[int, tuple[float, tuple[str, ...]]] = {}
        for gc, (score, path) in states.items():
            for codon, weight in candidates[amino_acid]:
                new_gc = gc + _gc_count(codon)
                value = score + math.log(weight)
                current = following.get(new_gc)
                if current is None or value > current[0]:
                    following[new_gc] = (value, (*path, codon))
        states = following

    feasible = {gc: value for gc, value in states.items() if low <= gc <= high}
    if feasible:
        centre = (low + high) / 2
        best_gc, best = max(
            feasible.items(),
            key=lambda item: (item[1][0], -abs(item[0] - centre)),
        )
        return best[1], best_gc, True

    def distance(gc: int) -> int:
        return max(low - gc, gc - high, 0)

    best_gc, best = min(
        states.items(),
        key=lambda item: (distance(item[0]), -item[1][0]),
    )
    return best[1], best_gc, False


def _block_ranges(length: int, block_nt: int) -> tuple[tuple[int, int], ...]:
    """Return the nucleotide ranges of every GC block."""
    return tuple(
        (start, min(start + block_nt, length)) for start in range(0, length, block_nt)
    )


def _block_window(
    policy: DesignPolicy,
    length: int,
) -> tuple[int, int]:
    """Return the allowed GC count range for a block of a given length."""
    if policy.gc_target is not None:
        centre = round(policy.gc_target * length)
        tolerance = max(1, round(policy.gc_block_tolerance_nt * length / policy.gc_block_nt))
        return _window(centre - tolerance, centre + tolerance, length)
    low = math.ceil(policy.gc_min * length)
    high = math.floor(policy.gc_max * length)
    return _window(low, high, length)


def _gc_in_range(cds: str, start: int, end: int) -> int:
    """Return the GC count of a nucleotide range."""
    return sum(1 for base in cds[start:end] if base in "GC")


def _swap(
    cds: str,
    codon_index: int,
    new_codon: str,
    policy: DesignPolicy,
) -> str:
    """Return the coding sequence with one codon replaced."""
    start = codon_index * 3
    return cds[:start] + new_codon + cds[start + 3 :]


def _block_of(codon_index: int, block_nt: int) -> int:
    """Return the index of the block a codon belongs to."""
    return (codon_index * 3) // block_nt


def _repair_motifs(
    cds: str,
    policy: DesignPolicy,
    candidates: Mapping[str, Sequence[tuple[str, float]]],
) -> str:
    """Remove forbidden motifs by changing codons where possible."""
    for _ in range(200):
        hits = find_motifs(cds, policy.motifs)
        if not hits:
            break
        hit = hits[0]
        first_codon = hit.start // 3
        last_codon = (hit.end - 1) // 3
        repaired = False
        for codon_index in range(first_codon, last_codon + 1):
            amino_acid = translate(cds[codon_index * 3 : codon_index * 3 + 3])
            for codon, _weight in candidates.get(amino_acid, ()):
                if codon == cds[codon_index * 3 : codon_index * 3 + 3]:
                    continue
                trial = _swap(cds, codon_index, codon, policy)
                if find_motifs(trial, policy.motifs):
                    continue
                if not _gc_not_worse(cds, trial, codon_index, policy):
                    continue
                cds = trial
                repaired = True
                break
            if repaired:
                break
        if not repaired:
            break
    return cds


def _repair_homopolymers(
    cds: str,
    policy: DesignPolicy,
    candidates: Mapping[str, Sequence[tuple[str, float]]],
) -> str:
    """Shorten homopolymer runs by changing codons where possible."""
    for _ in range(200):
        runs = homopolymer_runs(cds, policy.max_homopolymer)
        if not runs:
            break
        start, end = runs[0]
        first_codon = start // 3
        last_codon = (end - 1) // 3
        repaired = False
        for codon_index in range(first_codon, last_codon + 1):
            current = cds[codon_index * 3 : codon_index * 3 + 3]
            amino_acid = translate(current)
            for codon, _weight in candidates.get(amino_acid, ()):
                if codon == current:
                    continue
                trial = _swap(cds, codon_index, codon, policy)
                if len(homopolymer_runs(trial, policy.max_homopolymer)) >= len(runs):
                    continue
                if not _gc_not_worse(cds, trial, codon_index, policy):
                    continue
                cds = trial
                repaired = True
                break
            if repaired:
                break
        if not repaired:
            break
    return cds


def _polylysine_runs(cds: str, limit: int) -> tuple[tuple[int, int], ...]:
    """Find runs of consecutive AAA lysine codons longer than the limit."""
    codons = [cds[index : index + 3] for index in range(0, len(cds) - 2, 3)]
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for index, codon in enumerate([*codons, ""]):
        if codon == "AAA":
            if start is None:
                start = index
        elif start is not None:
            if index - start > limit:
                runs.append((start, index))
            start = None
    return tuple(runs)


def _repair_polylysine(
    cds: str,
    policy: DesignPolicy,
    candidates: Mapping[str, Sequence[tuple[str, float]]],
) -> str:
    """Break poly-lysine runs by changing an AAA codon to another lysine codon."""
    for _ in range(200):
        runs = _polylysine_runs(cds, policy.max_polylysine_codons)
        if not runs:
            break
        start, end = runs[0]
        repaired = False
        for codon_index in range(start, end):
            current = cds[codon_index * 3 : codon_index * 3 + 3]
            if current != "AAA":
                continue
            for codon, _weight in candidates.get("K", ()):
                if codon == current:
                    continue
                trial = _swap(cds, codon_index, codon, policy)
                if len(_polylysine_runs(trial, policy.max_polylysine_codons)) >= len(runs):
                    continue
                if not _gc_not_worse(cds, trial, codon_index, policy):
                    continue
                cds = trial
                repaired = True
                break
            if repaired:
                break
        if not repaired:
            break
    return cds


def _gc_not_worse(
    before: str,
    after: str,
    codon_index: int,
    policy: DesignPolicy,
) -> bool:
    """Check that a swap does not degrade the GC profile of its block.

    A swap must keep its block inside the GC band. If the block was already
    outside the band because the protein admits no other encoding, then the
    swap is allowed as long as it does not push the block further out. That
    matters: a forbidden motif has to be repairable even in a protein whose
    GC content cannot be fixed, and a motif left in place is a worse defect
    than a band that was already out of reach.
    """
    block_nt = policy.gc_block_nt
    block_index = _block_of(codon_index, block_nt)
    ranges = _block_ranges(len(before), block_nt)
    if block_index >= len(ranges):
        return True
    start, end = ranges[block_index]
    low, high = _block_window(policy, end - start)
    gc_before = _gc_in_range(before, start, end)
    gc_after = _gc_in_range(after, start, end)
    if low <= gc_before <= high:
        return low <= gc_after <= high
    return abs(gc_after - (low + high) / 2) <= abs(gc_before - (low + high) / 2)


def optimize_cds(
    protein: str,
    usage_table: Mapping[str, float],
    policy: DesignPolicy | None = None,
) -> OptimizationResult:
    """Design a coding sequence for a protein.

    Args:
        protein: The protein to encode, one-letter codes.
        usage_table: Codon usage table for the target host.
        policy: The design policy. Defaults to the human policy.

    Returns:
        The designed sequence and everything measured about it.

    Raises:
        ValueError: If the protein contains characters that are not amino
            acids.
    """
    active = policy or DesignPolicy()
    cleaned = "".join(protein.split()).upper()
    unknown = sorted(set(cleaned) - set(AA_TO_CODONS))
    if unknown:
        msg = f"optimize_cds() found non-standard amino acids {''.join(unknown)}"
        raise ValueError(msg)
    if not cleaned:
        msg = "optimize_cds() received an empty protein"
        raise ValueError(msg)

    candidates = _candidates(cleaned, usage_table)

    codons: list[str] = []
    block_gc: list[float] = []
    violations: list[str] = []
    optimal = True
    for start, end in _block_ranges(len(cleaned) * 3, active.gc_block_nt):
        residues = cleaned[start // 3 : end // 3]
        low, high = _block_window(active, end - start)
        block_codons, gc_count, feasible = _optimize_block(
            residues, candidates, low, high,
        )
        codons.extend(block_codons)
        block_gc.append(gc_count / (end - start))
        if not feasible:
            optimal = False
            band = f"{active.gc_min:.2f} to {active.gc_max:.2f}"
            violations.append(
                f"block {len(block_gc)} (nt {start + 1} to {end}) could only reach "
                f"GC {gc_count / (end - start):.3f}, band is {band}",
            )

    cds = "".join(codons)
    cds = _repair_polylysine(cds, active, candidates)
    cds = _repair_motifs(cds, active, candidates)
    cds = _repair_homopolymers(cds, active, candidates)

    measured_blocks: list[float] = []
    measured_violations: list[str] = []
    for start, end in _block_ranges(len(cds), active.gc_block_nt):
        length = end - start
        gc = _gc_in_range(cds, start, end)
        measured_blocks.append(gc / length)
        low, high = _block_window(active, length)
        if not low <= gc <= high:
            optimal = False
            measured_violations.append(
                f"block {len(measured_blocks)} (nt {start + 1} to {end}) sits at GC "
                f"{gc / length:.3f}, band is {active.gc_min:.2f} to {active.gc_max:.2f}",
            )

    windows = [
        gc / 50
        for gc in (
            sum(1 for base in cds[index : index + 50] if base in "GC")
            for index in range(0, max(1, len(cds) - 49))
        )
    ]

    return OptimizationResult(
        protein=cleaned,
        cds=cds,
        cai=codon_adaptation_index(cds, usage_table),
        gc_total=sum(1 for base in cds if base in "GC") / len(cds),
        block_gc=tuple(measured_blocks),
        block_violations=tuple(violations + measured_violations),
        window_gc_min=min(windows) if windows else 0.0,
        window_gc_max=max(windows) if windows else 0.0,
        motif_hits=find_motifs(cds, active.motifs),
        homopolymer_violations=homopolymer_runs(cds, active.max_homopolymer),
        polylysine_violations=_polylysine_runs(cds, active.max_polylysine_codons),
        cai_is_optimal=optimal,
    )


def repair_sequence(
    cds: str,
    usage_table: Mapping[str, float],
    policy: DesignPolicy | None = None,
) -> RepairResult:
    """Repair a coding sequence without changing the protein it encodes.

    This is the repair stage of the optimiser, exposed on its own so that a
    sequence that came from somewhere else, for example a collaborator's
    design, can be put through the same checks.

    Args:
        cds: The coding sequence to repair.
        usage_table: Codon usage table for the target host.
        policy: The design policy. Defaults to the human policy.

    Returns:
        The repaired sequence and the problems that remain.

    Raises:
        ValueError: If the sequence length is not a multiple of three.
    """
    active = policy or DesignPolicy()
    cleaned = "".join(cds.split()).upper()
    if len(cleaned) % 3:
        msg = f"repair_sequence() needs a length divisible by 3, got {len(cleaned)}"
        raise ValueError(msg)

    protein = translate(cleaned).rstrip("*")
    candidates = _candidates(protein, usage_table)
    repaired = _repair_polylysine(cleaned, active, candidates)
    repaired = _repair_motifs(repaired, active, candidates)
    repaired = _repair_homopolymers(repaired, active, candidates)

    return RepairResult(
        cds=repaired,
        motif_hits=find_motifs(repaired, active.motifs),
        homopolymer_violations=homopolymer_runs(repaired, active.max_homopolymer),
        polylysine_violations=_polylysine_runs(repaired, active.max_polylysine_codons),
    )


def naive_back_translation(
    protein: str,
    usage_table: Mapping[str, float],
) -> str:
    """Encode a protein with the single most used codon for each amino acid.

    This is what a designer gets with no optimiser: the most frequent codon
    of each amino acid, with no attention to GC content, to motifs, or to
    anything else. It is the baseline the optimiser has to beat.

    Args:
        protein: The protein to encode.
        usage_table: Codon usage table for the target host.

    Returns:
        The naive coding sequence.
    """
    candidates = _candidates(protein, usage_table)
    return "".join(candidates[amino_acid][0][0] for amino_acid in protein.upper())


def codon_agreement(first: str, second: str) -> float:
    """Measure how many codons two coding sequences choose identically.

    Args:
        first: The first coding sequence.
        second: The second coding sequence.

    Returns:
        The fraction of positions where the two sequences use the same codon,
        between 0 and 1.

    Raises:
        ValueError: If the two sequences have different lengths.
    """
    left = "".join(first.split()).upper()
    right = "".join(second.split()).upper()
    if len(left) != len(right):
        msg = f"codon_agreement() needs equal lengths, got {len(left)} and {len(right)}"
        raise ValueError(msg)
    if not left:
        msg = "codon_agreement() received an empty sequence"
        raise ValueError(msg)
    matches = sum(
        1
        for index in range(0, len(left), 3)
        if left[index : index + 3] == right[index : index + 3]
    )
    return matches / (len(left) // 3)
