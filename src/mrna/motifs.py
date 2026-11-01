"""Forbidden sequence motifs for mRNA constructs.

A motif is a short sequence that must not appear inside a designed coding
region. Each one carries the reason it is forbidden and a citation, and each
one is either verified or explicitly marked as a candidate. Nothing enters
the default policy without a citation, which is the project rule: a number or
a pattern with no source is a guess wearing the costume of a rule.

The distinction between verified and candidate matters here. Two of the
motifs below are established in the literature and are safe to enforce. The
cloning sites depend on which DNA template a manufacturing partner will use,
which this project does not yet know, so they stay candidates until someone
is told which template it is.
"""

from __future__ import annotations

import re
from dataclasses import dataclass

__all__ = [
    "CANDIDATE_MOTIFS",
    "DEFAULT_MOTIFS",
    "Motif",
    "MotifHit",
    "find_motifs",
    "homopolymer_runs",
]


@dataclass(frozen=True)
class Motif:
    """A forbidden sequence pattern.

    Attributes:
        name: Short name used in reports.
        pattern: Regular expression that matches the motif.
        reason: Why the motif is a problem.
        citation: Where the reason is established.
        verified: True when the citation supports enforcing the motif.
    """

    name: str
    pattern: str
    reason: str
    citation: str
    verified: bool = True

    def compiled(self) -> re.Pattern[str]:
        """Return the compiled regular expression for this motif."""
        return re.compile(self.pattern)


@dataclass(frozen=True)
class MotifHit:
    """One occurrence of a forbidden motif in a sequence.

    Attributes:
        motif: The motif that matched.
        start: Zero-based index of the first matching base.
        end: Zero-based index one past the last matching base.
    """

    motif: Motif
    start: int
    end: int

    def describe(self) -> str:
        """Return a human-readable description of the occurrence."""
        return f"{self.motif.name} at {self.start + 1} to {self.end}"


#: Motifs that are enforced by default.
DEFAULT_MOTIFS: tuple[Motif, ...] = (
    Motif(
        name="polyadenylation signal",
        pattern="AATAAA",
        reason=(
            "AAUAAA is the canonical polyadenylation signal of vertebrate "
            "pre-mRNAs, sitting 10 to 30 nucleotides upstream of the cleavage "
            "site. A copy inside a coding region can be used as a cryptic "
            "polyadenylation site, which yields a truncated transcript with no "
            "in-frame stop codon, and a ribosome that then reads into the "
            "poly(A) tail and stalls on a poly-lysine tract."
        ),
        citation=(
            "Proudfoot NJ, Brownlee BG. 1976. AATAAA is a signal for synthesis "
            "of polyadenylated RNA. Nucleic Acids Res 3(7):1183-1190. "
            "Fitzgerald PC, Shenk T. 1981. The sequence AAUAAA is required for "
            "accurate and efficient polyadenylation of the 5' long control "
            "region of SV40. Proc Natl Acad Sci USA 78(7):4182-4185. "
            "For the observation that polyadenylation inside a coding region "
            "produces a transcript without an in-frame stop codon, see "
            "'Signals for pre-mRNA cleavage and polyadenylation', "
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4451228/"
        ),
    ),
)

#: Motifs that are documented but not enforced, because the constraint they
#: represent depends on information this project does not have yet.
CANDIDATE_MOTIFS: tuple[Motif, ...] = (
    Motif(
        name="NotI site",
        pattern="GCGGCCGC",
        reason=(
            "A cloning site used in many linear DNA templates. A copy inside "
            "the coding region has to be removed by changing a codon, not the "
            "protein."
        ),
        citation=(
            "Template dependent. The restriction sites actually present in the "
            "chosen IVT template must be confirmed with the manufacturing "
            "partner before this motif can be enforced."
        ),
        verified=False,
    ),
    Motif(
        name="ApaI site",
        pattern="GGGCCC",
        reason="A cloning site used in some linear DNA templates.",
        citation=(
            "Template dependent, see the NotI entry above."
        ),
        verified=False,
    ),
    Motif(
        name="SapI site",
        pattern="GCTCTTC",
        reason=(
            "A type IIS cloning site used in some linear DNA templates. It is "
            "also the target of the guide-directed nuclease used in some "
            "template workflows."
        ),
        citation=(
            "Template dependent, see the NotI entry above."
        ),
        verified=False,
    ),
)


def find_motifs(sequence: str, motifs: tuple[Motif, ...]) -> tuple[MotifHit, ...]:
    """Find every occurrence of the given motifs in a sequence.

    Args:
        sequence: The nucleotide sequence to scan.
        motifs: The motifs to look for.

    Returns:
        The occurrences, ordered by position then by motif name.
    """
    cleaned = "".join(sequence.split()).upper()
    hits: list[MotifHit] = []
    for motif in motifs:
        for match in motif.compiled().finditer(cleaned):
            hits.append(MotifHit(motif=motif, start=match.start(), end=match.end()))
    hits.sort(key=lambda hit: (hit.start, hit.motif.name))
    return tuple(hits)


def homopolymer_runs(sequence: str, max_run: int) -> tuple[tuple[int, int], ...]:
    """Find runs of one identical base longer than a limit.

    Long homopolymer tracts are a manufacturing liability: they destabilise
    plasmid templates and they make T7 RNA polymerase slip, which broadens the
    length distribution of the transcript.

    Args:
        sequence: The nucleotide sequence to scan.
        max_run: Longest tolerated run. Runs of exactly this length are
            allowed, runs one longer are reported.

    Returns:
        The runs as ``(start, end)`` pairs, zero-based and end-exclusive.
    """
    cleaned = "".join(sequence.split()).upper()
    runs: list[tuple[int, int]] = []
    if not cleaned:
        return tuple()
    start = 0
    for index in range(1, len(cleaned) + 1):
        if index == len(cleaned) or cleaned[index] != cleaned[start]:
            if index - start > max_run:
                runs.append((start, index))
            start = index
    return tuple(runs)
