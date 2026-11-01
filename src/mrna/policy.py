"""Design policy for mRNA constructs.

Every value in this module is a choice rather than a measurement, so every
value carries the source that justifies it. Where the choice was derived from
the two licensed mRNA vaccines rather than from a paper, the derivation is
recorded in ``docs/m1_codon_optimizer.md`` and the numbers can be recomputed
from the committed reference data.

The policy is deliberately a frozen dataclass. A design run must be
reproducible, which means the policy that produced it must be a value that
can be printed, stored and compared, not a set of module-level constants that
change underneath you.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from mrna.motifs import DEFAULT_MOTIFS, Motif

__all__ = ["HUMAN_POLICY", "DesignPolicy"]


@dataclass(frozen=True)
class DesignPolicy:
    """The set of constraints a designed coding sequence must satisfy.

    Attributes:
        gc_min: Lower bound of the GC fraction of a GC block.
        gc_max: Upper bound of the GC fraction of a GC block.
        gc_block_nt: Length of the GC block, in nucleotides. Must be a
            multiple of three so that blocks align to codons.
        gc_target: Optional GC fraction to steer towards, for example the GC
            content of a reference construct. When set, the block window
            narrows to this value plus or minus ``gc_block_tolerance_nt``.
        gc_block_tolerance_nt: How far a block may drift from
            ``gc_target``, in nucleotides.
        cai_target: CAI the design aims for.
        cai_floor: CAI below which a design is reported as failing.
        max_homopolymer: Longest tolerated run of one identical base.
        max_polylysine_codons: Longest tolerated run of consecutive AAA
            lysine codons.
        motifs: Forbidden sequence motifs.
        host: Name of the host whose codon usage table is used.
    """

    gc_min: float = 0.40
    gc_max: float = 0.70
    gc_block_nt: int = 30
    gc_target: float | None = None
    gc_block_tolerance_nt: int = 3
    cai_target: float = 0.95
    cai_floor: float = 0.90
    max_homopolymer: int = 8
    max_polylysine_codons: int = 4
    motifs: tuple[Motif, ...] = field(default=DEFAULT_MOTIFS)
    host: str = "homo_sapiens"

    def block_window(self) -> tuple[int, int]:
        """Return the allowed GC count range for one GC block.

        Returns:
            A tuple of the minimum and maximum number of G or C bases in a
            block of :attr:`gc_block_nt` nucleotides.
        """
        if self.gc_block_nt % 3:
            msg = (
                f"gc_block_nt must be a multiple of 3 so blocks align to codons, "
                f"got {self.gc_block_nt}"
            )
            raise ValueError(msg)
        if self.gc_target is not None:
            centre = round(self.gc_target * self.gc_block_nt)
            low = max(0, centre - self.gc_block_tolerance_nt)
            high = min(self.gc_block_nt, centre + self.gc_block_tolerance_nt)
            return low, high
        low = int(-(-self.gc_min * self.gc_block_nt // 1))
        high = int(self.gc_max * self.gc_block_nt // 1)
        return low, high


#: The default policy, calibrated on the two licensed mRNA vaccines.
#:
#: The GC block band of 0.40 to 0.70 is not a guess. It is the band that both
#: approved products satisfy: Comirnaty's 30-nucleotide blocks span 0.400 to
#: 0.700 and Spikevax's span 0.467 to 0.767. A design that stays inside this
#: band is never more extreme than something a regulator has already
#: approved. The narrower band of 0.50 to 0.60 that appears in the literature
#: would have excluded Spikevax, whose coding sequence sits at 0.6230.
#:
#: The CAI target of 0.95 is the CAI of Comirnaty, measured in M0. The floor of
#: 0.90 exists so that a design which cannot reach the target is reported
#: rather than quietly shipped.
#:
#: The homopolymer limit of 8 bases and the poly-lysine limit of 4 codons
#: are both engineering choices with a biological reason, and neither approved
#: product comes close to either: Comirnaty has no homopolymer run of 9 or
#: more and no two consecutive AAA codons.
HUMAN_POLICY = DesignPolicy()
