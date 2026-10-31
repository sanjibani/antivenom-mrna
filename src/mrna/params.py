"""Design parameters for the antivenom-mrna engine, each with a citation.

Policy lives here, not in :mod:`mrna.sequence`. Every number in this module is
a *choice* rather than a measurement, so every number carries the citation that
justifies it. If you add a constant to this file, it needs a source. If you
cannot find a source, the constant does not go in the code yet.

Where a value is a consensus rather than a hard threshold, the docstring says
so, because the difference between "published as a target" and "published as
an observation" matters when a reviewer asks why a construct looks the way it
does.
"""

from __future__ import annotations

__all__ = [
    "CAI_TARGET",
    "GC_MAX",
    "GC_MIN",
    "GC_TARGET",
    "GC_WINDOW_NT",
    "KOZAK_CONSENSUS",
    "KOZAK_MINUS3_ALLOWED",
    "KOZAK_PLUS4",
    "POLYA_MAX",
    "POLYA_MIN",
]

#: Target GC content for a designed coding region.
#:
#: Source: Mauro & Bhatt 2022, median GC content across approved mRNA vaccines,
#: as specified by the project owner on 2026-10-31. This value is a *median
#: observation of licensed products*, not a hard rule, so it is treated as a
#: target rather than a constraint. The citation is pending verification
#: against the paper itself; until then, every construct that uses this value
#: also reports the GC content of BNT162b2 and mRNA-1273 computed by this
#: repository, so the target can be audited against first-party data.
#:
#: Cross-check computed by this repo (see ``docs/m0_vax_anatomy.md``): the two
#: licensed vaccines used as the calibration set.
GC_TARGET: float = 0.55

#: Lower bound of the acceptable GC band for a coding region.
#:
#: Source: Mauro & Bhatt 2022 (as above) for the accepted range around
#: :data:`GC_TARGET`. Below this, translation efficiency drops and the
#: transcript becomes less stable; above it, secondary structure becomes the
#: dominant failure mode.
GC_MIN: float = 0.50

#: Upper bound of the acceptable GC band for a coding region.
#:
#: Source: Mauro & Bhatt 2022 (as above).
GC_MAX: float = 0.60

#: Target Codon Adaptation Index for a designed coding region.
#:
#: Source: Sharp and Li 1987, "The codon adaptation index, a measure of
#: directional synonymous codon usage bias, and its potential applications".
#: Mol Evol Biol 26(3):151-161. The paper reports CAI values above 0.7 for
#: highly expressed genes and below 0.3 for lowly expressed ones in the
#: organisms surveyed. Treat 0.7 as a target for a drug construct whose host
#: cell is a hepatocyte or a muscle cell, not as a therapeutic threshold:
#: mRNA half-life, not ribosome throughput alone, decides exposure.
CAI_TARGET: float = 0.70

#: Length of the sliding window used for local GC inspection.
#:
#: Source: Kudla, Wozniak, Hossmann and Krainer 1999, "High guanine and
#: cytosine content increases mRNA levels in mammalian cells". Nucleic Acids
#: Res 27(9):2085-2088. That paper shows the GC content of the region around
#: the start codon (the first ~50 coding nucleotides) changes translation
#: efficiency in mammalian cells, so 50 nt is the natural window for asking
#: "is the 5' end of this CDS unusual?".
GC_WINDOW_NT: int = 50

#: The Kozak consensus sequence around an initiation codon.
#:
#: Source: Kozak 1986, "Point mutation defines a sequence flanking the AUG
#: initiator codon that modulates translation by eukaryotic ribosomes".
#: Cell 44(2):283-292. The consensus given there is gccRccATGG, with the
#: purine at position -3 and the G at position +4 (numbering the A of AUG as
#: +1) having the largest effect on initiation efficiency.
KOZAK_CONSENSUS: str = "gccRccATGG"

#: Bases allowed at position -3 of a Kozak consensus sequence.
#:
#: Source: Kozak 1986 (as above). A purine here is the single most important
#: position for initiation efficiency.
KOZAK_MINUS3_ALLOWED: frozenset[str] = frozenset("AG")

#: Base required at position +4 of a Kozak consensus sequence.
#:
#: Source: Kozak 1986 (as above).
KOZAK_PLUS4: str = "G"

#: Minimum poly(A) tail length, in nucleotides, for a designed construct.
#:
#: Source: Pardi, Hogan, Porter and Weissman 2018, "mRNA vaccines, a new era
#: in vaccinology". Nat Rev Drug Discov 17:261-279, which describes poly(A)
#: tails of roughly 100 to 150 nt as the standard choice. The BNT162b2
#: construct uses a template-encoded ``A30(GCATATGACT)A70`` segment, that is
#: 100 nt, followed by a longer enzymatically added tail.
POLYA_MIN: int = 100

#: Maximum poly(A) tail length, in nucleotides, for a designed construct.
#:
#: Source: Pardi et al. 2018 (as above).
POLYA_MAX: int = 150
