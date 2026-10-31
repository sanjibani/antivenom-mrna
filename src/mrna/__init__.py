"""Public API of the antivenom-mrna design engine.

The engine is deliberately small at this stage. M0 ships a primitive sequence
layer (:mod:`mrna.sequence`) and a measurement layer (:mod:`mrna.anatomy`).
Design policy lives in :mod:`mrna.params`.
"""

from __future__ import annotations

from mrna.anatomy import (
    ConstructReport,
    KozakContext,
    Orf,
    Variant,
    analyze_cds,
    codon_adaptation_index,
    codon_counts,
    find_orf,
    kozak_context,
    protein_variants,
    relative_adaptiveness,
    split_construct,
)
from mrna.fasta import read_fasta, write_fasta
from mrna.sequence import (
    AA_TO_CODONS,
    CODON_TABLE,
    STOP_CODONS,
    back_translate,
    gc_content,
    translate,
)

__version__ = "0.1.0"

__all__ = [
    "AA_TO_CODONS",
    "CODON_TABLE",
    "STOP_CODONS",
    "ConstructReport",
    "KozakContext",
    "Orf",
    "Variant",
    "analyze_cds",
    "back_translate",
    "codon_adaptation_index",
    "codon_counts",
    "find_orf",
    "gc_content",
    "kozak_context",
    "protein_variants",
    "read_fasta",
    "relative_adaptiveness",
    "split_construct",
    "translate",
    "write_fasta",
]
