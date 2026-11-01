"""Public API of the antivenom-mrna design engine.

The engine is deliberately small at this stage. M0 shipped a primitive sequence
layer (:mod:`mrna.sequence`) and a measurement layer (:mod:`mrna.anatomy`).
M1 added the codon optimiser (:mod:`mrna.optimize`) and the cited design
policy it obeys (:mod:`mrna.policy`). M2 added the UTR library
(:mod:`mrna.utr`) and the full-construct assembler (:mod:`mrna.assemble`).
M3 added the signal peptide library (:mod:`mrna.leader`). M4 added the
candidate table (:mod:`mrna.candidates`) and the tracker that validates it
(:mod:`mrna.tracker`). M5 added the antibody payloads a construct can carry
(:mod:`mrna.payload`) and the cocktail planner that decides which of them a
Big Four cocktail needs (:mod:`mrna.cocktail`).

The design policy constants live in :mod:`mrna.params`, and every one of them
carries the citation that justifies it.
"""

from __future__ import annotations

from mrna.anatomy import (
    ConstructReport,
    KozakContext,
    LayoutCheck,
    LayoutElement,
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
from mrna.assemble import (
    A110_POLYA,
    CAP1,
    Cap,
    Construct,
    ConstructPart,
    ElementMotifHit,
    PolyATail,
    UpstreamStart,
    assemble_construct,
    full_construct,
)
from mrna.candidates import (
    BIG_FOUR,
    CANDIDATES,
    CLINICAL_SET,
    EVIDENCE_LEVELS,
    MRNA_SUITABILITY,
    SPECIES_FAMILIES,
    SPECIES_TOXIN_FAMILIES,
    TOXIN_FAMILIES,
    Candidate,
)
from mrna.cocktail import (
    EXCLUDED_FAMILIES,
    TIERS,
    Cell,
    CocktailPlan,
    PayloadGroup,
    Tier,
    excluded_cells,
    plan_all_tiers,
    plan_cocktail,
    target_cells,
)
from mrna.fasta import read_fasta, write_fasta
from mrna.leader import HSA_SIGNAL_PEPTIDE, Leader
from mrna.motifs import CANDIDATE_MOTIFS, DEFAULT_MOTIFS, Motif, MotifHit
from mrna.optimize import (
    OptimizationResult,
    RepairResult,
    codon_agreement,
    naive_back_translation,
    optimize_cds,
    repair_sequence,
)
from mrna.payload import PAYLOAD_SPECS, Payload, PayloadSpec, load_payload, load_payloads
from mrna.policy import HUMAN_POLICY, DesignPolicy
from mrna.reference import load_codon_usage, load_single_sequence, reference_dir
from mrna.sequence import (
    AA_TO_CODONS,
    CODON_TABLE,
    STOP_CODONS,
    back_translate,
    gc_content,
    translate,
)
from mrna.tracker import CandidateTracker, Gap, derive_mrna_suitability
from mrna.utr import (
    AES_ELEMENT,
    AES_MTRNR1_3UTR,
    ALPHA_GLOBIN_5UTR,
    MTRNR1_ELEMENT,
    Utr,
)

__version__ = "0.5.0"

__all__ = [
    "A110_POLYA",
    "AA_TO_CODONS",
    "AES_ELEMENT",
    "AES_MTRNR1_3UTR",
    "ALPHA_GLOBIN_5UTR",
    "BIG_FOUR",
    "CANDIDATES",
    "CANDIDATE_MOTIFS",
    "CAP1",
    "CLINICAL_SET",
    "CODON_TABLE",
    "DEFAULT_MOTIFS",
    "EVIDENCE_LEVELS",
    "EXCLUDED_FAMILIES",
    "HSA_SIGNAL_PEPTIDE",
    "HUMAN_POLICY",
    "MRNA_SUITABILITY",
    "MTRNR1_ELEMENT",
    "PAYLOAD_SPECS",
    "SPECIES_FAMILIES",
    "SPECIES_TOXIN_FAMILIES",
    "STOP_CODONS",
    "TIERS",
    "TOXIN_FAMILIES",
    "Candidate",
    "CandidateTracker",
    "Cap",
    "Cell",
    "CocktailPlan",
    "Construct",
    "ConstructPart",
    "ConstructReport",
    "DesignPolicy",
    "ElementMotifHit",
    "Gap",
    "KozakContext",
    "LayoutCheck",
    "LayoutElement",
    "Leader",
    "Motif",
    "MotifHit",
    "OptimizationResult",
    "Orf",
    "Payload",
    "PayloadGroup",
    "PayloadSpec",
    "PolyATail",
    "RepairResult",
    "Tier",
    "UpstreamStart",
    "Utr",
    "Variant",
    "analyze_cds",
    "assemble_construct",
    "back_translate",
    "codon_adaptation_index",
    "codon_agreement",
    "codon_counts",
    "derive_mrna_suitability",
    "excluded_cells",
    "find_orf",
    "full_construct",
    "gc_content",
    "kozak_context",
    "load_codon_usage",
    "load_payload",
    "load_payloads",
    "load_single_sequence",
    "naive_back_translation",
    "optimize_cds",
    "plan_all_tiers",
    "plan_cocktail",
    "protein_variants",
    "read_fasta",
    "reference_dir",
    "relative_adaptiveness",
    "repair_sequence",
    "split_construct",
    "target_cells",
    "translate",
    "write_fasta",
]
