"""N-terminal leaders: the signal peptides that route a payload to the outside.

Biological background
---------------------
A therapeutic antibody made from mRNA has to leave the cell that made it.
The mRNA is delivered into muscle, the payload is translated in the cytosol,
and the protein has to travel into the endoplasmic reticulum, through the
Golgi and out of the cell. A signal peptide is the address label that sends
it: a short hydrophobic stretch at the N-terminus that is recognised by the
signal recognition particle, docks at the translocon, and is then cut off by
signal peptidase inside the ER lumen.

Two things follow from that, and both matter for a design.

First, the signal peptide is a *choice*, and the right one is the one that
has already worked in a construct someone else got approved. The human serum
albumin signal peptide is one of the two most commonly used in mRNA vaccine
designs, which makes it the safe default for a project that has no wet lab.

Second, signal peptidase removes the leader, so the mature protein starts
exactly at the residue after it. That is why the leader is kept as a
separate element of the construct rather than being folded into the payload:
the coordinate table of a licensed product lists the leader and the antigen
as different elements, and so does this engine.

The leader is not a UTR. It is translated, it is cleaved off, and it never
appears in the mature protein. :mod:`mrna.utr` holds the parts of the message
that are never translated at all.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["HSA_SIGNAL_PEPTIDE", "Leader"]


@dataclass(frozen=True)
class Leader:
    """An N-terminal leader, with the source that justifies it.

    The leader is a stretch of protein, not of DNA. It is described in amino
    acids because that is what it is: a signal peptidase reads the protein
    and cuts it, and the DNA that encodes it is a codon-optimiser's
    business, not the library's.

    Attributes:
        name: Name used for the element in the element table.
        sequence: The leader, in the one-letter amino acid alphabet.
        description: What the leader is, in the words of the source.
        citation: Where that description comes from.
    """

    name: str
    sequence: str
    description: str
    citation: str

    @property
    def length(self) -> int:
        """Length of the leader in amino acids."""
        return len(self.sequence)


#: The 18-residue human serum albumin signal peptide.
#:
#: This is the leader of the human serum albumin precursor, residues 1 to 18.
#: Albumin is made by the liver and is the most abundant protein in human
#: plasma, which is why its leader is a good choice for a secreted payload:
#: the cell already knows how to read it.
#:
#: Note the sequence: ``MKWVTFISLLFLFSSAYS``. The h-region is ``VTFISLLFLF``,
#: with two phenylalanines rather than three leucines. That is what the
#: UniProt record says and what the structural literature quotes, and it is
#: worth writing down rather than typing from memory.
HSA_SIGNAL_PEPTIDE: Leader = Leader(
    name="sig",
    sequence="MKWVTFISLLFLFSSAYS",
    description=(
        "Human serum albumin signal peptide, which directs the nascent "
        "polypeptide into the endoplasmic reticulum for secretion."
    ),
    citation=(
        "Human serum albumin precursor, UniProt P02768, signal peptide, "
        "residues 1 to 18. For its use in mRNA vaccine designs: 'The human "
        "serum albumin (HSA) (MKWVTFISLLFLFSSAYS) and H1 HAWAII SPs are two "
        "commonly used SPs in mRNA vaccine designs', Structure-Based "
        "Optimization of Pathogen Signal Sequences for Enhanced Antigen "
        "Expression in Humans for Vaccine Designs, bioRxiv 2025, "
        "doi:10.1101/2025.11.25.690495 (preprint). For the efficiency of "
        "this leader in mammalian cells: Attallah C, Etcheverrigaray M, "
        "Kratje R, Oggero M. 2017. A highly efficient modified human serum "
        "albumin signal peptide to secrete proteins in cells derived from "
        "different mammalian species. Protein Expr Purif 132:27-33. "
        "doi:10.1016/j.pep.2017.01.003."
    ),
)
