"""The P2A ribosome-skipping element, the sourced way to put two payloads on one transcript.

Biological background
---------------------
A 2A peptide is a short viral sequence that makes the ribosome stop forming
the peptide bond between the last glycine and the following proline of its
own encoding. Translation continues on the messenger RNA, so one open
reading frame yields two proteins: an upstream product that carries the 2A
residues at its C-terminus, ending on that glycine, and a downstream product
that begins with the proline. Nothing is cleaved: the bond was never made.

For this repository the element answers the question M5 left open. A
cocktail is several antibodies, and several antibodies on several
transcripts is several drugs. One transcript carrying two payloads is one
drug, and the 2A element is the published, characterised way to build it.

What the element costs
----------------------
The skip is not free of consequence, and the cost is asymmetric. The
upstream product keeps a C-terminal scar: every residue of the spacer and
the 2A peptide up to the skip glycine. The downstream product begins with
one extra residue, the proline, and when the downstream payload is secreted
that proline is followed by its signal peptide, which signal peptidase
removes together with the proline. The downstream payload therefore comes
out exactly as published, and the upstream payload does not.

That asymmetry is a design rule: the antibody whose published data matters
more belongs downstream, where the transcript reproduces it residue for
residue.

Provenance
----------
The sequence is the porcine teschovirus-1 2A peptide characterised by Kim
and colleagues. It is recorded here exactly as it appears in two structures
deposited by different groups, in both cases preceded by the same glycine-
serine-glycine spacer the characterisation recommends:

* PDB entry 31FW, entity 4, residues 51 to 69, in a human tRNA ligase
  complex expression construct.
* PDB entry 8E99, entity 2, residues 4 to 22, in an NMDA receptor
  expression construct.

Both raw entries are committed under ``data/reference/raw/`` and the test
suite extracts the peptide from each and checks the two agree with this
constant, which is the two-source rule of M0 applied to a new element.

Sources: Kim JH, Lee SR, Li LH, et al. 2011. High cleavage efficiency of a
2A peptide derived from porcine teschovirus-1 in human cell lines, zebrafish
and mice. PLoS ONE 6(4):e18556. doi:10.1371/journal.pone.0018556. Donnelly
ML, Hughes LE, Luke G, et al. 2001. Analysis of the aphthovirus 2A/2B
polyprotein 'cleavage' mechanism indicates not a proteolytic reaction, but a
novel translational effect: a polyprotein 'CIS-ACTING ELEMENT' (2A) and
sequence requirements. J Gen Virol 82:1027-1041.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["GSG_P2A_LINKER", "P2A", "P2A_DOWNSTREAM_PROLINE", "TwoA"]


@dataclass(frozen=True)
class TwoA:
    """A ribosome-skipping element, with the sources that justify it.

    The element is a stretch of protein, not of DNA: the ribosome reads the
    peptide as it is being made, and the DNA that encodes it is a codon
    optimiser's business, not the library's.

    Attributes:
        name: Name used for the element in the element table.
        sequence: The peptide, in the one-letter amino acid alphabet, as
            deposited. The final proline is the first residue of the
            downstream product, and is part of the deposited 19-mer because
            that is how the constructs that carry it are sequenced.
        description: What the element is, in the words of the source.
        citation: Where the sequence comes from.
    """

    name: str
    sequence: str
    description: str
    citation: str

    @property
    def upstream_retained(self) -> str:
        """The residues that stay on the C-terminus of the upstream product.

        The skip happens between the final glycine and the proline, so the
        upstream product ends on the glycine: every residue of the deposited
        peptide except the final proline.
        """
        return self.sequence[:-1]

    @property
    def downstream_proline(self) -> str:
        """The single residue that begins the downstream product."""
        return self.sequence[-1]


#: The porcine teschovirus-1 2A peptide, from two deposited constructs.
P2A = TwoA(
    name="P2A",
    sequence="ATNFSLLKQAGDVEENPGP",
    description=(
        "The porcine teschovirus-1 2A ribosome-skipping peptide, preceded by "
        "the glycine-serine-glycine spacer both deposited constructs use"
    ),
    citation=(
        "PDB entry 31FW entity 4 residues 51-69 and PDB entry 8E99 entity 2 "
        "residues 4-22, two expression constructs deposited by different "
        "groups carrying the identical peptide; characterised by Kim JH, Lee "
        "SR, Li LH, et al. 2011. PLoS ONE 6(4):e18556. "
        "doi:10.1371/journal.pone.0018556"
    ),
)

#: The spacer and the peptide together, as the transcript encodes them.
#:
#: The glycine-serine-glycine spacer is the arrangement both source
#: structures use and the characterisation recommends. The upstream product
#: carries all of it: spacer plus the peptide up to the skip glycine, 21
#: residues in all.
GSG_P2A_LINKER: str = "GSG" + P2A.sequence[:-1]

#: The residue the downstream product begins with.
P2A_DOWNSTREAM_PROLINE: str = P2A.downstream_proline
