"""Assembly of a full mRNA construct, and the element table that describes it.

What a construct is
-------------------
An mRNA drug is one linear chain, read 5' to 3'::

    5'-cap | 5' UTR | signal peptide | coding region | 3' UTR | poly(A)- 3'

The cap is a chemical modification, not a stretch of nucleotides. The signal
peptide is translated but is not the payload. Everything else is a sequence
that was chosen. The job of this module is to join those choices in order,
assign every one of them a coordinate, and hand back a table.

Why the element table is the deliverable
-----------------------------------------
A regulatory document describes an mRNA as a table: a name for each element
and a coordinate range for each name. The WHO INN document for Comirnaty
publishes six such rows, and those rows are the only public description of
the geometry of a licensed mRNA drug. A construct that this engine assembles
has to be describable in the same format, because that is how a reader checks
it: line by line against a product a regulator has already approved.

So the output of :func:`assemble_construct` is not a string. It is a sequence
plus the table that explains it, and the table is the part worth reading.

Three things the table has to get right
---------------------------------------
*Coordinates.* Positions are 1-based and inclusive, and they are computed from
the parts rather than asserted, so a table cannot drift away from the
sequence it describes.

*The cap footprint.* The WHO document prints two nucleotides for the cap: the
7-methylguanosine at position 1 and the first backbone nucleotide at position
2, joined by a 5'-to-5' triphosphate. The document then numbers its 5'-UTR
element from position 3. The cap is therefore a chemical modification that
occupies the first two printed positions, and reproducing the document's
table means reproducing that convention exactly. See :data:`CAP1`.

*Stop codons.* A regulatory document counts every printed base, so its
coding-region row includes the stop codons. A coding sequence, by the usual
definition, stops at the first of them. Comirnaty prints two consecutive
stop codons, so its published coding row is three nucleotides longer than
the coding sequence. The assembler takes the coding element as it is printed
and reports the open reading frame separately, so the two never get confused.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from mrna.anatomy import LayoutElement, Orf, find_orf
from mrna.motifs import Motif, find_motifs
from mrna.params import KOZAK_MINUS3_ALLOWED, KOZAK_PLUS4
from mrna.sequence import translate
from mrna.utr import Utr

__all__ = [
    "A110_POLYA",
    "CAP1",
    "Cap",
    "Construct",
    "ConstructPart",
    "ElementMotifHit",
    "PolyATail",
    "UpstreamStart",
    "assemble_construct",
    "full_construct",
]


@dataclass(frozen=True)
class Cap:
    """The 5' cap of a transcript, and the nucleotides a document prints for it.

    A cap is not a sequence. It is a 7-methylguanosine joined to the first
    backbone nucleotide by a 5'-to-5' triphosphate, and a regulatory
    document prints two nucleotides for it so that its coordinates line up
    with the chemistry. The two printed nucleotides are the m7G at position 1
    and the first backbone nucleotide at position 2.

    Attributes:
        name: Name used for the element in the element table.
        prefix: The nucleotides printed for the cap.
        description: What the cap is, in the words of the source.
        citation: Where that description comes from.
    """

    name: str
    prefix: str
    description: str
    citation: str

    @property
    def length(self) -> int:
        """Length of the printed cap footprint, in nucleotides."""
        return len(self.prefix)


#: The cap of a Comirnaty-style transcript: a modified 5'-cap1 structure.
#:
#: The cap is 7-methylguanosine with a 2'-O-methyl on the first backbone
#: nucleotide, and the whole 5' end is written with N1-methylpseudouridine
#: elsewhere in the transcript.
CAP1: Cap = Cap(
    name="cap",
    prefix="GA",
    description=(
        "A modified 5'-cap1 structure (m7G+m3'-5'-ppp-5'-Am). The two "
        "nucleotides printed for the cap are the 7-methylguanosine at "
        "position 1 and the first backbone nucleotide at position 2, joined "
        "by a 5'-to-5' triphosphate."
    ),
    citation=(
        "WHO International Nonproprietary Names Programme document 11889 "
        "(Comirnaty), February 2021, section '5'-capping structure': 'cap "
        "G1A2 = m7G+m3'-5'-ppp-5'-Am2'-3'-p-'. Cap chemistry: Kariko K, "
        "Buckstein M, Ni H, Weissman D. 2005. Suppression of RNA recognition "
        "by Toll-like receptors: the impact of nucleoside modification and "
        "the evolutionary origin of RNA. Immunity 23(2):165-175. "
        "doi:10.1016/j.immuni.2005.06.008. Kariko K, Muramatsu H, Welsh FA, "
        "Ludwig J, Kato H, Akira S, Weissman D. 2008. Incorporation of "
        "pseudouridine into mRNA yields superior nonimmunogenic vector with "
        "increased translational capacity and biological stability. Mol Ther "
        "16(11):1833-1840. doi:10.1038/mt.2008.200."
    ),
)


@dataclass(frozen=True)
class PolyATail:
    """A poly(A) tail, and the reason it has the shape it has.

    The poly(A) tail is what an mRNA is held by, and it is also what a
    ribosome runs into if it reads past the end of the coding region. Its
    length is therefore a compromise between stability and the risk of
    stalling on a poly-lysine tract, which is why the tail is followed by a
    short linker rather than by more adenosines.

    Attributes:
        name: Name used for the element in the element table.
        sequence: The template-encoded tail.
        description: What the tail is, in the words of the source.
        citation: Where that description comes from.
    """

    name: str
    sequence: str
    description: str
    citation: str

    @property
    def length(self) -> int:
        """Length of the tail in nucleotides."""
        return len(self.sequence)


#: The 110-nucleotide template-encoded poly(A) tail of BNT162b2.
#:
#: Thirty adenosines, a ten-nucleotide linker, seventy more adenosines. The
#: linker is there so that a ribosome which misreads past the stop codon
#: runs into a different sequence rather than into a longer poly-lysine
#: tract.
A110_POLYA: PolyATail = PolyATail(
    name="poly(A)",
    sequence="A" * 30 + "GCATATGACT" + "A" * 70,
    description=(
        "A 110-nucleotide poly(A)-tail consisting of a stretch of 30 "
        "adenosine residues, followed by a 10-nucleotide linker sequence and "
        "another 70 adenosine residues."
    ),
    citation=(
        "WHO International Nonproprietary Names Programme document 11889 "
        "(Comirnaty), February 2021, feature table, poly(A) at positions "
        "4175-4284. For poly(A) tail lengths in general: Pardi N, Hogan MJ, "
        "Porter FW, Weissman D. 2018. mRNA vaccines, a new era in "
        "vaccinology. Nat Rev Drug Discov 17:261-279. "
        "doi:10.1038/nrd.2017.243."
    ),
)


@dataclass(frozen=True)
class ConstructPart:
    """One element of a construct.

    Attributes:
        name: Name used for the element in the element table.
        sequence: The nucleotides of the element.
        description: What the element is, in the words of the source.
        citation: Where that description comes from.
    """

    name: str
    sequence: str
    description: str = ""
    citation: str = ""


@dataclass(frozen=True)
class ElementMotifHit:
    """One forbidden motif, found inside one element of a construct.

    Attributes:
        element: Name of the element the motif sits in.
        motif: The motif that matched.
        start: Zero-based index of the first matching base in the construct.
        end: Zero-based index one past the last matching base in the construct.
    """

    element: str
    motif: Motif
    start: int
    end: int

    def describe(self) -> str:
        """Return a human-readable description of the occurrence."""
        return f"{self.motif.name} in {self.element} at {self.start + 1} to {self.end}"


@dataclass(frozen=True)
class UpstreamStart:
    """An ATG inside the 5' UTR, with its initiation context.

    An upstream ATG is a start codon the ribosome could use instead of the
    intended one. It matters most when it sits in a strong initiation
    context, because then the ribosome is likely to start there and produce
    a short, useless protein.

    Attributes:
        index: Zero-based index of the A of the ATG, in the construct.
        context: The nucleotides from position -3 to position +4.
        minus3: The base at position -3, or an empty string if the ATG is
            too close to the 5' end for the position to exist.
        plus4: The base at position +4, or an empty string if the ATG is
            too close to the 3' end for the position to exist.
    """

    index: int
    context: str
    minus3: str
    plus4: str

    @property
    def matches_consensus(self) -> bool:
        """Whether both key Kozak positions follow the consensus.

        Kozak 1986 found that the purine at position -3 and the G at
        position +4 are the two positions that matter most for initiation
        efficiency.
        """
        return self.minus3 in KOZAK_MINUS3_ALLOWED and self.plus4 == KOZAK_PLUS4

    @property
    def is_strong(self) -> bool:
        """Whether the context is strong enough to worry about.

        A context that matches the consensus at both key positions is one a
        ribosome will use in preference to a downstream start. A context that
        matches at neither is one a ribosome will usually ignore.
        """
        return self.matches_consensus


@dataclass(frozen=True)
class Construct:
    """A full mRNA construct and the table that describes it.

    Attributes:
        name: Name of the construct.
        parts: The elements, in 5'-to-3' order.
        elements: The published-style element table, with coordinates.
    """

    name: str
    parts: tuple[ConstructPart, ...]
    elements: tuple[LayoutElement, ...]

    @property
    def sequence(self) -> str:
        """The full construct sequence, 5' to 3'."""
        return "".join(part.sequence for part in self.parts)

    @property
    def length(self) -> int:
        """Total length of the construct in nucleotides."""
        return sum(len(part.sequence) for part in self.parts)

    def part(self, name: str) -> ConstructPart:
        """Return one element by name.

        Args:
            name: The name of the element.

        Returns:
            The element.

        Raises:
            KeyError: If the construct has no element with that name.
        """
        for part in self.parts:
            if part.name == name:
                return part
        msg = f"construct {self.name!r} has no element named {name!r}"
        raise KeyError(msg)

    def open_reading_frame(self) -> Orf:
        """Locate the open reading frame of the assembled construct.

        Returns:
            The located :class:`~mrna.anatomy.Orf`, from the start codon to
            the first in-frame stop codon.

        Raises:
            ValueError: If the construct has no complete open reading frame.
        """
        return find_orf(self.sequence)

    @property
    def protein(self) -> str:
        """The protein the construct encodes, without the stop marker."""
        return translate(self.open_reading_frame().cds).rstrip("*")

    def upstream_starts(self) -> tuple[UpstreamStart, ...]:
        """Find every ATG in the 5' UTR, with its initiation context.

        A start codon in the 5' UTR competes with the intended one. The
        ribosome scans from the 5' end and usually takes the first start
        codon it meets, so an upstream ATG in a strong context can divert
        translation and produce a short protein instead of the antigen.

        Returns:
            The upstream start codons, ordered by position. Empty when the
            construct has no 5' UTR.
        """
        if not self.parts or self.parts[0].name != "cap":
            return ()

        utr5 = self.parts[1] if len(self.parts) > 1 else None
        if utr5 is None or utr5.name != "5'-UTR":
            return ()

        offset = len(self.parts[0].sequence)
        starts: list[UpstreamStart] = []
        index = utr5.sequence.find("ATG")
        while index >= 0:
            left = utr5.sequence[max(0, index - 3) : index]
            right = utr5.sequence[index + 3 : index + 6]
            minus3 = left[0] if len(left) == 3 else ""
            plus4 = right[0] if right else ""
            context = left + "ATG" + right
            starts.append(
                UpstreamStart(
                    index=offset + index,
                    context=context,
                    minus3=minus3,
                    plus4=plus4,
                ),
            )
            index = utr5.sequence.find("ATG", index + 1)
        return tuple(starts)

    def find_motifs(self, motifs: Sequence[Motif]) -> tuple[ElementMotifHit, ...]:
        """Find forbidden motifs, reported per element.

        The scope matters. A polyadenylation signal inside a coding region
        can cause premature cleavage and a transcript with no in-frame stop
        codon, so it is a defect there. The same signal inside a 5' UTR has
        no cleavage site upstream of it and is harmless. Reporting the
        elements separately is what keeps the difference visible.

        Args:
            motifs: The motifs to look for.

        Returns:
            The occurrences, ordered by position in the construct.
        """
        hits: list[ElementMotifHit] = []
        for part, element in zip(self.parts, self.elements, strict=True):
            for hit in find_motifs(part.sequence, tuple(motifs)):
                hits.append(
                    ElementMotifHit(
                        element=part.name,
                        motif=hit.motif,
                        start=element.start - 1 + hit.start,
                        end=element.start - 1 + hit.end,
                    ),
                )
        return tuple(hits)

    def element_table(self) -> str:
        """Return the element table in the format a regulatory document uses.

        Returns:
            A markdown table with the columns ``Element``, ``Description``
            and ``Position``.
        """
        lines = ["| Element | Description | Position |", "| --- | --- | --- |"]
        for element in self.elements:
            lines.append(
                f"| {element.name} | {element.description} | "
                f"{element.start}-{element.end} |",
            )
        return "\n".join(lines)

    def fasta(self, *, width: int = 60) -> str:
        """Return the construct as a FASTA record.

        Args:
            width: Number of sequence characters per line.

        Returns:
            The FASTA record, ending in a newline.
        """
        lines = [f">{self.name}"]
        lines.extend(self.sequence[index : index + width] for index in range(0, self.length, width))
        return "\n".join(lines) + "\n"


def assemble_construct(name: str, parts: Sequence[ConstructPart]) -> Construct:
    """Join elements into a construct and compute its element table.

    The coordinates are computed from the parts, never supplied. A table
    that is written by hand can disagree with the sequence it describes; a
    table that is computed cannot.

    Args:
        name: Name of the construct.
        parts: The elements, in 5'-to-3' order.

    Returns:
        The assembled :class:`Construct`.

    Raises:
        ValueError: If there are no parts, if any part is empty, if a part
            contains a character outside the DNA alphabet, or if two parts
            share a name.
    """
    if not parts:
        msg = "assemble_construct() needs at least one part"
        raise ValueError(msg)

    seen: set[str] = set()
    elements: list[LayoutElement] = []
    offset = 0
    for part in parts:
        if part.name in seen:
            msg = f"assemble_construct() got two parts named {part.name!r}"
            raise ValueError(msg)
        seen.add(part.name)

        sequence = "".join(part.sequence.split()).upper()
        if not sequence:
            msg = f"construct part {part.name!r} is empty"
            raise ValueError(msg)
        invalid = sorted(set(sequence) - set("ACGT"))
        if invalid:
            msg = (
                f"construct part {part.name!r} contains characters outside "
                f"ACGT: {''.join(invalid)}"
            )
            raise ValueError(msg)

        length = len(sequence)
        elements.append(
            LayoutElement(
                name=part.name,
                start=offset + 1,
                end=offset + length,
                description=part.description,
            ),
        )
        offset += length

    normalised = tuple(
        ConstructPart(
            name=part.name,
            sequence="".join(part.sequence.split()).upper(),
            description=part.description,
            citation=part.citation,
        )
        for part in parts
    )
    return Construct(name=name, parts=normalised, elements=tuple(elements))


def full_construct(
    name: str,
    *,
    utr5: Utr,
    coding: str,
    utr3: Utr,
    poly_a: PolyATail,
    cap: Cap = CAP1,
    signal_peptide_nt: int = 0,
    coding_name: str = "antigen",
    coding_description: str = "",
) -> Construct:
    """Assemble a construct from the parts a design actually chooses.

    The coding element is taken as it will be printed, which means it may
    carry stop codons beyond the first one. That is what a regulatory
    document counts, and reproducing a licensed product means reproducing
    its table exactly.

    Args:
        name: Name of the construct.
        utr5: The 5' UTR element.
        coding: The coding element: signal peptide, then antigen, then the
            stop codons as they will be printed.
        utr3: The 3' UTR element.
        poly_a: The poly(A) tail.
        cap: The 5' cap.
        signal_peptide_nt: Length of the signal peptide at the start of the
            coding element, in nucleotides. Zero when the antigen is not
            secreted.
        coding_name: Name for the coding element in the element table.
        coding_description: Description for the coding element.

    Returns:
        The assembled :class:`Construct`.

    Raises:
        ValueError: If the signal peptide is longer than the coding element,
            or if the parts are not assemblable.
    """
    cleaned = "".join(coding.split()).upper()
    if signal_peptide_nt < 0:
        msg = f"full_construct() got a negative signal peptide length {signal_peptide_nt}"
        raise ValueError(msg)
    if signal_peptide_nt > len(cleaned):
        msg = (
            f"full_construct() got a signal peptide of {signal_peptide_nt} nt "
            f"for a coding element of only {len(cleaned)} nt"
        )
        raise ValueError(msg)

    parts = [
        ConstructPart(
            name=cap.name,
            sequence=cap.prefix,
            description=cap.description,
            citation=cap.citation,
        ),
        ConstructPart(
            name="5'-UTR",
            sequence=utr5.sequence,
            description=utr5.description,
            citation=utr5.citation,
        ),
    ]

    if signal_peptide_nt:
        parts.append(
            ConstructPart(
                name="sig",
                sequence=cleaned[:signal_peptide_nt],
                description=(
                    "Signal peptide, which guides translocation of the "
                    "nascent polypeptide chain into the endoplasmic reticulum."
                ),
                citation="",
            ),
        )

    parts.append(
        ConstructPart(
            name=coding_name,
            sequence=cleaned[signal_peptide_nt:],
            description=coding_description,
            citation="",
        ),
    )
    parts.append(
        ConstructPart(
            name="3'-UTR",
            sequence=utr3.sequence,
            description=utr3.description,
            citation=utr3.citation,
        ),
    )
    parts.append(
        ConstructPart(
            name=poly_a.name,
            sequence=poly_a.sequence,
            description=poly_a.description,
            citation=poly_a.citation,
        ),
    )

    return assemble_construct(name, parts)
