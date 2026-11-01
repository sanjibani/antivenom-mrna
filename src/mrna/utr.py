"""A library of untranslated regions, every entry carrying its source.

Biological background
---------------------
Two stretches of an mRNA are never translated. The 5' UTR is the stretch the
ribosome scans through on its way to the start codon, so it is engineered to
be easy to scan and hard to misread. The 3' UTR is never reached by the
ribosome at all, and it is the part that decides how long the message
survives in the cytoplasm and therefore how much protein is made.

Both are borrowed from nature rather than invented. The rule in this project
is that a borrowed element arrives with the citation that justifies it, and
that the same rule applies to the elements this library holds: they are
published, regulator-accepted sequences, not plausible-sounding guesses.

The two entries here are the ones Comirnaty uses, which makes them the
safest possible starting point for an antivenom construct. They are not
claimed to be optimal. They are claimed to be the elements a regulator has
already accepted inside a licensed product, which is a different claim and a
smaller one, and the only one this library can support.

A note on lengths
-----------------
The same 5' UTR is quoted as 52, 53 or 54 nucleotides by different sources,
and all three describe the same molecule. The WHO INN document numbers the
7-methylguanosine cap as position 1 and the first backbone nucleotide as
position 2, and starts its 5'-UTR element at position 3, which makes the
element 52 nt. Counting the first backbone nucleotide as part of the UTR
makes it 53 nt, and counting the cap guanosine as well makes it 54 nt. This
library uses the document's numbering, because the document's coordinates are
what a design has to match.
"""

from __future__ import annotations

from dataclasses import dataclass

from mrna.sequence import gc_content

__all__ = [
    "AES_ELEMENT",
    "AES_MTRNR1_3UTR",
    "ALPHA_GLOBIN_5UTR",
    "MTRNR1_ELEMENT",
    "Utr",
]


@dataclass(frozen=True)
class Utr:
    """One untranslated region, with the source that justifies it.

    Attributes:
        name: Short name used in reports.
        kind: Where the region sits, for example ``5' UTR`` or ``3' UTR``.
        sequence: The region itself, in the DNA alphabet.
        description: What the region is, in the words of the source.
        citation: Where the sequence and the description come from.
        note: Anything a reader would otherwise get wrong.
    """

    name: str
    kind: str
    sequence: str
    description: str
    citation: str
    note: str = ""

    @property
    def length(self) -> int:
        """Length of the region in nucleotides."""
        return len(self.sequence)

    def gc_content(self) -> float:
        """Return the GC fraction of the region.

        Returns:
            The fraction of bases that are G or C.

        Raises:
            ValueError: If the region is empty.
        """
        return gc_content(self.sequence)

    def contains(self, motif_sequence: str) -> bool:
        """Report whether a literal sequence appears inside the region.

        Args:
            motif_sequence: The sequence to look for, in the DNA alphabet.

        Returns:
            True when the sequence occurs in the region.
        """
        needle = "".join(motif_sequence.split()).upper()
        return needle in self.sequence


#: The 136-nucleotide AES element, from the human AES/TLE5 gene.
#:
#: AES stands for amino-terminal enhancer of split, and AES is also known as
#: TLE5. The element was found by a cell-based screen of human UTR fragments
#: and was one of the two that made the Pfizer/BioNTech 3' UTR.
#:
#: Source: Orlandini von Niessen AG, Poleganov MA, Rechner C, Plaschke A,
#: Kranz LM, Fesser M, Diken M, Loewer M, Vallazza B, Beissert T, Bukur V,
#: Kuhn AN, Tureci O, Sahin U. 2019. Improving mRNA-based therapeutic gene
#: delivery by expression-augmenting 3' UTRs identified by cellular library
#: screening. Mol Ther 27(4):824-836. doi:10.1016/j.ymthe.2018.12.011.
AES_ELEMENT: Utr = Utr(
    name="AES element",
    kind="3' UTR element",
    sequence=(
        "CTGGTACTGCATGCACGCAATGCTAGCTGCCCCTTTCCCGTCCTGGGTACCCCGAGTCTCCCCCGACCTCGGGT"
        "CCCAGGTATGCTCCCACCTCCACCTGCCCCACTCACCACCTCTGCTAGTTCCAGACACCTCC"
    ),
    description=(
        "The AES element, from the amino-terminal enhancer of split mRNA, used "
        "as the 5' half of the 3' UTR of BNT162b2."
    ),
    citation=(
        "Orlandini von Niessen AG, Poleganov MA, Rechner C, Plaschke A, Kranz "
        "LM, Fesser M, Diken M, Loewer M, Vallazza B, Beissert T, Bukur V, Kuhn "
        "AN, Tureci O, Sahin U. 2019. Improving mRNA-based therapeutic gene "
        "delivery by expression-augmenting 3' UTRs identified by cellular "
        "library screening. Mol Ther 27(4):824-836. "
        "doi:10.1016/j.ymthe.2018.12.011. The element is also printed in WHO "
        "International Nonproprietary Names Programme document 11889 "
        "(Comirnaty), February 2021, as part of the 3'-UTR at positions "
        "3880-4174."
    ),
    note=(
        "AES is a member of the Groucho/Transducin-like enhancer of split gene "
        "family. It regulates androgen receptor transcriptional activity and "
        "Notch and Wnt signalling, and it acts as a tumour suppressor. The "
        "stability function this project relies on was discovered by the "
        "screen above, not assumed from its role in transcription."
    ),
)


#: The 142-nucleotide mitochondrial 12S rRNA element.
#:
#: mtRNR1 is the mitochondrial-encoded 12S ribosomal RNA, which is involved in
#: translating the 13 protein-coding genes of mitochondrial DNA. It was the
#: second of the two elements that made the Pfizer/BioNTech 3' UTR.
#:
#: Source: Orlandini von Niessen AG et al. 2019, as above.
MTRNR1_ELEMENT: Utr = Utr(
    name="mitochondrial 12S rRNA element",
    kind="3' UTR element",
    sequence=(
        "CAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCCACACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCA"
        "ATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACACC"
    ),
    description=(
        "The mitochondrial 12S ribosomal RNA element, used as the 3' half of "
        "the 3' UTR of BNT162b2."
    ),
    citation=(
        "Orlandini von Niessen AG, Poleganov MA, Rechner C, Plaschke A, Kranz "
        "LM, Fesser M, Diken M, Loewer M, Vallazza B, Beissert T, Bukur V, Kuhn "
        "AN, Tureci O, Sahin U. 2019. Improving mRNA-based therapeutic gene "
        "delivery by expression-augmenting 3' UTRs identified by cellular "
        "library screening. Mol Ther 27(4):824-836. "
        "doi:10.1016/j.ymthe.2018.12.011. The element is also printed in WHO "
        "International Nonproprietary Names Programme document 11889 "
        "(Comirnaty), February 2021, as part of the 3'-UTR at positions "
        "3880-4174."
    ),
    note=(
        "Mutations in mtRNR1 are associated with deafness, which is a "
        "different matter from the stability function measured by the screen "
        "that found this element."
    ),
)


#: The 295-nucleotide 3' UTR of BNT162b2: AES, then the 12S rRNA element.
#:
#: This is the published 3' UTR, taken from the primary regulatory document
#: and decomposed into the two elements the document names.
#:
#: Source: WHO International Nonproprietary Names Programme document 11889
#: (Comirnaty), February 2021, feature table, 3'-UTR at positions 3880-4174.
AES_MTRNR1_3UTR: Utr = Utr(
    name="AES plus 12S rRNA 3'UTR",
    kind="3' UTR",
    sequence=(
        "CTCGAG"
        "CTGGTACTGCATGCACGCAATGCTAGCTGCCCCTTTCCCGTCCTGGGTACCCCGAGTCTCCCCCGACCTCGGGT"
        "CCCAGGTATGCTCCCACCTCCACCTGCCCCACTCACCACCTCTGCTAGTTCCAGACACCTCC"
        "CAAGCACGCAGCAATGCAGCTCAAAACGCTTAGCCTAGCCACACCCCCACGGGAAACAGCAGTGATTAACCTTTAGCA"
        "ATAAACGAAAGTTTAACTAAGCTATACTAACCCCAGGGTTGGTCAATTTCGTGCCAGCCACACC"
        "CTGGAGCTAGC"
    ),
    description=(
        "A 3' untranslated region comprising two sequence elements derived "
        "from the amino-terminal enhancer of split (AES) mRNA and the "
        "mitochondrial encoded 12S ribosomal RNA, to confer RNA stability and "
        "high total protein expression. As used in BNT162b2."
    ),
    citation=(
        "WHO International Nonproprietary Names Programme document 11889 "
        "(Comirnaty), February 2021, feature table, 3'-UTR at positions "
        "3880-4174. The two constituent elements were found by Orlandini von "
        "Niessen AG et al. 2019, Mol Ther 27(4):824-836, "
        "doi:10.1016/j.ymthe.2018.12.011."
    ),
    note=(
        "The published 295 nt decomposes into a 6-nt prefix, the 136-nt AES "
        "element, the 142-nt mtRNR1 element and an 11-nt suffix. The prefix "
        "CTCGAG is the recognition sequence of the restriction endonuclease "
        "XhoI; the source does not say whether it is a deliberate cloning scar "
        "or part of the element. The order in the source's own sentence, 'AES "
        "mRNA and the mitochondrial encoded 12S ribosomal RNA', is the order "
        "in the sequence."
    ),
)


#: The 52-nucleotide 5' UTR of BNT162b2, from human alpha-globin RNA.
#:
#: Source: WHO International Nonproprietary Names Programme document 11889
#: (Comirnaty), February 2021, feature table, 5'-UTR at positions 3-54.
ALPHA_GLOBIN_5UTR: Utr = Utr(
    name="alpha-globin 5'UTR",
    kind="5' UTR",
    sequence="GAATAAACTAGTATTCTTCTGGTCCCCACAGACTCAGAGAGAACCCGCCACC",
    description=(
        "5'-untranslated region derived from human alpha-globin RNA with an "
        "optimized Kozak sequence. As used in BNT162b2."
    ),
    citation=(
        "WHO International Nonproprietary Names Programme document 11889 "
        "(Comirnaty), February 2021, feature table, 5'-UTR at positions 3-54. "
        "The same 52 nucleotides are printed independently in Fang E, Liu X, "
        "Li M, Zhang Z, Song L, Zhu B, Wu X, Liu J. 2022. Advances in COVID-19 "
        "mRNA vaccine development. Signal Transduct Target Ther 7(1):94, Table "
        "1, doi:10.1038/s41392-022-00950-y."
    ),
    note=(
        "The last six nucleotides, GCCACC, are the Kozak sequence. Kozak M. "
        "1986. Point mutation defines a sequence flanking the AUG initiator "
        "codon that modulates translation by eukaryotic ribosomes. Cell "
        "44(2):283-292. The element ends immediately before the start codon, "
        "so this UTR and the coding region form a single continuous "
        "initiation context rather than two separate parts."
    ),
)
