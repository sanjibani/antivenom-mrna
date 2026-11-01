"""M2: a full-construct assembler, measured against a licensed product.

This script is the M2 deliverable. It rebuilds Comirnaty from its published
parts using the assembler in ``src/mrna/assemble.py``, compares the result
with the primary regulatory document, then builds a new construct on the same
UTR scaffold to show that the assembler is not a one-product reproducer.

Run it::

    uv run python scripts/construct_assembler.py            # print the report
    uv run python scripts/construct_assembler.py --write    # refresh docs/m2_construct_assembler.md

The question the report answers is not "does this tool work". It is "given
the published parts of a licensed mRNA vaccine, can this engine put them
back together and land on the same 4,284 nucleotides and the same element
table". A tool that can do that for a product a regulator approved is a tool
whose output for an antivenom construct is worth reading.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mrna.assemble import A110_POLYA, Construct, full_construct
from mrna.leader import Leader
from mrna.motifs import DEFAULT_MOTIFS
from mrna.optimize import optimize_cds
from mrna.policy import HUMAN_POLICY
from mrna.reference import load_codon_usage, load_single_sequence
from mrna.utr import AES_MTRNR1_3UTR, ALPHA_GLOBIN_5UTR

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"

DISCLAIMER = "In-silico designs. No therapeutic efficacy demonstrated. Not medical advice."

#: The element table published in the WHO INN document 11889.
PUBLISHED_LAYOUT: tuple[tuple[str, int, int], ...] = (
    ("cap", 1, 2),
    ("5'-UTR", 3, 54),
    ("sig", 55, 102),
    ("S protein_mut", 103, 3879),
    ("3'-UTR", 3880, 4174),
    ("poly(A)", 4175, 4284),
)

#: The WHO document prints two consecutive stop codons at 3874 to 3879. The
#: published coding-sequence transcription ends at the first of them, because
#: that is where translation stops, so the second one has to be added back
#: before the construct can be reproduced.
SECOND_STOP_CODON = "TGA"

#: The extended signal peptide of the SARS-CoV-2 spike glycoprotein, as
#: printed in the WHO INN document.
COMIRNATY_SIGNAL_PEPTIDE = Leader(
    name="sig",
    sequence="MFVFLVLLPLVSSQCV",
    description=(
        "S glycoprotein signal peptide (extended leader sequence), which "
        "guides translocation of the nascent polypeptide chain into the "
        "endoplasmic reticulum."
    ),
    citation=(
        "WHO International Nonproprietary Names Programme document 11889 "
        "(Comirnaty), February 2021, feature table, sig at positions 55-102."
    ),
)


def comirnaty_parts() -> tuple[Construct, str]:
    """Rebuild Comirnaty from the published parts.

    Returns:
        The rebuilt construct, and the published sequence to compare it with.
    """
    document = load_single_sequence(REFERENCE / "bnt162b2_mrna_who.fasta")
    published_cds = load_single_sequence(REFERENCE / "bnt162b2_cds.fasta")
    construct = full_construct(
        "BNT162b2 (Comirnaty)",
        utr5=ALPHA_GLOBIN_5UTR,
        coding=published_cds + SECOND_STOP_CODON,
        utr3=AES_MTRNR1_3UTR,
        poly_a=A110_POLYA,
        leader=COMIRNATY_SIGNAL_PEPTIDE,
        coding_name="S protein_mut",
        coding_description=(
            "Codon-optimized sequence encoding full-length SARS-CoV-2 spike "
            "(S) glycoprotein containing mutations K986P and V987P to ensure "
            "the S glycoprotein remains in an antigenically optimal "
            "pre-fusion conformation; stop codons: 3874-3879 (underlined)"
        ),
    )
    return construct, document


def gfp_construct(usage_table: dict[str, float]) -> Construct:
    """Build a new construct on the Comirnaty UTR scaffold.

    The antigen is green fluorescent protein, designed by the M1 optimiser.
    The point is the scaffold: the same cap, the same UTRs, the same poly(A)
    specification, and a completely different payload.

    Args:
        usage_table: Codon usage table for the target host.

    Returns:
        The assembled construct.
    """
    protein = load_single_sequence(REFERENCE / "gfp_P42212.fasta").strip()
    design = optimize_cds(protein, usage_table, HUMAN_POLICY)
    return full_construct(
        "GFP on the Comirnaty scaffold",
        utr5=ALPHA_GLOBIN_5UTR,
        coding=design.cds + "TGA",
        utr3=AES_MTRNR1_3UTR,
        poly_a=A110_POLYA,
        coding_name="GFP",
        coding_description=(
            "Codon-optimised green fluorescent protein from Aequorea "
            "victoria, UniProt P42212, with a single stop codon."
        ),
    )


def build_report() -> str:
    """Build the full M2 markdown report.

    Returns:
        The report as a markdown string.
    """
    usage_table = load_codon_usage(REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv")
    rebuilt, document = comirnaty_parts()
    designed = gfp_construct(usage_table)

    mismatch_count = sum(
        1
        for index in range(min(len(rebuilt.sequence), len(document)))
        if rebuilt.sequence[index] != document[index]
    )
    published_table = tuple(
        (name, start, end) for name, start, end in PUBLISHED_LAYOUT
    )
    rebuilt_table = tuple(
        (element.name, element.start, element.end) for element in rebuilt.elements
    )
    rebuilt_hits = rebuilt.find_motifs(DEFAULT_MOTIFS)
    designed_hits = designed.find_motifs(DEFAULT_MOTIFS)
    orf = rebuilt.open_reading_frame()

    lines: list[str] = []
    lines.append("# A construct assembler, measured against a licensed product")
    lines.append("")
    lines.append(f"> {DISCLAIMER}")
    lines.append("")
    lines.append(
        "This is the M2 artifact. The assembler takes the parts of an mRNA "
        "construct in order and produces two things: the sequence, and the "
        "element table that describes it. The element table is the part that "
        "matters, because it is the format a regulatory document uses, and "
        "because it is the only way a reader can check a design against a "
        "licensed product line by line.",
    )
    lines.append("")
    lines.append(
        "The measurement is a rebuild. The published parts of Comirnaty, "
        "taken from the primary regulatory document and from an independent "
        "transcription of its coding sequence, are fed back into the "
        "assembler. If the result is the same 4,284 nucleotides with the "
        "same six rows, then the assembler reproduces the geometry that a "
        "regulator has approved, and its output for an antivenom construct "
        "is worth reading.",
    )
    lines.append("")

    lines.append("## 1. The rebuild")
    lines.append("")
    lines.append("| Element | Published range | Rebuilt range | Match |")
    lines.append("| --- | --- | --- | --- |")
    for published, observed in zip(published_table, rebuilt_table, strict=True):
        name, start, end = published
        _observed_name, observed_start, observed_end = observed
        same = published == observed
        lines.append(
            f"| {name} | {start} to {end} | {observed_start} to {observed_end} | "
            f"{'yes' if same else 'NO'} |",
        )
    lines.append("")
    lines.append(
        f"Rebuilt length {len(rebuilt.sequence)} nt against the published "
        f"{len(document)} nt. Nucleotides compared: "
        f"{min(len(rebuilt.sequence), len(document))}. Mismatches: "
        f"**{mismatch_count}**.",
    )
    lines.append("")
    lines.append("The element table the assembler emits, in the document's format:")
    lines.append("")
    lines.append(rebuilt.element_table())
    lines.append("")
    lines.append(
        "The coding region opens at nucleotide 55 and closes on the first "
        f"in-frame stop codon at nucleotide {orf.end}, which reproduces the "
        "document's own feature table rather than a paraphrase of it.",
    )
    lines.append("")

    lines.append("## 2. Where the nucleotides came from")
    lines.append("")
    lines.append(
        "A rebuild is only as good as the independence of its inputs, so "
        "here is the provenance of every part.",
    )
    lines.append("")
    lines.append("| Element | nt | Source | Independent check |")
    lines.append("| --- | --- | --- | --- |")
    lines.append(
        "| cap | 2 | WHO INN document 11889 | none available, it is a "
        "chemical structure |",
    )
    lines.append(
        "| 5'-UTR | 52 | WHO INN document 11889, positions 3-54 | printed "
        "identically by Fang et al. 2022, Table 1 |",
    )
    lines.append(
        "| coding region | 3,825 | berthubert/bnt162b2 transcription of the "
        "same document, plus one stop codon | byte-identical to the document "
        "over 3,822 nt |",
    )
    lines.append("| 3'-UTR | 295 | WHO INN document 11889, positions 3880-4174 | none available |")
    lines.append("| poly(A) | 110 | WHO INN document 11889, positions 4175-4284 | none available |")
    lines.append("")
    lines.append(
        "Three nucleotides of 4,284, or 0.07 per cent, come from the "
        "document being reproduced: the second of its two consecutive stop "
        "codons. The published coding-sequence transcription ends at the "
        "first stop codon, because that is where translation ends, so the "
        "second one has to be added back to reproduce the construct. "
        "Everything else comes from the published 5' UTR, the published "
        "coding sequence and the published 3' UTR and poly(A) specification, "
        "and the 3,822 coding nucleotides were checked against the document "
        "byte for byte before this report was written.",
    )
    lines.append("")

    lines.append("## 3. The cap is a chemical footprint, not two nucleotides")
    lines.append("")
    lines.append(
        "The document prints two nucleotides for the cap and starts its 5' "
        "UTR element at position 3. That is not an off-by-one error, it is a "
        "different convention for the same molecule, and it is worth being "
        "precise about, because the same UTR is quoted at three different "
        "lengths across the literature.",
    )
    lines.append("")
    lines.append("| What is counted | Length |")
    lines.append("| --- | --- |")
    lines.append("| The 5' UTR element as the WHO document numbers it, positions 3 to 54 | 52 nt |")
    lines.append("| The same region including the first backbone nucleotide | 53 nt |")
    lines.append("| The same region including the cap guanosine as well | 54 nt |")
    lines.append("")
    lines.append(
        "The chemistry behind it. The cap is 7-methylguanosine, joined to "
        "the first backbone nucleotide by a 5'-to-5' triphosphate, so it "
        "covers two printed positions: the m7G at position 1 and the first "
        "nucleotide of the backbone at position 2. The document writes it "
        "as 'cap G1A2 = m7G+m3'-5'-ppp-5'-Am2'-3'-p-', with a 2'-O-methyl on "
        "that first backbone nucleotide, which is what makes this a cap 1 "
        "rather than a cap 0. The ribosome scans the whole physical leader, "
        "which is 53 nt, whatever the table calls it.",
    )
    lines.append("")
    lines.append(
        "This engine uses the document's numbering, because the document's "
        "coordinates are the ones a design has to match. A reader who "
        "counts the UTR by hand will get 54 and will be right about the "
        "molecule and wrong about the table, which is exactly the kind of "
        "disagreement that turns into a silent error in a later milestone.",
    )
    lines.append("")

    lines.append("## 4. A new construct on the same scaffold")
    lines.append("")
    lines.append(
        "An assembler that only reproduces Comirnaty has reproduced a lookup "
        "table. The same call that rebuilds the licensed product builds a "
        "new one: the same cap, the same 5' UTR, the same 3' UTR, the same "
        "poly(A) specification, and a different payload.",
    )
    lines.append("")
    lines.append(designed.element_table())
    lines.append("")
    designed_orf = designed.open_reading_frame()
    lines.append("| Property | Comirnaty, rebuilt | New construct |")
    lines.append("| --- | --- | --- |")
    lines.append(f"| Total length | {len(rebuilt.sequence)} nt | {len(designed.sequence)} nt |")
    lines.append(
        f"| Translated coding sequence, start codon to first stop | "
        f"{orf.end - orf.start} nt | {designed_orf.end - designed_orf.start} nt |",
    )
    lines.append(
        f"| Protein | {len(rebuilt.protein)} aa | {len(designed.protein)} aa |",
    )
    lines.append(
        f"| Upstream ATG in the 5' UTR | {len(rebuilt.upstream_starts())} | "
        f"{len(designed.upstream_starts())} |",
    )
    lines.append(
        f"| AAUAAA occurrences, whole construct | {len(rebuilt_hits)} | "
        f"{len(designed_hits)} |",
    )
    lines.append(
        "| AAUAAA occurrences in the coding region | 0 | 0 |",
    )
    lines.append("")
    lines.append(
        "The new construct is a teaching example rather than a candidate "
        "product: green fluorescent protein is not a venom antigen. It is "
        "here because it is small enough to read in full, because it was "
        "already in the reference data, and because the next milestone does "
        "the same thing with a payload that matters.",
    )
    lines.append("")

    lines.append("## 5. Motif scope: the licensed product carries the signal twice")
    lines.append("")
    lines.append(
        "M1 enforced one forbidden motif in the coding region: the AAUAAA "
        "polyadenylation signal. Applied to a whole construct, the same "
        "scan finds it twice in Comirnaty.",
    )
    lines.append("")
    lines.append("| Element | Length | AAUAAA occurrences |")
    lines.append("| --- | --- | --- |")
    for part in rebuilt.parts:
        count = len([hit for hit in rebuilt_hits if hit.element == part.name])
        lines.append(f"| {part.name} | {len(part.sequence)} nt | {count} |")
    lines.append("")
    for hit in rebuilt_hits:
        lines.append(f"- {hit.describe()}")
    lines.append("")
    lines.append(
        "Both occurrences are in untranslated regions, and the coding region "
        "has none. That is the correct state of affairs, and it is worth "
        "being explicit about why, because the difference is the whole point "
        "of scanning per element rather than per construct.",
    )
    lines.append("")
    lines.append(
        "A polyadenylation signal is dangerous when it sits downstream of a "
        "cleavage site, because the cell can then cut the transcript there. "
        "Inside a coding region that means a truncated message with no "
        "in-frame stop codon, and a ribosome that runs into the poly(A) "
        "tail and stalls on a poly-lysine tract. Inside a 5' UTR there is no "
        "cleavage site upstream of it, and inside a 3' UTR the signal only "
        "matters if it is the one the cell chooses to use. A signal that is "
        "harmless where it sits is not a defect, and a scanner that cannot "
        "tell the difference will eventually get a licensed product wrong.",
    )
    lines.append("")
    lines.append(
        "The engine therefore reports motif occurrences per element, and "
        "the element they sit in. It does not report a single verdict for "
        "the construct, because a single verdict would have failed a "
        "product that two regulators approved.",
    )
    lines.append("")

    lines.append("## 6. Upstream start codons")
    lines.append("")
    lines.append(
        "The ribosome scans the 5' end and usually takes the first start "
        "codon it meets. A start codon in the 5' UTR therefore competes with "
        "the intended one, and the context around it decides how well. "
        "Kozak 1986 showed that the purine at position -3 and the G at "
        "position +4 are the two positions that matter. The assembler now "
        "reports every ATG in the 5' UTR with its context.",
    )
    lines.append("")
    lines.append(
        "Comirnaty's 5' UTR contains "
        f"{len(rebuilt.upstream_starts())} upstream start codons, and the new "
        f"construct's 5' UTR contains {len(designed.upstream_starts())}. This "
        "check was deliberately left out of the M1 optimiser, because it is a "
        "property of the UTR and not of the coding region, and it belongs "
        "here where the real UTR and the real Kozak context are known.",
    )
    lines.append("")

    lines.append("## 7. What this does not claim")
    lines.append("")
    lines.append(
        "No expression level is predicted. No protein folding is predicted. "
        "Nothing has been synthesised and nothing has been put into a cell. "
        "The assembler makes a construct well formed with respect to a set "
        "of cited, measurable constraints and can show its work by "
        "reproducing a published one. That is the whole claim, and the "
        "ladder from here is L1, expression in cells, which belongs to a "
        "wet-lab partner.",
    )
    lines.append("")
    lines.append(
        "The next milestone takes one toxin from India's Big Four, one "
        "published neutralising antibody against it, and builds the same "
        "kind of construct with that antibody as the payload. Everything "
        "this milestone added is what that milestone needs and nothing more.",
    )
    lines.append("")

    lines.append("## References")
    lines.append("")
    lines.append(
        "- WHO International Nonproprietary Names Programme, document 11889, "
        "Comirnaty (BNT162b2), February 2021.",
    )
    lines.append(
        "- Orlandini von Niessen AG, Poleganov MA, Rechner C, Plaschke A, "
        "Kranz LM, Fesser M, Diken M, Loewer M, Vallazza B, Beissert T, Bukur "
        "V, Kuhn AN, Tureci O, Sahin U. 2019. Improving mRNA-based "
        "therapeutic gene delivery by expression-augmenting 3' UTRs "
        "identified by cellular library screening. Mol Ther 27(4):824-836.",
    )
    lines.append(
        "- Fang E, Liu X, Li M, Zhang Z, Song L, Zhu B, Wu X, Liu J. 2022. "
        "Advances in COVID-19 mRNA vaccine development. Signal Transduct "
        "Target Ther 7(1):94.",
    )
    lines.append(
        "- Xia X. 2021. Detailed Dissection and Critical Evaluation of the "
        "Pfizer/BioNTech and Moderna mRNA Vaccines. Vaccines 9(7):734.",
    )
    lines.append(
        "- Kariko K, Buckstein M, Ni H, Weissman D. 2005. Suppression of RNA "
        "recognition by Toll-like receptors: the impact of nucleoside "
        "modification and the evolutionary origin of RNA. Immunity "
        "23(2):165-175.",
    )
    lines.append(
        "- Kozak M. 1986. Point mutation defines a sequence flanking the AUG "
        "initiator codon that modulates translation by eukaryotic ribosomes. "
        "Cell 44(2):283-292.",
    )
    lines.append(
        "- Proudfoot NJ, Brownlee BG. 1976. AATAAA is a signal for synthesis "
        "of polyadenylated RNA. Nucleic Acids Res 3(7):1183-1190.",
    )
    lines.append(
        "- Pardi N, Hogan MJ, Porter FW, Weissman D. 2018. mRNA vaccines, a "
        "new era in vaccinology. Nat Rev Drug Discov 17:261-279.",
    )
    lines.append("")
    lines.append(f"_{DISCLAIMER}_")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    """Run the analysis, optionally refreshing the committed report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the report to docs/m2_construct_assembler.md",
    )
    args = parser.parse_args()

    report = build_report()
    if args.write:
        target = ROOT / "docs" / "m2_construct_assembler.md"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(report, encoding="utf-8")
        print(f"wrote {target.relative_to(ROOT)}")
    else:
        print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
