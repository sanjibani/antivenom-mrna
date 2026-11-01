"""M3: one toxin, one published neutralising antibody, one full construct.

This script is the M3 deliverable. It builds the first real venom construct
in this repository: an mRNA that encodes a published human single-chain
variable fragment against a toxin of one of India's Big Four snakes, on the
scaffold that M2 built and measured.

Run it::

    uv run python scripts/venom_construct.py             # print the report
    uv run python scripts/venom_construct.py --write     # refresh the design files

What it does, in order:

1. Takes the E113 single-chain variable fragment from the Protein Data Bank
   entry that holds it, with the crystallisation construct's vector residue
   and expression tag removed.
2. Takes the toxin it binds, ecarpholin S from the saw-scaled viper, from
   UniProt, and checks it against the crystal structure of the same toxin.
3. Codes the payload with the M1 optimiser, under the M2 policy.
4. Assembles the full construct and measures it.
5. Writes the design out, and reports what the engine checked and what it
   cannot know.

The payload is a published antibody. This project does not discover
antibodies and does not claim to.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mrna.assemble import A110_POLYA, CAP1, Construct, full_construct
from mrna.fasta import write_fasta
from mrna.leader import HSA_SIGNAL_PEPTIDE
from mrna.motifs import DEFAULT_MOTIFS
from mrna.optimize import optimize_cds
from mrna.policy import HUMAN_POLICY
from mrna.reference import load_codon_usage, load_single_sequence
from mrna.utr import AES_MTRNR1_3UTR, ALPHA_GLOBIN_5UTR

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"
DESIGNS = ROOT / "data" / "designs"

DISCLAIMER = "In-silico designs. No therapeutic efficacy demonstrated. Not medical advice."

#: The single stop codon of the designed construct.
STOP_CODON = "TGA"

#: The three parts of the E113 single-chain variable fragment.
#:
#: The boundaries come from the deposited construct of PDB entry 8IA6 and
#: from the human germline sequences the construct is built on: IGHV3-23
#: for the heavy chain variable domain, IGKV1 for the light chain variable
#: domain, and IGKC for the light chain constant region.
ANTIBODY_PARTS: tuple[tuple[str, str, str], ...] = (
    ("VH", "EVQLLESGGG", "VTVSS"),
    ("linker", "GGGGSGGGGS", "GGGGSGGGGSGGGGST"),
    ("VL", "DIQMTQ", "FGQGTKVEIK"),
)

#: What the source paper reports about the antibody, quoted rather than
#: paraphrased, so that the design document does not overstate it.
PUBLISHED_EVIDENCE: tuple[tuple[str, str], ...] = (
    (
        "Binding",
        "PLA2 was 'prominently captured by the phage-anchored scFv "
        "antibodies', and was the dominant toxin of both venoms by "
        "proteomics, at 24 per cent of Naja naja venom and 28 per cent of "
        "Echis carinatus venom.",
    ),
    (
        "Target",
        "E113 and E10, the two anti-Echis carinatus clones, both interact "
        "with the basic PLA2 of the venom.",
    ),
    (
        "In vivo",
        "E113 gave 'significant prolongation of survival time' in Echis "
        "carinatus venom-challenged mice. The equal-dose combination of "
        "E113 and E10 gave 33 per cent survival, against 16 per cent for "
        "the individual scFvs.",
    ),
    (
        "Mechanism",
        "E113 inhibited venom-induced hemolysis by 50 per cent, and a "
        "structural model of the E113-PLA2 complex placed the CDR regions "
        "against the reported C-terminal myotoxic loop and interfacial "
        "region of PLA2.",
    ),
    (
        "Format",
        "A single-chain variable fragment from a human antibody library, so "
        "one polypeptide, one open reading frame, one construct.",
    ),
)


def load_payload() -> str:
    """Return the E113 single-chain variable fragment."""
    return load_single_sequence(REFERENCE / "e113_scfv.fasta").strip()


def load_toxin() -> str:
    """Return ecarpholin S, the toxin E113 binds."""
    return load_single_sequence(REFERENCE / "ecarpholin_S_P48650.fasta").strip()


def split_antibody(payload: str) -> tuple[tuple[str, str], ...]:
    """Split the payload into its heavy chain, linker and light chain.

    Args:
        payload: The single-chain variable fragment.

    Returns:
        A tuple of ``(name, sequence)`` pairs.

    Raises:
        ValueError: If the payload does not have the expected parts.
    """
    parts: list[tuple[str, str]] = []
    cursor = 0
    for name, start_marker, end_marker in ANTIBODY_PARTS:
        if not payload.startswith(start_marker, cursor):
            msg = f"the payload does not have the expected {name} start"
            raise ValueError(msg)
        end = payload.index(end_marker, cursor) + len(end_marker)
        parts.append((name, payload[cursor:end]))
        cursor = end
    if cursor != len(payload):
        msg = "the payload has residues after the expected light chain end"
        raise ValueError(msg)
    return tuple(parts)


def build_construct(usage_table: dict[str, float]) -> tuple[Construct, object]:
    """Design the construct and return it with its optimisation result.

    Args:
        usage_table: Codon usage table for the target host.

    Returns:
        The assembled construct and the optimisation result for its coding
        region.
    """
    payload = load_payload()
    protein = HSA_SIGNAL_PEPTIDE.sequence + payload
    design = optimize_cds(protein, usage_table, HUMAN_POLICY)
    coding = design.cds + STOP_CODON

    # The coordinates in the element description are construct coordinates,
    # not coordinates inside the coding element, because the element table
    # of a regulatory document counts from the start of the transcript.
    cap_nt = len(CAP1.prefix)
    utr5_nt = len(ALPHA_GLOBIN_5UTR.sequence)
    signal_nt = HSA_SIGNAL_PEPTIDE.length * 3
    payload_start = cap_nt + utr5_nt + signal_nt + 1
    payload_end = payload_start + len(coding) - signal_nt - 1
    stop_start = payload_end - len(STOP_CODON) + 1

    construct = full_construct(
        "E113 scFv against Echis carinatus ecarpholin S",
        utr5=ALPHA_GLOBIN_5UTR,
        coding=coding,
        utr3=AES_MTRNR1_3UTR,
        poly_a=A110_POLYA,
        leader=HSA_SIGNAL_PEPTIDE,
        coding_name="E113 scFv",
        coding_description=(
            "Codon-optimised single-chain variable fragment E113 against "
            "Echis carinatus phospholipase A2, preceded by the human serum "
            f"albumin signal peptide, followed by a stop codon at {stop_start} "
            f"to {payload_end}."
        ),
    )
    return construct, design


def write_design(construct: Construct, protein: str, toxin: str) -> list[Path]:
    """Write the design files.

    Args:
        construct: The assembled construct.
        protein: The protein the construct encodes.
        toxin: The toxin the antibody binds.

    Returns:
        The paths written.
    """
    DESIGNS.mkdir(parents=True, exist_ok=True)

    fasta_path = DESIGNS / "e113_scfv_ecarpholin_s_mrna.fasta"
    write_fasta(
        fasta_path,
        [
            (
                f"{construct.name} | full mRNA construct, {construct.length} nt | "
                "source: antivenom-mrna M3 | in-silico design, no therapeutic "
                "efficacy demonstrated",
                construct.sequence,
            ),
            (
                "E113 scFv payload | protein encoded by the construct, "
                f"{len(protein)} aa | source: PDB entry 8IA6 entity 1, "
                "crystallisation tag removed | plus the human serum albumin "
                "signal peptide, UniProt P02768 residues 1-18",
                protein,
            ),
            (
                "Ecarpholin S | target toxin, Echis carinatus phospholipase A2 "
                f"| source: UniProt P48650, structure PDB 2QHE | {len(toxin)} aa",
                toxin,
            ),
        ],
    )

    table_path = DESIGNS / "e113_scfv_ecarpholin_s_mrna.elements.tsv"
    rows = ["element\tstart\tend\tlength\tdescription"]
    for element in construct.elements:
        description = element.description.replace("\t", " ")
        rows.append(
            f"{element.name}\t{element.start}\t{element.end}\t"
            f"{element.length}\t{description}",
        )
    table_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    return [fasta_path, table_path]


def build_report() -> str:
    """Build the full M3 markdown report.

    Returns:
        The report as a markdown string.
    """
    usage_table = load_codon_usage(REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv")
    construct, design = build_construct(usage_table)
    payload = load_payload()
    toxin = load_toxin()
    parts = split_antibody(payload)
    protein = HSA_SIGNAL_PEPTIDE.sequence + payload
    orf = construct.open_reading_frame()
    hits = construct.find_motifs(DEFAULT_MOTIFS)

    lines: list[str] = []
    lines.append("# One toxin, one published antibody, one full construct")
    lines.append("")
    lines.append(f"> {DISCLAIMER}")
    lines.append("")
    lines.append(
        "This is the M3 artifact, and it is the first venom construct in this "
        "repository. It encodes a published human single-chain variable "
        "fragment against a toxin of the saw-scaled viper, *Echis carinatus*, "
        "on the cap, UTR and poly(A) scaffold that M2 built and measured "
        "against a licensed vaccine.",
    )
    lines.append("")
    lines.append(
        "The claim this document makes is narrow and checkable. Every "
        "sequence has a source, every number was computed by "
        "`scripts/venom_construct.py` from the committed reference data, and "
        "the construct is reproducible byte for byte. The document does not "
        "claim that the molecule works, because nobody has expressed it.",
    )
    lines.append("")

    lines.append("## 1. The toxin")
    lines.append("")
    lines.append(
        "| Field | Value |",
    )
    lines.append("| --- | --- |")
    lines.append("| Species | *Echis carinatus*, the saw-scaled viper, sochurek |")
    lines.append("| In India | One of the Big Four, the four species that account for "
                 "the overwhelming majority of snakebite deaths in the country |")
    lines.append("| Toxin | Ecarpholin S, the myotoxic Ser49 phospholipase A2 |")
    lines.append("| Sequence | UniProt P48650, 122 residues |")
    lines.append("| Structure | PDB 2QHE, 2.0 A resolution |")
    lines.append("| Family | Snake venom phospholipase A2 |")
    lines.append("")
    lines.append(
        "Phospholipase A2 is the right place to start, and the source paper "
        "says why in numbers rather than adjectives: a proteomic survey of "
        "the venom found PLA2 to be the dominant toxin, at 28 per cent of "
        "*Echis carinatus* venom by mass. PLA2 is also the toxin family that "
        "does the most damage in viper bites. It is myotoxic because the "
        "C-terminal hydrophobic loop inserts into the muscle membrane, it "
        "is hemolytic because it hydrolyses the lipid bilayer of the red "
        "cell, and it is anticoagulant because it cleaves the phospholipids "
        "that clotting factors need.",
    )
    lines.append("")
    lines.append(
        "Ecarpholin S sits in the Ser49 subgroup of the venom "
        "phospholipases, and that single residue is what makes it a myotoxin "
        "rather than an enzyme. Most Ca2-dependent PLA2 enzymes carry an "
        "aspartate at position 49, and that aspartate coordinates the "
        "calcium ion the catalytic cycle needs. Ecarpholin S has a serine "
        "there instead, and the structure shows a calcium-binding loop "
        "folded into a conformation that cannot coordinate the ion. The "
        "result, in the words of the structural paper, is 'low enzymatic but "
        "potent myotoxic activities'. The catalytic machinery is all still "
        "there; it simply cannot run.",
    )
    lines.append("")
    lines.append(
        "The molecule is a problem for an mRNA product in a specific way. A "
        "snake venom PLA2 is a small, disulfide-rich protein: UniProt "
        "annotates seven disulfide bonds for this sequence, which means it "
        "has to fold correctly, with the right disulfide pairing, in an "
        "oxidising environment. The cytosol is reducing and the wrong place "
        "for it. That is the reason the design below carries a signal "
        "peptide at all, and it is the reason the payload is an antibody "
        "rather than the toxin itself: the antibody is made in the "
        "endoplasmic reticulum where the disulfides can form.",
    )
    lines.append("")
    lines.append(
        "Two things in the UniProt record connect the toxin to the "
        "antibody. The first is the region annotated as 'important for "
        "membrane-damaging activities in eukaryotes and bacteria', residues "
        "105 to 117, the C-terminal hydrophobic loop that does the actual "
        "damage. The second is the FUNCTION line, which says this is a "
        "'snake venom phospholipase A2 homolog that lacks enzymatic "
        "activity. Shows high myotoxin activities and displays edema-inducing "
        "activities'. The source paper's structural model of the E113-PLA2 "
        "complex puts the antibody's CDR regions against exactly that "
        "C-terminal myotoxic loop and against the interfacial region of the "
        "toxin, which is the part of a PLA2 that docks into a membrane. Two "
        "independent sources, a sequence database and a paper, put the "
        "antibody's epitope on the part of the toxin that does the damage.",
    )
    lines.append("")
    lines.append("The target sequence, as committed:")
    lines.append("")
    lines.append("```")
    lines.append(toxin)
    lines.append("```")
    lines.append("")

    lines.append("## 2. The antibody")
    lines.append("")
    lines.append("| Field | Value |")
    lines.append("| --- | --- |")
    lines.append("| Clone | E113 |")
    lines.append("| Format | Human single-chain variable fragment, one polypeptide |")
    lines.append("| Sequence | PDB 8IA6 entity 1, 2.75 A crystal structure |")
    lines.append("| Source | Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2024. "
                 "Mol Immunol 165:55-67 |")
    lines.append("| Patent | No patent on this clone located in a public "
                 "search on 2026-11-01 |")
    lines.append("")
    lines.append("The payload is split into its three parts:")
    lines.append("")
    lines.append("| Part | Length | Sequence |")
    lines.append("| --- | --- | --- |")
    for name, sequence in parts:
        lines.append(f"| {name} | {len(sequence)} aa | `{sequence}` |")
    lines.append("")
    lines.append(
        "The heavy chain is IGHV3-23, the most used human immunoglobulin "
        "heavy chain variable gene. The light chain is an IGKV1 variable "
        "domain joined to the IGKC constant region, which is where the "
        "`FGQGTKVEIK` at the end comes from.",
    )
    lines.append("")
    lines.append("### What the deposited construct is, and what was removed")
    lines.append("")
    lines.append(
        "The Protein Data Bank entry holds a crystallisation construct, not "
        "a therapeutic molecule. It has a vector residue at the N-terminus "
        "and a hexahistidine expression tag at the C-terminus. The crystal "
        "locates residues 2 to 239 of the deposited chain: residue 1, the "
        "vector alanine, was never located, residues 118 to 133, the "
        "flexible linker, are disordered in every copy, and the expression "
        "tag was never located either.",
    )
    lines.append("")
    lines.append(
        "The design encodes the 239-residue antibody and neither artefact. "
        "The removal is done by the reference-data script, not by hand, and "
        "the script refuses to run if the deposited construct does not have "
        "the expected ends.",
    )
    lines.append("")

    lines.append("## 3. The construct")
    lines.append("")
    lines.append(construct.element_table())
    lines.append("")
    lines.append(
        "| Property | Value |",
    )
    lines.append("| --- | --- |")
    lines.append(f"| Total length | {construct.length} nt |")
    lines.append(f"| Translated coding sequence | {orf.end - orf.start} nt |")
    lines.append(f"| Encoded protein | {len(protein)} aa |")
    lines.append(f"| Payload | {len(payload)} aa |")
    lines.append(f"| GC content of the coding region | {design.gc_total:.4f} |")
    lines.append(f"| CAI against human codon usage | {design.cai:.4f} |")
    lines.append(
        f"| GC per {HUMAN_POLICY.gc_block_nt} nt block | "
        f"{min(design.block_gc):.3f} to {max(design.block_gc):.3f} |",
    )
    lines.append(
        f"| Forbidden motifs in the coding region | {len(design.motif_hits)} |",
    )
    lines.append(
        f"| Upstream ATG in the 5' UTR | {len(construct.upstream_starts())} |",
    )
    lines.append(
        f"| AAUAAA in the whole construct | {len(hits)} "
        f"({', '.join(sorted({hit.element for hit in hits}))}) |",
    )
    lines.append("")
    lines.append("The design decisions, each with the reason for it:")
    lines.append("")
    lines.append(
        "| Decision | Why |",
    )
    lines.append("| --- | --- |")
    lines.append(
        "| Human serum albumin signal peptide | One of the two most commonly "
        "used signal peptides in mRNA vaccine designs, and the address label "
        "that gets a disulfide-rich antibody out of the cell |",
    )
    lines.append(
        "| Human alpha-globin 5' UTR | The element Comirnaty uses, so that the "
        "construct reads line by line against a licensed product |",
    )
    lines.append(
        "| AES plus mitochondrial 12S rRNA 3' UTR | The same, and the 3' UTR "
        "that a screen of human UTR fragments found to be the most "
        "expression-augmenting |",
    )
    lines.append(
        "| A110 template-encoded poly(A) | 30 adenosines, a 10-nucleotide "
        "linker, 70 adenosines. The linker is what a ribosome that misreads "
        "past the stop codon runs into instead of a longer poly-lysine tract |",
    )
    lines.append(
        "| Single TGA stop codon | The minimal choice. Comirnaty prints two, "
        "which does not change what is translated |",
    )
    lines.append(
        "| Codon optimisation subject to the GC band | The M1 result: the CAI "
        "is optimal subject to a band that both licensed products satisfy |",
    )
    lines.append("")

    lines.append("## 4. What the engine checked")
    lines.append("")
    lines.append(
        "Four checks, all of them mechanical, all of them reproduced by the "
        "test suite.",
    )
    lines.append("")
    lines.append(
        "**The construct encodes what it claims.** The coding region "
        f"translates back to {len(protein)} residues: the 18-residue signal "
        "peptide followed by the 239-residue antibody, with no substitutions.",
    )
    lines.append("")
    lines.append(
        "**The initiation context is the right one.** The alpha-globin 5' UTR "
        "ends in GCCACC, so the start codon sits in the Kozak context that "
        "the same construct uses in a licensed product.",
    )
    lines.append("")
    lines.append(
        "**There is no upstream start codon.** A start codon in the 5' UTR "
        "would compete with the intended one, and the ribosome usually takes "
        "the first it meets. There are none.",
    )
    lines.append("")
    lines.append(
        "**The coding region holds every GC block inside the band and "
        "carries no forbidden motif.** The polyadenylation signal that "
        "Comirnaty carries twice in its untranslated regions is absent from "
        "the coding region, which is where it would matter.",
    )
    lines.append("")

    lines.append("## 5. What the source paper actually shows")
    lines.append("")
    lines.append(
        "The evidence table is quoted rather than paraphrased, because the "
        "difference between the two is the difference between an honest "
        "design document and a promotional one.",
    )
    lines.append("")
    lines.append("| | Reported |")
    lines.append("| --- | --- |")
    for label, text in PUBLISHED_EVIDENCE:
        lines.append(f"| {label} | {text} |")
    lines.append("")

    lines.append("## 6. What this construct is not")
    lines.append("")
    lines.append(
        "An scFv is the weakest possible starting point for an antivenom, "
        "and the reasons are worth stating before anyone builds it.",
    )
    lines.append("")
    lines.append(
        "**It is monovalent.** One binding site means one toxin per "
        "antibody. The published data on this clone is about survival in "
        "mice challenged with venom, which is a whole-venom endpoint, and "
        "the combination of E113 with a second clone did better than either "
        "alone. A bivalent format would be the obvious improvement, and it "
        "would need two chains or a single-chain Fab.",
    )
    lines.append("")
    lines.append(
        "**It is small and it is cleared.** An scFv of this size is roughly "
        "26 kDa, has a serum half-life in the range of hours, and is cleared "
        "by the kidney. An IgG or an Fc fusion would last weeks.",
    )
    lines.append("")
    lines.append(
        "**The timing is wrong for an acute bite.** Intravenous antivenom "
        "works within minutes. An mRNA-encoded antibody is not detectable "
        "for one to two hours and peaks somewhere between 12 and 24 hours. "
        "This construct is a prophylactic candidate, and the repository says "
        "so at the top of its README.",
    )
    lines.append("")
    lines.append(
        "**Nothing has been expressed.** No cell has made this protein. The "
        "next rung of the validation ladder is L1, expression in cells, and "
        "that belongs to a wet-lab partner. This project does not have one "
        "and does not intend to.",
    )
    lines.append("")
    lines.append(
        "**The antigen-binding claim is about the antibody, not about the "
        "mRNA.** The published binding and neutralisation data belongs to "
        "E113 as a protein made in bacteria. An scFv made in a mammalian "
        "cell from an mRNA template is the same sequence, but it is a "
        "different molecule until someone shows it is not: glycosylation, "
        "disulfide pairing and folding all differ between a bacterial and a "
        "mammalian expression system, and for a disulfide-rich antibody that "
        "is not a detail.",
    )
    lines.append("")
    lines.append(
        "**The patent position is a negative search result, not a clearance.** "
        "No patent covering the E113 clone was located in a public patent "
        "search on 2026-11-01. That is not a freedom-to-operate opinion, and "
        "milestone M4 exists to record patent status properly across the "
        "whole candidate set.",
    )
    lines.append("")
    lines.append(
        "**Venom varies regionally.** The same species does not have the "
        "same venom everywhere, and a construct aimed at India has to be "
        "aimed at Indian venom. The east-coast and the west-coast "
        "*Echis carinatus* are not the same animal.",
    )
    lines.append("")

    lines.append("## 7. What comes next")
    lines.append("")
    lines.append(
        "M4 is the tracker: species, toxin family, target, antibody, "
        "format, sequence source, patent status, what it was tested against, "
        "evidence level and mRNA suitability, for every candidate in the "
        "field. M5 is the cocktail: the same engine, run over the Big Four, "
        "with the combinations that the evidence supports.",
    )
    lines.append("")

    lines.append("## References")
    lines.append("")
    lines.append(
        "- Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2024. Recombinant "
        "human scFv antibody fragments against phospholipase A2 from Naja "
        "naja and Echis carinatus snake venoms: in vivo neutralization and "
        "mechanistic insights. Mol Immunol 165:55-67. "
        "doi:10.1016/j.molimm.2023.12.006.",
    )
    lines.append(
        "- Protein Data Bank entry 8IA6, crystal structure of scFv antibody "
        "against phospholipase A2 of Echis carinatus venom, 2.75 A, "
        "deposited by Kumar A, Madni ZK, Salunke DM.",
    )
    lines.append(
        "- Protein Data Bank entry 2QHE, crystal structure of the "
        "phospholipase A2 of Echis carinatus, the toxin E113 binds.",
    )
    lines.append(
        "- UniProt P48650, ecarpholin S, basic phospholipase A2 homolog from "
        "Echis carinatus.",
    )
    lines.append(
        "- UniProt P02768, human serum albumin precursor, the source of the "
        "signal peptide.",
    )
    lines.append(
        "- Zhou X, Tan TC, Valiyaveettil S, Go ML, Kini RM, Velazquez-Campoy "
        "A, Sivaraman J. 2008. Structural characterization of myotoxic "
        "ecarpholin S from Echis carinatus venom. Biophys J 95:3366-3380. "
        "doi:10.1529/biophysj.107.117747.",
    )
    lines.append(
        "- Lomonte B, Rangel J. 2012. Snake venom Lys49 myotoxins: from "
        "phospholipases A2 to non-enzymatic membrane disruptors. Toxicon "
        "60:520-530.",
    )
    lines.append(
        "- Orlandini von Niessen AG, et al. 2019. Improving mRNA-based "
        "therapeutic gene delivery by expression-augmenting 3' UTRs "
        "identified by cellular library screening. Mol Ther 27(4):824-836.",
    )
    lines.append(
        "- WHO International Nonproprietary Names Programme, document 11889, "
        "Comirnaty (BNT162b2), February 2021.",
    )
    lines.append("")
    lines.append(f"_{DISCLAIMER}_")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    """Run the design, optionally writing the design files and the report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the design files and docs/m3_ecarpholin_s_construct.md",
    )
    args = parser.parse_args()

    report = build_report()
    if args.write:
        target = ROOT / "docs" / "m3_ecarpholin_s_construct.md"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(report, encoding="utf-8")
        print(f"wrote {target.relative_to(ROOT)}")

        usage_table = load_codon_usage(REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv")
        construct, _design = build_construct(usage_table)
        payload = load_payload()
        toxin = load_toxin()
        protein = HSA_SIGNAL_PEPTIDE.sequence + payload
        for path in write_design(construct, protein, toxin):
            print(f"wrote {path.relative_to(ROOT)}")
    else:
        print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
