"""M6: two payloads, one transcript.

This script is the M6 deliverable. It answers the design question M5
deferred: whether several payloads can share one transcript. The answer is
yes, on the record, at the standard this repository holds everything else
to. The 2A route that M5 could not source is sourced now: the porcine
teschovirus-1 2A peptide is read out of two deposited constructs from two
different groups, and the two agree with each other and with the
characterisation paper.

Run it::

    uv run python scripts/cocktail_transcript.py             # print the report
    uv run python scripts/cocktail_transcript.py --write     # refresh the design files

What it does, in order:

1. Reads the P2A peptide from the two deposited constructs in
   ``data/reference/raw/`` and checks both against the engine's constant.
2. Assembles the ORF protein: the human serum albumin signal peptide, the
   TPL1158_01_C09 nanobody, the GSG-P2A skip element, a second signal
   peptide preceded by the carry-over proline, and the E113 single-chain
   variable fragment.
3. Codes the whole open reading frame with the M1 optimiser, under the M2
   policy, and slices it back into elements by residue counts.
4. Assembles the full construct, measures it, and checks the two mature
   products the cell will actually make: the upstream nanobody with its
   P2A scar, and the downstream scFv, which comes out exactly as published.
5. Writes the design out, and reports what the engine checked and what it
   cannot know.

The order of the payloads is the design decision. The ribosome skips
between the last glycine of the 2A element and the proline that follows it,
so the upstream product carries a 21-residue C-terminal scar that the
published TPL1158_01_C09 data does not describe, while the downstream
product starts with one proline that sits in front of a signal peptide,
which signal peptidase removes together with the proline. The downstream
payload is therefore reproduced residue for residue, and the upstream
payload is not. E113 carries the stronger published evidence, in vivo
neutralisation at L4 against binding only at L2 for the nanobody, so E113
goes downstream.

The licence position is unchanged from M5 and restated in the report: a
public patent search finds an application covering TPL1158_01_C09, and an
open programme built on that sequence is not free to exist. This
construct is a design artifact, not a freedom-to-operate opinion.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mrna.assemble import (
    A110_POLYA,
    CAP1,
    Construct,
    ConstructPart,
    assemble_construct,
)
from mrna.fasta import read_fasta, write_fasta
from mrna.leader import HSA_SIGNAL_PEPTIDE
from mrna.motifs import DEFAULT_MOTIFS
from mrna.optimize import OptimizationResult, optimize_cds
from mrna.p2a import GSG_P2A_LINKER, P2A, P2A_DOWNSTREAM_PROLINE
from mrna.payload import load_payloads
from mrna.policy import HUMAN_POLICY
from mrna.reference import load_codon_usage
from mrna.sequence import translate
from mrna.utr import AES_MTRNR1_3UTR, ALPHA_GLOBIN_5UTR

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"
DESIGNS = ROOT / "data" / "designs"
RAW = REFERENCE / "raw"

REPORT_PATH = ROOT / "docs" / "m6_single_transcript.md"
FASTA_PATH = DESIGNS / "tpl1158_01_c09_p2a_e113_cocktail_mrna.fasta"
TABLE_PATH = DESIGNS / "tpl1158_01_c09_p2a_e113_cocktail.elements.tsv"

DISCLAIMER = (
    "In-silico designs. No therapeutic efficacy demonstrated. "
    "Not medical advice."
)

STOP_CODON = "TGA"

#: The raw deposited constructs that carry the P2A peptide. The test suite
#: extracts the peptide from each of these files and checks the two agree
#: with the engine's constant, which is the M0 two-source rule applied to
#: a new element.
P2A_RAW_SOURCES: tuple[str, ...] = (
    "rcsb_31FW_p2a_construct.fasta",
    "rcsb_8E99_p2a_construct.fasta",
)

#: The elements of the open reading frame, in order, named as they appear
#: in the element table. The optimiser codes the ORF as one protein; the
#: slicing back into elements is what turns one coding sequence into a
#: regulatory-format table.
ORF_ELEMENTS: tuple[str, ...] = (
    "signal peptide",
    "TPL1158_01_C09 VHH",
    "P2A",
    "signal peptide 2",
    "E113 scFv",
)


def read_p2a_from_raw() -> dict[str, str]:
    """Extract the P2A peptide from each deposited construct.

    Returns:
        Mapping of raw file name to the peptide found in it.

    Raises:
        ValueError: If an entry does not carry the peptide, or carries a
            different one.
    """
    found: dict[str, str] = {}
    for filename in P2A_RAW_SOURCES:
        records = read_fasta(RAW / filename)
        peptide = next(
            (seq for _, seq in records if P2A.sequence in seq),
            None,
        )
        if peptide is None:
            msg = f"{filename} holds no record carrying the P2A peptide"
            raise ValueError(msg)
        found[filename] = P2A.sequence
    return found


def build_construct(
    usage_table: dict[str, float],
) -> tuple[Construct, OptimizationResult, dict[str, str]]:
    """Design the two-payload transcript and measure it.

    Args:
        usage_table: Codon usage table for the target host.

    Returns:
        The assembled construct, the optimisation result for its open
        reading frame, and the mature products the cell makes: the
        upstream nanobody with its scar, and the downstream scFv exactly
        as published.

    Raises:
        ValueError: If a design check fails.
    """
    payloads = {payload.name: payload for payload in load_payloads()}
    tpl = payloads["TPL1158_01_C09"]
    e113 = payloads["E113"]
    if not tpl.is_single_polypeptide or not e113.is_single_polypeptide:
        msg = "the cocktail transcript needs two single-chain payloads"
        raise ValueError(msg)

    # The open reading frame, as the ribosome reads it: leader, upstream
    # payload, the skip element up to and including the skip glycine, the
    # proline that begins the downstream product, the downstream leader,
    # and the downstream payload.
    upstream_mature = tpl.sequence + GSG_P2A_LINKER
    downstream_nascent = P2A_DOWNSTREAM_PROLINE + HSA_SIGNAL_PEPTIDE.sequence
    downstream_mature = e113.sequence
    orf_protein = (
        HSA_SIGNAL_PEPTIDE.sequence
        + upstream_mature
        + downstream_nascent
        + downstream_mature
    )

    design = optimize_cds(orf_protein, usage_table, HUMAN_POLICY)
    if design.motif_hits:
        msg = "the optimiser left forbidden motifs in the coding sequence"
        raise ValueError(msg)
    if translate(design.cds).rstrip("*") != orf_protein:
        msg = "the designed coding sequence does not translate to the ORF protein"
        raise ValueError(msg)

    # Slice the coding sequence back into elements, in residues.
    residue_counts = {
        "signal peptide": len(HSA_SIGNAL_PEPTIDE.sequence),
        "TPL1158_01_C09 VHH": len(tpl.sequence),
        "P2A": len(GSG_P2A_LINKER),
        "signal peptide 2": len(downstream_nascent),
        "E113 scFv": len(downstream_mature),
    }
    if sum(residue_counts.values()) != len(orf_protein):
        msg = "the element residue counts do not sum to the ORF protein"
        raise ValueError(msg)

    sequences: dict[str, str] = {}
    offset = 0
    for name in ORF_ELEMENTS:
        codons = residue_counts[name] * 3
        sequences[name] = design.cds[offset : offset + codons]
        offset += codons

    # The mature products are what the cell makes, and both are checked
    # against the published molecules by translation, not trusted.
    products = {
        "upstream_mature": (
            translate(sequences["TPL1158_01_C09 VHH"]).rstrip("*")
            + translate(sequences["P2A"]).rstrip("*")
        ),
        "downstream_mature": translate(sequences["E113 scFv"]).rstrip("*"),
    }
    if products["upstream_mature"] != upstream_mature:
        msg = "the upstream mature product does not reproduce the designed protein"
        raise ValueError(msg)
    if products["downstream_mature"] != downstream_mature:
        msg = "the downstream mature product is not the published E113 sequence"
        raise ValueError(msg)

    parts = [
        ConstructPart(
            name=CAP1.name,
            sequence=CAP1.prefix,
            description=CAP1.description,
            citation=CAP1.citation,
        ),
        ConstructPart(
            name="5'-UTR",
            sequence=ALPHA_GLOBIN_5UTR.sequence,
            description=ALPHA_GLOBIN_5UTR.description,
            citation=ALPHA_GLOBIN_5UTR.citation,
        ),
        ConstructPart(
            name="signal peptide",
            sequence=sequences["signal peptide"],
            description=HSA_SIGNAL_PEPTIDE.description,
            citation=HSA_SIGNAL_PEPTIDE.citation,
        ),
        ConstructPart(
            name="TPL1158_01_C09 VHH",
            sequence=sequences["TPL1158_01_C09 VHH"],
            description=(
                "Codon-optimised camelid VHH TPL1158_01_C09 against "
                "alpha-cobratoxin, carrying the GSG-P2A scar at its "
                "C-terminus: 21 residues the published binding data does "
                "not describe"
            ),
            citation="",
        ),
        ConstructPart(
            name="P2A",
            sequence=sequences["P2A"],
            description=P2A.description,
            citation=P2A.citation,
        ),
        ConstructPart(
            name="signal peptide 2",
            sequence=sequences["signal peptide 2"],
            description=(
                "The proline the ribosome carries over from the skip, "
                "followed by a second human serum albumin signal peptide. "
                "Signal peptidase removes both together, so the mature "
                "downstream product starts exactly at the published E113 "
                "sequence"
            ),
            citation=HSA_SIGNAL_PEPTIDE.citation,
        ),
        ConstructPart(
            name="E113 scFv",
            sequence=sequences["E113 scFv"] + STOP_CODON,
            description=(
                "Codon-optimised human scFv E113 against Echis carinatus "
                "phospholipase A2, reproduced residue for residue from the "
                "published molecule, followed by a stop codon"
            ),
            citation="",
        ),
        ConstructPart(
            name="3'-UTR",
            sequence=AES_MTRNR1_3UTR.sequence,
            description=AES_MTRNR1_3UTR.description,
            citation=AES_MTRNR1_3UTR.citation,
        ),
        ConstructPart(
            name=A110_POLYA.name,
            sequence=A110_POLYA.sequence,
            description=A110_POLYA.description,
            citation=A110_POLYA.citation,
        ),
    ]
    construct = assemble_construct(
        "TPL1158_01_C09 and E113 on one transcript, joined by the "
        "porcine teschovirus-1 P2A ribosome-skipping element",
        parts,
    )
    return construct, design, products


def write_design(
    construct: Construct,
    design: OptimizationResult,
    products: dict[str, str],
) -> list[Path]:
    """Write the design files.

    Args:
        construct: The assembled construct.
        design: The optimisation result for the open reading frame.
        products: The mature products, as :func:`build_construct` returns
            them.

    Returns:
        The paths written.
    """
    DESIGNS.mkdir(parents=True, exist_ok=True)

    write_fasta(
        FASTA_PATH,
        [
            (
                f"{construct.name} | full mRNA construct, {construct.length} nt "
                "| source: antivenom-mrna M6 | in-silico design, no "
                "therapeutic efficacy demonstrated",
                construct.sequence,
            ),
            (
                "Upstream mature product | TPL1158_01_C09 VHH with the "
                f"GSG-P2A C-terminal scar, {len(products['upstream_mature'])} "
                "aa | source: PDB entry 9GCN entity 1 plus the GSG-P2A "
                "element of src/mrna/p2a.py | the scar is 21 residues the "
                "published binding data does not describe",
                products["upstream_mature"],
            ),
            (
                "Downstream mature product | E113 scFv, "
                f"{len(products['downstream_mature'])} aa | source: PDB "
                "entry 8IA6 entity 1, crystallisation tag removed | "
                "reproduced residue for residue: the carry-over proline and "
                "the second signal peptide are removed by signal peptidase",
                products["downstream_mature"],
            ),
            (
                "Open reading frame protein | what the transcript encodes, "
                f"{len(design.protein)} aa | leader, upstream payload, "
                "GSG-P2A, carry-over proline, second leader, downstream "
                "payload",
                design.protein,
            ),
        ],
    )

    rows = ["element\tstart\tend\tlength\tdescription"]
    for element in construct.elements:
        description = element.description.replace("\t", " ")
        rows.append(
            f"{element.name}\t{element.start}\t{element.end}\t"
            f"{element.length}\t{description}",
        )
    TABLE_PATH.write_text("\n".join(rows) + "\n", encoding="utf-8")

    return [FASTA_PATH, TABLE_PATH]


def build_report() -> str:
    """Build the full M6 markdown report.

    Returns:
        The report as a markdown string.
    """
    usage_table = load_codon_usage(REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv")
    read_p2a_from_raw()
    construct, design, products = build_construct(usage_table)
    payloads = {payload.name: payload for payload in load_payloads()}
    tpl = payloads["TPL1158_01_C09"]
    d2 = payloads["D2"]
    hits = construct.find_motifs(DEFAULT_MOTIFS)
    coding_elements = set(ORF_ELEMENTS)
    coding_hits = [hit for hit in hits if hit.element in coding_elements]
    utr_hits = [hit for hit in hits if hit.element not in coding_elements]

    lines: list[str] = []
    lines.append("# M6: two payloads, one transcript")
    lines.append("")
    lines.append(f"> {DISCLAIMER}")
    lines.append("")

    lines.append("## 1. The design question M5 left open")
    lines.append("")
    lines.append(
        "M5 built one construct per payload and said why: a single "
        "transcript carrying several antibodies needs a 2A peptide or an "
        "internal ribosome entry site, and neither had been sourced to the "
        "standard this repository holds sequences to. The attempt to "
        "source a 2A peptide failed on every endpoint, and the decision "
        "that killed the route became a rule: no sequence enters this "
        "repository unless it has been located in the database that "
        "publishes it and recorded with a checksum.",
    )
    lines.append("")
    lines.append(
        "The rule held, and the route is open now. The porcine "
        "teschovirus-1 2A peptide appears, identical over all 19 residues, "
        "in the expression constructs of two deposited structures from two "
        "different groups: PDB entry 31FW, entity 4, residues 51 to 69, "
        "and PDB entry 8E99, entity 2, residues 4 to 22, in both cases "
        "preceded by the glycine-serine-glycine spacer the "
        "characterisation recommends. Both entries are committed under "
        "`data/reference/raw/` with checksums, and the test suite extracts "
        "the peptide from each and checks the two agree with the engine's "
        "constant. The characterisation itself is Kim and colleagues, 2011.",
    )
    lines.append("")

    lines.append("## 2. What the skip costs, and which payload pays")
    lines.append("")
    lines.append(
        "A 2A element is not a proteolytic cleavage site. The ribosome "
        "simply does not make the peptide bond between the last glycine of "
        "the element and the proline that follows it: translation pauses, "
        "the upstream product is released with the element still attached, "
        "and the downstream product begins with that proline. Nothing is "
        "cut, so nothing is clean.",
    )
    lines.append("")
    lines.append(
        f"The upstream product carries a {len(GSG_P2A_LINKER)}-residue "
        "C-terminal scar: the spacer and the 2A peptide up to the skip "
        "glycine. The published binding data for TPL1158_01_C09 belongs to "
        f"the {len(tpl.sequence)}-residue nanobody, not to the "
        f"{len(tpl.sequence) + len(GSG_P2A_LINKER)}-residue scarred "
        "version this transcript encodes, and the report says so rather "
        "than quietly assuming the scar is benign.",
    )
    lines.append("")
    lines.append(
        "The downstream product starts with one proline, and when the "
        "downstream payload is secreted that proline sits in front of a "
        "signal peptide, which signal peptidase removes together with the "
        "proline. The mature downstream product is therefore the "
        "published molecule exactly.",
    )
    lines.append("")
    lines.append(
        "That asymmetry is the design rule, and it decides the order: the "
        "payload whose published data matters more goes downstream. E113 "
        "has reached L4 on the validation ladder, in vivo neutralisation, "
        "against L2, binding only, for the nanobody, so E113 goes "
        "downstream and the scar goes on the nanobody whose data was "
        "binding-only to begin with.",
    )
    lines.append("")

    lines.append("## 3. The construct")
    lines.append("")
    lines.append(f"The construct is {construct.length} nt:")
    lines.append("")
    lines.append(construct.element_table())
    lines.append("")
    lines.append("| Measurement | Value |")
    lines.append("| --- | --- |")
    lines.append(f"| Coding region length | {len(design.cds)} nt |")
    lines.append(f"| GC content of the coding region | {design.gc_total:.4f} |")
    lines.append(f"| CAI against human codon usage | {design.cai:.4f} |")
    lines.append(
        f"| GC range across 30-nucleotide blocks | {min(design.block_gc):.3f} "
        f"to {max(design.block_gc):.3f} |",
    )
    lines.append(
        f"| Forbidden motifs in the coding elements | {len(coding_hits)} |",
    )
    lines.append(
        f"| Forbidden motifs in the UTR elements | {len(utr_hits)} |",
    )
    lines.append(
        f"| CAI is optimal subject to the GC band | "
        f"{'yes' if design.cai_is_optimal else 'no'} |",
    )
    lines.append("")

    lines.append("## 4. What the engine checked")
    lines.append("")
    lines.append(
        "- **The P2A peptide is the same molecule in two deposited "
        "constructs.** The engine reads both raw files and checks each "
        "against the constant in `src/mrna/p2a.py`; the test suite does "
        "the same against the committed raw downloads.",
    )
    lines.append(
        "- **The open reading frame translates to the designed protein.** "
        f"{len(design.protein)} residues, checked by translation, not by "
        "counting codons.",
    )
    lines.append(
        "- **The mature upstream product is the published nanobody plus "
        f"its scar, {len(products['upstream_mature'])} residues**, checked "
        "by translation of the corresponding element slices.",
    )
    lines.append(
        "- **The mature downstream product is the published E113 sequence, "
        f"{len(products['downstream_mature'])} residues, reproduced "
        "exactly.** This is the check that makes the design claim true "
        "rather than asserted: the carry-over proline and the second "
        "signal peptide are gone, and what remains is the molecule the "
        "2024 paper measured.",
    )
    lines.append(
        f"- **The coding elements carry no forbidden motif.** The policy "
        f"enforces the polyadenylation signal out of every coding element, "
        f"and this transcript carries {len(coding_hits)} occurrences there; "
        f"the {len(utr_hits)} the construct does carry are both in the "
        "UTRs, where the licensed product carries them too, and the M2 "
        "report explains why they are not defects.",
    )
    lines.append("")

    lines.append("## 5. What this transcript is not")
    lines.append("")
    lines.append(
        "- **It is not two published antibodies.** It is one scarred "
        "nanobody and one exact scFv, on one message. The scar is stated "
        "in the element table, in the FASTA header of the upstream "
        "product, and in this paragraph, three times, because a design "
        "whose payload differs from the published molecule has to say so "
        "where a reader will look.",
    )
    lines.append(
        "- **It is not a free cocktail.** Both payloads bind in vitro and "
        "the scFv has neutralisation data on its own; nothing published "
        "measures the two molecules together, on one message or on two. "
        "The M5 rule that a plan is not a molecule applies here with more "
        "force: co-expression is not co-formulation, and the ratio of the "
        "two products depends on skip efficiency, which the "
        "characterisation paper reports as high for P2A but which this "
        "project cannot measure.",
    )
    lines.append(
        "- **It is not licensed.** TPL1158_01_C09 is named 71 times in "
        "WO2025078529A1, published by Danmarks Tekniske Universitet, and "
        "the position is unchanged from M5: an open programme built on "
        "that sequence is not free to exist, and this construct is a "
        "design artifact rather than a freedom-to-operate opinion.",
    )
    lines.append("")

    lines.append("## 6. The third payload, and what it waits for")
    lines.append("")
    lines.append(
        f"The refresh added a third transcribable payload, D2, the "
        f"{len(d2.sequence)}-residue designed nanobody against "
        "alpha-cobratoxin whose sequences sit in an openly downloadable "
        "preprint supplement. It is committed reference data now, with "
        "the extraction pinned to the deposited row by the 9GCN framework "
        "prefix. It does not join this transcript, for two reasons the "
        "report states rather than resolves: a three-payload transcript "
        "has two scars and two ratios where this one has one, and D2 is a "
        "designed sequence whose liabilities the preprint reports but "
        "whose molecule nobody has neutralised in vivo against whole "
        "venom at the level E113 has. The second P2A site a "
        "three-payload version needs is one optimiser run away.",
    )
    lines.append("")

    lines.append("## 7. What comes next")
    lines.append("")
    lines.append(
        "The regional-variant strand: the venom variation inside one "
        "species, mapped for the Indian spectacled cobra by Laxme and "
        "colleagues in 2021, treated as a design input rather than a "
        "footnote. The tracker names the two gaps a variant would fill if "
        "a sequence ever lands: krait phospholipase A2, and the "
        "metalloproteinases of Indian Russell's viper.",
    )
    lines.append("")

    lines.append("## References")
    lines.append("")
    lines.append(
        "- Kim JH, Lee SR, Li LH, et al. 2011. High cleavage efficiency of "
        "a 2A peptide derived from porcine teschovirus-1 in human cell "
        "lines, zebrafish and mice. PLoS ONE 6(4):e18556. "
        "doi:10.1371/journal.pone.0018556.",
    )
    lines.append(
        "- Donnelly ML, Hughes LE, Luke G, et al. 2001. The 2A/2B "
        "polyprotein 'cleavage' mechanism. J Gen Virol 82:1027-1041.",
    )
    lines.append(
        "- Szymczak AL, Workman CJ, Wang Y, et al. 2004. Correction of "
        "multi-gene deficiency in vivo using a single 'self-cleaving' 2A "
        "peptide-based vector. Nat Biotechnol 22:589-594.",
    )
    lines.append(
        "- PDB entries 31FW and 8E99, the two deposited constructs the "
        "peptide is read from.",
    )
    lines.append(
        "- RCSB Protein Data Bank entries 9GCN and 8IA6, the payloads; "
        "Overath et al. 2026, bioRxiv 2026.09.01.748349, the third "
        "payload.",
    )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    """Print or write the M6 report and design files."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the report and the design files",
    )
    args = parser.parse_args()
    report = build_report()
    if args.write:
        construct, design, products = build_construct(
            load_codon_usage(REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv"),
        )
        paths = write_design(construct, design, products)
        REPORT_PATH.write_text(report, encoding="utf-8")
        print(f"wrote {REPORT_PATH.relative_to(ROOT)}")
        for path in paths:
            print(f"wrote {path.relative_to(ROOT)}")
    else:
        print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
