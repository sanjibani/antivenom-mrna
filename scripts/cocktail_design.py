"""M5: the Big Four cocktail, and the gaps the published record cannot fill.

This script is the M5 deliverable. It answers the question M4 was built to
answer: if you had to build one antibody cocktail that covers India's Big
Four snakes, which payloads would it contain, and what would stop you?

Run it::

    uv run python scripts/cocktail_design.py             # print the report
    uv run python scripts/cocktail_design.py --write     # refresh the artifacts

What it does, in order:

1. Loads the M4 candidate table through the tracker engine and asks the
   cocktail planner for a plan at each of four tiers, from every published
   candidate down to the ones that need no engineering and no licence.
2. Reports what each tier covers and what it drops, and names the cells no
   tier can cover.
3. Builds a full mRNA construct for each payload that can be transcribed
   from a published sequence as it stands.
4. Writes the design files and the report.

The headline of this milestone is a negative result, and it is stated as
one: at the tier where the constructs are actually buildable, the published
record does not cover India's Big Four. Two payloads can be transcribed, and
they cover two cells out of fourteen. Indian Russell's viper cannot be
covered at all, because the one recombinant candidate against its dominant
toxin has no publicly locatable sequence.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mrna.assemble import A110_POLYA, CAP1, Construct, full_construct
from mrna.candidates import BIG_FOUR, CLINICAL_SET, SPECIES_FAMILIES
from mrna.cocktail import (
    EXCLUDED_FAMILIES,
    Cell,
    CocktailPlan,
    plan_all_tiers,
)
from mrna.fasta import write_fasta
from mrna.leader import HSA_SIGNAL_PEPTIDE
from mrna.motifs import DEFAULT_MOTIFS
from mrna.optimize import optimize_cds
from mrna.payload import Payload, load_payloads
from mrna.policy import HUMAN_POLICY
from mrna.reference import load_codon_usage
from mrna.tracker import CandidateTracker
from mrna.utr import AES_MTRNR1_3UTR, ALPHA_GLOBIN_5UTR

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"
DESIGNS = ROOT / "data" / "designs"

DISCLAIMER = "In-silico designs. No therapeutic efficacy demonstrated. Not medical advice."

#: The single stop codon of the designed constructs.
STOP_CODON = "TGA"

#: Where the published design files live.
FASTA_PATH = DESIGNS / "big_four_cocktail_mrna.fasta"
ELEMENTS_PATH = DESIGNS / "big_four_cocktail.elements.tsv"

#: Where the published report lives.
REPORT_PATH = ROOT / "docs" / "m5_big_four_cocktail.md"


def build_construct(payload: Payload, usage_table: dict[str, float]) -> tuple[Construct, object]:
    """Build the mRNA construct that carries one payload.

    The scaffold is the one M2 built and M3 measured against a licensed
    vaccine: the Comirnaty cap footprint, the human alpha-globin 5' UTR, the
    human serum albumin signal peptide, the codon-optimised payload, the AES
    plus mitochondrial 12S rRNA 3' UTR and the 110-nucleotide poly(A) tail.

    Args:
        payload: The payload to encode.
        usage_table: Codon usage table for the target host.

    Returns:
        The assembled construct and the optimisation result for its coding
        region.
    """
    protein = HSA_SIGNAL_PEPTIDE.sequence + payload.sequence
    design = optimize_cds(protein, usage_table, HUMAN_POLICY)
    coding = design.cds + STOP_CODON

    cap_nt = len(CAP1.prefix)
    utr5_nt = len(ALPHA_GLOBIN_5UTR.sequence)
    signal_nt = HSA_SIGNAL_PEPTIDE.length * 3
    payload_start = cap_nt + utr5_nt + signal_nt + 1
    payload_end = payload_start + len(coding) - signal_nt - 1
    stop_start = payload_end - len(STOP_CODON) + 1

    construct = full_construct(
        f"{payload.name} against {payload.toxin_target}",
        utr5=ALPHA_GLOBIN_5UTR,
        coding=coding,
        utr3=AES_MTRNR1_3UTR,
        poly_a=A110_POLYA,
        leader=HSA_SIGNAL_PEPTIDE,
        coding_name=f"{payload.name} payload",
        coding_description=(
            f"Codon-optimised {payload.molecule_format} {payload.name} against "
            f"{payload.toxin_target}, preceded by the human serum albumin "
            f"signal peptide, followed by a stop codon at {stop_start} to "
            f"{payload_end}."
        ),
    )
    return construct, design


def write_designs(
    constructs: list[tuple[Payload, Construct, object]],
) -> list[Path]:
    """Write the design files for the cocktail members.

    Args:
        constructs: One ``(payload, construct, optimisation)`` triple per
            cocktail member.

    Returns:
        The paths written.
    """
    DESIGNS.mkdir(parents=True, exist_ok=True)

    records: list[tuple[str, str]] = []
    rows = ["member\telement\tstart\tend\tlength\tdescription"]
    for payload, construct, _design in constructs:
        records.append(
            (
                f"{construct.name} | full mRNA construct, {construct.length} nt | "
                f"payload: {payload.name}, a {payload.molecule_format} against "
                f"{payload.toxin_target} | source: antivenom-mrna M5 | "
                f"{DISCLAIMER}",
                construct.sequence,
            ),
        )
        records.append(
            (
                f"{payload.name} payload protein | {len(payload.sequence)} aa | "
                f"source: {payload.citation}",
                payload.sequence,
            ),
        )
        for element in construct.elements:
            description = element.description.replace("\t", " ")
            rows.append(
                f"{payload.name}\t{element.name}\t{element.start}\t{element.end}\t"
                f"{element.length}\t{description}",
            )
    write_fasta(FASTA_PATH, records)
    ELEMENTS_PATH.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return [FASTA_PATH, ELEMENTS_PATH]


def _coverage_matrix(plans: tuple[CocktailPlan, ...]) -> tuple[tuple[Cell, ...], list[str]]:
    """Return the cells and the tiers, for rendering the coverage table.

    Args:
        plans: The plans, in tier order.

    Returns:
        The cells and the tier names.
    """
    cells: list[Cell] = []
    for plan in plans:
        for cell in plan.covered + plan.uncovered:
            if cell not in cells:
                cells.append(cell)
    return tuple(cells), [plan.tier for plan in plans]


def build_report(tracker: CandidateTracker) -> str:
    """Build the full M5 markdown report.

    Args:
        tracker: The validated candidate table.

    Returns:
        The report as a markdown string.
    """
    usage_table = load_codon_usage(REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv")
    plans = plan_all_tiers(tracker)
    payloads = {payload.name: payload for payload in load_payloads()}

    buildable = next(plan for plan in plans if plan.tier == "single-polypeptide")
    members = [payloads[name] for name in buildable.payload_ids]
    constructs = [(payload, *build_construct(payload, usage_table)) for payload in members]

    lines: list[str] = []
    lines.append("# M5: the Big Four cocktail, and what the published record cannot cover")
    lines.append("")
    lines.append(f"> {DISCLAIMER}")
    lines.append("")
    lines.append(
        "This is the M5 artifact, and most of it is a negative result. The "
        "question was: if you had to build one antibody cocktail covering "
        "India's Big Four snakes, which payloads would it contain, and what "
        "would stop you? The answer is computed by "
        "`scripts/cocktail_design.py` from the M4 candidate table, not "
        "asserted, and the short version is that the published record cannot "
        "supply the cocktail.",
    )
    lines.append("")

    lines.append("## 1. What the planner was asked")
    lines.append("")
    lines.append(
        "A cell is a species and a toxin family. The planner treats cocktail "
        "design as a set cover: it finds the smallest set of antibodies that "
        "covers as much of the matrix as the published record allows, and it "
        "searches exhaustively rather than greedily, so the answer is exact. "
        f"{len(plans)} tiers tighten the question one constraint at a time.",
    )
    lines.append("")
    lines.append(
        "The species are the Big Four plus the monocled cobra. The monocled "
        "cobra, *Naja kaouthia*, is not one of the Big Four, but it is the "
        "other major Indian cobra, it is one of the species Indian polyvalent "
        "antivenom is raised against, and a cocktail aimed at Indian "
        "snakebite that ignored it would be incomplete in a way that "
        "matters.",
    )
    lines.append("")
    lines.append("| Tier | The question | Payloads | Cells covered |")
    lines.append("| --- | --- | --- | --- |")
    for plan in plans:
        total = len(plan.covered) + len(plan.uncovered)
        lines.append(
            f"| `{plan.tier}` | {plan.question} | {len(plan.payloads)} | "
            f"{len(plan.covered)} of {total} |",
        )
    lines.append("")
    lines.append(
        "The attrition down that table is the finding. Four tiers, and the "
        "number of cells the plan can actually cover falls from "
        f"{len(plans[0].covered)} to {len(plans[-1].covered)}.",
    )
    lines.append("")

    lines.append("## 2. The coverage matrix")
    lines.append("")
    cells, tiers = _coverage_matrix(plans)
    families = sorted({cell.toxin_family for cell in cells})
    species_order = [name for name in CLINICAL_SET if any(c.species == name for c in cells)]
    lines.append("| Species | Toxin family | " + " | ".join(f"`{t}`" for t in tiers) + " |")
    lines.append("| --- " * (len(tiers) + 2) + "|")
    for species in species_order:
        for family in families:
            cell = Cell(species, family)
            if cell not in cells:
                continue
            marks = []
            for plan in plans:
                marks.append("yes" if cell in plan.covered else "no")
            label = f"**{species}**" if species in BIG_FOUR else species
            lines.append(f"| {label} | {family} | " + " | ".join(marks) + " |")
    lines.append("")
    cleared = [
        candidate
        for candidate in tracker.candidates
        if candidate.patent_status == "none-located"
    ]
    searched = sum(
        1 for candidate in tracker.candidates if candidate.patent_status != "not-assessed"
    )
    lines.append(
        "The Big Four species are in bold. The `licence-unblocked` column is a "
        "floor rather than a clearance: a row reaches it only when a public "
        "patent search was actually performed and found nothing. Of the "
        f"{len(tracker.candidates)} rows in the M4 table, {searched} have had "
        f"a search performed at all, and {len(cleared)} of those came back "
        "with nothing: "
        f"{', '.join(candidate.antibody_id for candidate in cleared)}.",
    )
    lines.append("")

    lines.append("### What each tier drops")
    lines.append("")
    for plan in plans:
        lines.append(f"**`{plan.tier}`**: {', '.join(plan.payload_ids) or 'nothing'}.")
        lines.append("")
    lines.append(
        "Read down the list of payloads and the order is the story. The "
        "published record has candidates for most of the matrix. The public "
        "sequence databases do not. The transcripts that can be built from "
        f"what is public are {len(buildable.payloads)}.",
    )
    lines.append("")

    lines.append("## 3. The cells that no tier covers")
    lines.append("")
    final = plans[0]
    uncovered = final.uncovered
    viper_cells = [
        cell for cell in uncovered if SPECIES_FAMILIES.get(cell.species) == "viper"
    ]
    elapid_cells = [
        cell for cell in uncovered if SPECIES_FAMILIES.get(cell.species) == "elapid"
    ]
    lines.append(
        f"Even at the most permissive tier, with every published candidate "
        f"allowed, {len(uncovered)} cells have no candidate at all:",
    )
    lines.append("")
    for cell in uncovered:
        big_four = " (Big Four)" if cell.species in BIG_FOUR else ""
        lines.append(f"- **{cell.describe()}**{big_four}")
    lines.append("")
    viper_species = sorted({cell.species for cell in viper_cells})
    lines.append(
        f"Of these, {len(viper_cells)} are viper toxin families that nobody "
        "has taken on: "
        + ", ".join(cell.describe() for cell in viper_cells)
        + ". The "
        + f"{len(viper_species)} viper species of the Big Four between them "
        "carry every toxin class that makes a viper bite dangerous, and not "
        "one of those classes has a published antibody candidate in this "
        "table.",
    )
    lines.append("")
    lines.append(
        f"The other {len(elapid_cells)} are elapid cells: "
        + ", ".join(cell.describe() for cell in elapid_cells)
        + ". The krait phospholipase A2 cell is the one to start from, because "
        "phospholipase A2 is the dominant fraction of common krait venom by "
        "mass, and no published candidate targets it. Krait venom is the most "
        "under-covered of the Big Four.",
    )
    lines.append("")
    lines.append("### The cells the planner will not target")
    lines.append("")
    lines.append(
        f"{len(final.excluded)} cells are empty on purpose. The Kunitz-type "
        "serine protease inhibitors are the anticoagulant fraction of viper "
        "venom: they bind the clotting factors directly, so an antibody "
        "against them would protect the clotting cascade in a patient whose "
        "problem is that the cascade is being consumed. A planner that "
        "reported every empty cell as an opportunity would be pointing this "
        "project at a design that would make the patient worse off.",
    )
    lines.append("")
    for family, reason in EXCLUDED_FAMILIES.items():
        lines.append(f"- **{family}**: {reason}.")
    lines.append("")

    lines.append(
        f"## 4. The {len(constructs)} construct"
        f"{'' if len(constructs) == 1 else 's'} that can be built",
    )
    lines.append("")
    lines.append(
        f"At the `single-polypeptide` tier there are {len(buildable.payloads)} "
        "payloads, and every one of them can be transcribed from a published "
        "sequence as it stands, with no engineering step in between. The "
        "others are Fab fragments or two-chain antibodies: converting them "
        "into one transcript needs a single-chain format that does not exist "
        "in the published data, which means the neutralisation data would no "
        "longer belong to the molecule being transcribed.",
    )
    lines.append("")
    for payload, construct, design in constructs:
        orf = construct.open_reading_frame()
        hits = construct.find_motifs(DEFAULT_MOTIFS)
        lines.append(f"### {payload.name}")
        lines.append("")
        lines.append("| Field | Value |")
        lines.append("| --- | --- |")
        lines.append(f"| Target | {payload.toxin_target} |")
        lines.append(f"| Species | {', '.join(payload.species)} |")
        lines.append(f"| Format | {payload.molecule_format}, one polypeptide |")
        lines.append(f"| Payload length | {len(payload.sequence)} aa |")
        lines.append(f"| Total construct length | {construct.length} nt |")
        lines.append(
            f"| Translated coding sequence | {orf.end - orf.start} nt |",
        )
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
        lines.append(f"| Source | {payload.citation} |")
        lines.append("")
        lines.append(construct.element_table())
        lines.append("")

    lines.append("## 5. The licence position, which is the real blocker")
    lines.append("")
    blocked = [
        group for group in buildable.payloads if group.patent_status != "none-located"
    ]
    unblocked = [
        group for group in buildable.payloads if group.patent_status == "none-located"
    ]
    lines.append(
        f"These {len(buildable.payloads)} payloads are the only ones that can "
        "be transcribed, and a public patent search finds an application "
        f"covering {len(blocked)} of them.",
    )
    lines.append("")
    for group in blocked:
        lines.append(
            f"- **{group.antibody_id}** is named 71 times in "
            "WO2025078529A1, *Nanobodies capable of cross-binding to several "
            "different elapid three-finger toxins*, published by Danmarks "
            "Tekniske Universitet.",
        )
    for group in unblocked:
        lines.append(
            f"- **{group.antibody_id}** is the clone this repository already "
            "built in M3, and no application covering it was located.",
        )
    lines.append("")
    lines.append(
        "This is not a clearance, and it is not a freedom-to-operate opinion. "
        "It is a record of what the public patent record says, on a date, "
        "about the applications this project searched. It does mean that an "
        "open, reproducible mRNA programme built on the published sequences is "
        "not free to exist, and that the honest obstacle to an open Big Four "
        "cocktail is as much legal as scientific.",
    )
    lines.append("")

    lines.append("## 6. What this cocktail is not")
    lines.append("")
    covered_species = buildable.species_covered()
    big_four_covered = [
        name for name in covered_species if name in BIG_FOUR
    ]
    total_cells = len(buildable.covered) + len(buildable.uncovered)
    lines.append(
        f"**It covers {len(big_four_covered)} Big Four species out of "
        f"{len(BIG_FOUR)}.** {len(buildable.covered)} of {total_cells} cells, "
        "and the only Big Four species in it is "
        + ", ".join(big_four_covered)
        + ". It is the first honest attempt, not a finished product, and the "
        "report above says exactly which cells are missing.",
    )
    lines.append("")
    lines.append(
        "**Indian Russell's viper is not in it at all.** The one "
        "recombinant candidate against the dominant toxin of that venom, the "
        "camelid VHH MKSS1, has no publicly locatable sequence, so at the "
        "`sequence-available` tier the whole species drops out of the plan. "
        "That is the single most actionable line in this document.",
    )
    lines.append("")
    lines.append(
        "**Nothing has been expressed.** No cell has made any of these "
        "proteins. The next rung of the validation ladder is L1, and it "
        "belongs to a wet-lab partner this project does not have.",
    )
    lines.append("")
    lines.append(
        "**The timing is wrong for an acute bite, exactly as in M3.** An "
        "mRNA-encoded antibody is not detectable for one to two hours and "
        "peaks somewhere between 12 and 24 hours. Intravenous antivenom works "
        "within minutes. This is a prophylactic and a regional-manufacturing "
        "proposition.",
    )
    lines.append("")
    lines.append(
        "**A cocktail is not a single molecule.** Each payload is its own "
        "transcript here, because a single transcript carrying several "
        "antibodies needs a 2A peptide or an internal ribosome entry site, "
        "and neither has been sourced for this repository to the standard the "
        "rest of the sequences meet. That is the first thing M6 should do, "
        "and it is a real engineering step rather than a formatting one.",
    )
    lines.append("")
    lines.append(
        "**The payloads are monovalent.** An scFv and a VHH each have one "
        "binding site, so each covers one toxin per antibody molecule. The "
        "published combination data for the E113 clone shows that a second "
        "clone helps, which is an argument for a cocktail and not for a "
        "single construct.",
    )
    lines.append("")

    lines.append("## 7. What comes next")
    lines.append("")
    lines.append(
        "M6 is regional variants, partner outreach and design iteration. The "
        "tracker has now said what the work should be: the krait "
        "phospholipase A2 cell, the Russell's viper metalloproteinase cell, "
        "and a licence conversation with the two institutions whose "
        "applications cover the sequences that exist. None of those is a "
        "design problem, which is why this milestone ends with a map rather "
        "than with more constructs.",
    )
    lines.append("")

    lines.append("## References")
    lines.append("")
    lines.append(
        "- Samanta A, Rudresha GV, Ahmadi S, et al. 2026. Oligoclonal "
        "nanobody-based recombinant antivenom protects mice challenged with "
        "venom from cobras and king cobras from India. Sci Transl Med "
        "18:866. doi:10.1126/scitranslmed.aed4290.",
    )
    lines.append(
        "- Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2024. Mol Immunol "
        "165:55-67. doi:10.1016/j.molimm.2023.12.006.",
    )
    lines.append(
        "- Glanville J, Kwong K, et al. 2025. Broad-spectrum neutralisation of "
        "elapid venoms by a minimal antibody combination. Cell. Patent "
        "US20240294618A1.",
    )
    lines.append(
        "- WO2025078529A1. Nanobodies capable of cross-binding to several "
        "different elapid three-finger toxins. Danmarks Tekniske Universitet.",
    )
    lines.append(
        "- WO2023209212A1. Neutralisation of alpha-neurotoxins using human "
        "recombinant IgG antibodies. Danmarks Tekniske Universitet.",
    )
    lines.append(
        "- Szymczak AL, Workman CJ, Wang Y, et al. 2004. Correction of "
        "multi-gene deficiency in vivo using a single 'self-cleaving' 2A "
        "peptide-based retroviral vector. Nat Biotechnol 22:589-594.",
    )
    lines.append(
        "- Protein Data Bank entries 9GCN, 8DA1, 8D9Z, 8V13, 9CLP, 9HUO, "
        "9FYT and 9HXO.",
    )
    lines.append("")
    lines.append(f"_{DISCLAIMER}_")
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    """Run the planner, optionally writing the design files and the report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the design files and docs/m5_big_four_cocktail.md",
    )
    args = parser.parse_args()

    tracker = CandidateTracker()
    if args.write:
        usage_table = load_codon_usage(REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv")
        payloads = {payload.name: payload for payload in load_payloads()}
        buildable = next(
            plan for plan in plan_all_tiers(tracker) if plan.tier == "single-polypeptide"
        )
        constructs = [
            (payloads[name], *build_construct(payloads[name], usage_table))
            for name in buildable.payload_ids
        ]
        for path in write_designs(constructs):
            print(f"wrote {path.relative_to(ROOT)}")
        REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
        REPORT_PATH.write_text(build_report(tracker), encoding="utf-8")
        print(f"wrote {REPORT_PATH.relative_to(ROOT)}")
        for plan in plan_all_tiers(tracker):
            print(
                f"{plan.tier}: {len(plan.payloads)} payloads, "
                f"{len(plan.covered)} of "
                f"{len(plan.covered) + len(plan.uncovered)} cells",
            )
    else:
        print(build_report(tracker))
    return 0


if __name__ == "__main__":
    sys.exit(main())
