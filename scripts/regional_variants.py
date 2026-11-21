"""M7: regional venom variation as a design input.

The species-level coverage matrix treats a species' venom as one thing.
The two Sunagar-group surveys of 2021 print regional composition data for
the spectacled cobra and for Russell's viper that say otherwise, and the
commercial antivenom's regional failures are the clinical consequence.
This script turns the printed fractions into the question a matrix cannot
ask: for each surveyed population, how much of the venom's mass falls in
cells the cocktail plan can and cannot cover, and which of the remainder
the matrix does not even model.

Run it::

    uv run python scripts/regional_variants.py             # print the report
    uv run python scripts/regional_variants.py --write     # refresh the report
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mrna.cocktail import plan_cocktail
from mrna.regional import (
    CATEGORY_TO_FAMILY,
    REGIONAL_VENOM_CITATIONS,
    regional_coverage,
)
from mrna.tracker import CandidateTracker

ROOT = Path(__file__).resolve().parent.parent
REPORT_PATH = ROOT / "docs" / "m7_regional_variants.md"

DISCLAIMER = (
    "In-silico designs. No therapeutic efficacy demonstrated. "
    "Not medical advice."
)

#: The tier names, in the order the planner asks them. A family covered at
#: a strict tier is covered at every looser one; the report records the
#: strictest tier at which each mapped family is still covered.
TIER_ORDER = ("published", "sequence-available", "single-polypeptide", "licence-unblocked")


def build_report() -> str:
    """Build the full M7 markdown report.

    Returns:
        The report as a markdown string.
    """
    tracker = CandidateTracker()
    records = regional_coverage(tracker)
    buildable = plan_cocktail(tracker, "single-polypeptide")
    buildable_ids = [group.antibody_id for group in buildable.payloads]

    lines: list[str] = []
    lines.append("# M7: regional venom variation as a design input")
    lines.append("")
    lines.append(f"> {DISCLAIMER}")
    lines.append("")

    lines.append("## 1. Why the species-level matrix is not enough")
    lines.append("")
    lines.append(
        "The coverage matrix of M4 and the planner of M5 ask their "
        "questions per species. A species is not one venom. The spectacled "
        "cobra of Punjab carries about 80 per cent neurotoxic three-finger "
        "toxins by venom mass; the desert population of Rajasthan spreads "
        "its venom across cytotoxins, about 42 per cent, and "
        "phospholipase A2, about 20 per cent, and the commercial "
        "antivenom fails against it completely at five times the lethal "
        "dose. The Russell's viper of the Gangetic Plain carries about 39 "
        "per cent Kunitz-type inhibitors, the Western Ghats population 59 "
        "per cent phospholipase A2. The same species, different venoms, "
        "and a cocktail sized for the average is sized for neither.",
    )
    lines.append("")
    lines.append(
        "This report does one thing with the printed regional data: it "
        "weights the coverage matrix by venom mass, per population. The "
        "species-level question is which cells a plan covers. The regional "
        "question is which cells a plan covers *in the venom a patient in "
        "a given state is actually carrying*.",
    )
    lines.append("")

    lines.append("## 2. What the surveys measured, and what they did not")
    lines.append("")
    lines.append(
        "Both surveys are open access and both print their mass "
        "spectrometry for three populations per species, chosen for their "
        "distinct chromatographic and toxicity profiles. The other "
        "populations of both surveys were characterised by chromatography "
        "alone, and no composition is recorded for them here, because a "
        "fraction nobody printed is not a fraction anybody can check.",
    )
    lines.append("")
    lines.append(
        "The surveys' toxin categories are not the matrix's toxin "
        "families, and the correspondence this report uses is an "
        "approximation it names. The cobra survey's neurotoxic "
        "three-finger toxins include both long-chain and short-chain "
        "alpha-neurotoxins; the matrix separates 3FTx-L, the long-chain "
        "neurotoxins, from 3FTx-S, the cytotoxins. This report maps the "
        "survey's neurotoxic fraction onto 3FTx-L and its cytotoxic "
        "fraction onto 3FTx-S. A survey category with no matrix family at "
        "all, such as the snake venom serine proteases and C-type lectins "
        "of the viper survey, is recorded as unmapped rather than squeezed "
        "into the nearest family. The unmapped fraction is itself a "
        "finding: it is the part of the venom the matrix does not model.",
    )
    lines.append("")

    lines.append("## 3. The populations, weighted tier by tier")
    lines.append("")
    lines.append(
        "Every mapped family below is checked against the plan at all four "
        "tiers, and the table records the strictest tier at which the "
        "family is still covered. The tiers are the planner's: published, "
        "sequence-available, single-polypeptide, licence-unblocked. The "
        "single-polypeptide tier is the one where a transcript can carry "
        f"the payload, and its plan is {', '.join(buildable_ids)}. A family "
        "the planner excludes on purpose, the anticoagulant Kunitz-type "
        "fraction, is recorded as excluded, not as uncovered: the planner "
        "refuses it because an antibody against it would make the "
        "coagulopathy worse, which is the M5 rule and it applies region by "
        "region.",
    )
    lines.append("")

    for record in records:
        profile = record["profile"]
        lines.append(
            f"### {profile.population} ({profile.state}, {profile.zone}) "
            f"- {profile.species}",
        )
        lines.append("")
        lines.append("| Survey category | Venom mass | Matrix family | Status |")
        lines.append("| --- | --- | --- | --- |")
        mapping_status = {
            family_record["family"]: family_record
            for family_record in record["families"]
        }
        for category, fraction in profile.fractions.items():
            family = CATEGORY_TO_FAMILY[profile.species].get(category)
            if family is None:
                lines.append(
                    f"| {category} | {fraction:.2f} | none | "
                    "unmapped: the matrix does not model this family |",
                )
                continue
            family_record = mapping_status[family]
            if family_record["excluded"]:
                status = "excluded by design"
            elif family_record["covered_through"] is None:
                status = "uncovered at every tier"
            else:
                status = f"covered through `{family_record['covered_through']}`"
            lines.append(
                f"| {category} | {fraction:.2f} | {family} | {status} |",
            )
        if profile.notes:
            lines.append("")
            lines.append(f"Note: {profile.notes}")
        lines.append("")

    lines.append("## 4. What the weighting changes")
    lines.append("")
    lines.append(
        "At the species level, the spectacled cobra's three cells all have "
        "candidates at the published tier, so the species reads as covered. "
        "Weighted by population and read down the tiers, the picture is "
        "different, and the difference is the point:",
    )
    lines.append("")
    lines.append(
        "- **Every cobra cell is covered only by candidates that cannot be "
        "transcribed as they stand, so the whole species goes uncovered at "
        "the tier where a transcript can carry it.** The two neurotoxin "
        "cells are held by Fabs, which have sequences but two-chain "
        "formats; the phospholipase A2 cell is held by an scFv whose "
        "sequence was never located. In Punjab, about 80 per cent of the "
        "venom mass sits behind the 3FTx-L cell when it drops. In "
        "Rajasthan it is about 42 per cent behind the cytotoxin cell and "
        "20 per cent behind phospholipase A2, and the antivenom fails "
        "against this venom completely at five times the lethal dose. The "
        "gap is not a row in a matrix; it is the venom a patient carries.",
    )
    lines.append(
        "- **For the viper, the biggest mapped fraction in two of three "
        "populations is the family the planner excludes on purpose.** The "
        "Gangetic Plain venom is about 39 per cent Kunitz-type inhibitors, "
        "and no antibody against them belongs in a design, because they "
        "are the anticoagulant fraction of a venom whose patient is "
        "already bleeding. The biggest actionable mapped fraction in every "
        "viper population is phospholipase A2, covered only at the "
        "published tier, because its one recombinant candidate has no "
        "publicly locatable sequence: the same data-availability problem "
        "M4 and M5 named, now with a mass attached, 59 per cent of the "
        "Western Ghats venom, 45 per cent of the Deccan Plateau's, 35 per "
        "cent of the Gangetic Plain's.",
    )
    lines.append(
        "- **The viper matrix under-models the venom.** Snake venom serine "
        "proteases are 8 to 27 per cent of the mapped populations' venom "
        "by mass and C-type lectins up to 20 per cent, and neither family "
        "has a column in the gap matrix. Adding them is a decision for the "
        "next milestone, and this report records the need rather than "
        "silently widening the matrix underneath M4's published artifacts.",
    )
    lines.append("")

    lines.append("## 5. What this report is not")
    lines.append("")
    lines.append(
        "- **It is not a regional cocktail design.** Six populations, "
        "three per species, with printed composition, is the survey's "
        "record, not India's map. The bioclimatic modelling of Sarangi and "
        "colleagues in 2025 extends the coverage; the composition data "
        "behind this report is deliberately limited to what was printed.",
    )
    lines.append(
        "- **It is not an antivenom evaluation.** The regional failures "
        "the surveys report, the below-claim potency in Punjab, the "
        "complete failure against the Rajasthan cobra at five times the "
        "lethal dose, the unneutralised clotting against the Gangetic "
        "Plain and coastal vipers, are the surveys' findings and this "
        "report's motivation. An mRNA construct is not the answer to an "
        "acute bite regardless of region, and the README's caveat "
        "applies here with full force.",
    )
    lines.append("")

    lines.append("## References")
    lines.append("")
    for key in sorted(REGIONAL_VENOM_CITATIONS):
        lines.append(f"- {REGIONAL_VENOM_CITATIONS[key]}")
    lines.append(
        "- Sarangi MK, et al. 2025. Predictive modelling of bioclimatic "
        "venom variation in Russell's viper. PLoS Negl Trop Dis. "
        "doi:10.1371/journal.pntd.0012949.",
    )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    """Print or write the M7 report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the report",
    )
    args = parser.parse_args()
    report = build_report()
    if args.write:
        REPORT_PATH.write_text(report, encoding="utf-8")
        print(f"wrote {REPORT_PATH.relative_to(ROOT)}")
    else:
        print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
