"""M4: the candidate tracker for India's Big Four snakes.

This script is the M4 deliverable. It maps every published antibody candidate
that has a claimed neutralising activity against a toxin, or against the whole
venom, of one of India's Big Four snakes, and it works out which cells of the
map are empty.

Run it::

    uv run python scripts/candidate_tracker.py             # print the report
    uv run python scripts/candidate_tracker.py --write     # refresh the artifacts

What it does, in order:

1. Loads the curated table from ``src/mrna/candidates.py`` through the
   tracker engine, which validates every row against its controlled
   vocabularies before anything else happens.
2. Works out which toxin families each candidate covers, and at what rung of
   the validation ladder.
3. Finds the species and toxin family cells that no candidate covers.
4. Writes ``data/candidates.tsv`` and ``docs/m4_candidate_tracker.md``.

The tracker is deliberately not a bibliography. A literature map that cannot be
checked is a literature map with opinions in it, so the engine recomputes the
mRNA suitability column from the format and the sequence source, and refuses
to load a row where the two disagree.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from mrna.tracker import CandidateTracker

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
DOCS = ROOT / "docs"

DISCLAIMER = "In-silico designs. No therapeutic efficacy demonstrated. Not medical advice."

#: Where the published table lives.
TSV_PATH = DATA / "candidates.tsv"

#: Where the published report lives.
REPORT_PATH = DOCS / "m4_candidate_tracker.md"


def build_report(tracker: CandidateTracker) -> str:
    """Build the M4 markdown report.

    Args:
        tracker: The validated tracker.

    Returns:
        The report as a markdown string.
    """
    return tracker.render_markdown()


def write_artifacts(tracker: CandidateTracker) -> list[Path]:
    """Write the published table and the published report.

    Args:
        tracker: The validated tracker.

    Returns:
        The paths written.
    """
    DATA.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    TSV_PATH.write_text(tracker.to_tsv(), encoding="utf-8")
    REPORT_PATH.write_text(build_report(tracker), encoding="utf-8")
    return [TSV_PATH, REPORT_PATH]


def main() -> int:
    """Run the tracker, optionally writing the artifacts."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write",
        action="store_true",
        help="write data/candidates.tsv and docs/m4_candidate_tracker.md",
    )
    args = parser.parse_args()

    tracker = CandidateTracker()
    if args.write:
        for path in write_artifacts(tracker):
            print(f"wrote {path.relative_to(ROOT)}")
        gaps = tracker.gaps()
        print(f"{len(tracker.candidates)} candidates, {len(gaps)} uncovered cells")
    else:
        print(build_report(tracker))
    return 0


if __name__ == "__main__":
    sys.exit(main())
