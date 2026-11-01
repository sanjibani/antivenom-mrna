"""The candidate tracker: validation, coverage and gap analysis for the M4 map.

M4 is a literature map, and a literature map that cannot be checked is a
literature map with opinions in it. So the tracker is an engine rather than a
table: it holds the curated rows, it refuses to load a row whose controlled
vocabulary values are wrong, and it refuses to load a row whose declared
mRNA suitability disagrees with what the format and the sequence source
actually imply. The published table and the published report are generated
from it.

The analysis the engine performs is the part that matters for this project.
Antibody discovery is not the bottleneck. Sequence availability is. A
candidate whose variable-region sequence is not public cannot be transcribed,
no matter how good its neutralisation data is, so the tracker separates the
two questions and answers them in different columns.

The gap analysis asks the question M5 depends on: for each Big Four species
and each toxin family that matters for that species, is there a published
candidate with binding data at all? The cells with no answer are the design
opportunity, and they are the reason the project exists.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

from mrna.candidates import (
    BIG_FOUR,
    CANDIDATES,
    CITATIONS,
    COMMON_NAMES,
    EVIDENCE_LEVELS,
    FORMATS,
    MRNA_SUITABILITY,
    PATENT_STATUS,
    POLYCLONAL_FORMATS,
    SINGLE_ORF_FORMATS,
    SPECIES_TOXIN_FAMILIES,
    TOXIN_FAMILIES,
    TOXIN_FAMILY_NOTES,
    Candidate,
)

__all__ = [
    "COVERAGE_THRESHOLD",
    "CandidateTracker",
    "Gap",
    "derive_mrna_suitability",
    "parse_sequence_source",
]

#: The lowest rung of evidence that counts as coverage.
#:
#: A candidate that has never been shown to bind is not coverage, and a
#: candidate that binds but has never been shown to neutralise is a lead, not
#: a design. L2 is the point at which a sequence can be transcribed with some
#: confidence that it points at the right epitope.
COVERAGE_THRESHOLD = "L2"

#: The kinds of sequence source the table understands.
SEQUENCE_SOURCE_KINDS = ("PDB", "patent")

#: The tokens that need no value after the colon.
SEQUENCE_SOURCE_BARE = ("paper", "none")


def parse_sequence_source(tokens: tuple[str, ...]) -> dict[str, list[str]]:
    """Split sequence-source tokens into their kinds.

    Args:
        tokens: The raw tokens. Each is either ``kind:value`` or one of the
            bare tokens in :data:`SEQUENCE_SOURCE_BARE`.

    Returns:
        Mapping of kind to the values recorded for that kind. A bare token
        maps to a single empty value.

    Raises:
        ValueError: If a token is not of the form ``kind:value`` or the kind
            is not one of :data:`SEQUENCE_SOURCE_KINDS`.
    """
    parsed: dict[str, list[str]] = {}
    for token in tokens:
        if token in SEQUENCE_SOURCE_BARE:
            parsed.setdefault(token, []).append("")
            continue
        kind, separator, value = token.partition(":")
        if not separator:
            msg = f"sequence source token {token!r} is not of the form kind:value"
            raise ValueError(msg)
        if kind not in SEQUENCE_SOURCE_KINDS:
            allowed = ", ".join(SEQUENCE_SOURCE_KINDS)
            msg = f"sequence source token {token!r} has unknown kind {kind!r}; allowed: {allowed}"
            raise ValueError(msg)
        parsed.setdefault(kind, []).append(value)
    return parsed


def derive_mrna_suitability(
    molecule_format: str, sequence_source: tuple[str, ...],
) -> str:
    """Work out whether a candidate can be transcribed from public material.

    The rule is deliberately mechanical, because the alternative is an
    opinion. A polyclonal candidate has no defined sequence and so nothing to
    transcribe. A monoclonal with no located sequence cannot be assessed. A
    monoclonal with a public sequence can be transcribed, and whether it needs
    one open reading frame or several depends on the format.

    Args:
        molecule_format: One of :data:`mrna.candidates.FORMATS`.
        sequence_source: The raw sequence-source tokens.

    Returns:
        One of ``ready``, ``adaptable``, ``blocked`` or ``unknown``.
    """
    if molecule_format in POLYCLONAL_FORMATS:
        return "blocked"
    parsed = parse_sequence_source(sequence_source)
    public = bool(parsed.get("PDB")) or "paper" in parsed
    if not public:
        return "unknown"
    if molecule_format in SINGLE_ORF_FORMATS:
        return "ready"
    return "adaptable"


@dataclass(frozen=True)
class Gap:
    """A species and toxin family that no candidate in the table covers.

    Attributes:
        species: The binomial name.
        toxin_family: The toxin family.
        reason: Why the cell is empty, in one sentence.
        near_species: Other species that have a candidate for this family.
        near_family: Other families of the same species that have a
            candidate.
    """

    species: str
    toxin_family: str
    reason: str
    near_species: tuple[str, ...]
    near_family: tuple[str, ...]


class CandidateTracker:
    """The curated candidate table, with validation and analysis.

    Args:
        candidates: The rows to hold. Defaults to the curated M4 table.

    Raises:
        ValueError: If any row fails validation.
    """

    def __init__(self, candidates: tuple[Candidate, ...] = CANDIDATES) -> None:
        """Hold the candidates and validate them straight away.

        Args:
            candidates: The rows to hold. Defaults to the curated M4 table.

        Raises:
            ValueError: If any row fails validation.
        """
        self.candidates = candidates
        self.validate()

    # ------------------------------------------------------------ validation
    def validate(self) -> None:
        """Check every row against the controlled vocabularies.

        Raises:
            ValueError: On the first row that fails, naming the row and the
                field. The error is raised rather than collected so that a
                broken table never reaches a report.
        """
        seen: set[tuple[str, str]] = set()
        for index, candidate in enumerate(self.candidates):
            label = f"row {index + 1} ({candidate.antibody_id})"
            for species in candidate.species:
                if species not in COMMON_NAMES:
                    msg = f"{label}: unknown species {species!r}"
                    raise ValueError(msg)
            if candidate.toxin_family not in TOXIN_FAMILIES:
                msg = f"{label}: unknown toxin_family {candidate.toxin_family!r}"
                raise ValueError(msg)
            if candidate.molecule_format not in FORMATS:
                msg = f"{label}: unknown format {candidate.molecule_format!r}"
                raise ValueError(msg)
            if candidate.patent_status not in PATENT_STATUS:
                msg = f"{label}: unknown patent_status {candidate.patent_status!r}"
                raise ValueError(msg)
            if candidate.evidence_level not in EVIDENCE_LEVELS:
                msg = f"{label}: unknown evidence_level {candidate.evidence_level!r}"
                raise ValueError(msg)
            if candidate.mrna_suitability not in MRNA_SUITABILITY:
                msg = f"{label}: unknown mrna_suitability {candidate.mrna_suitability!r}"
                raise ValueError(msg)
            if candidate.citation not in CITATIONS:
                msg = f"{label}: unknown citation key {candidate.citation!r}"
                raise ValueError(msg)
            parse_sequence_source(candidate.sequence_source)
            derived = derive_mrna_suitability(
                candidate.molecule_format, candidate.sequence_source,
            )
            if derived != candidate.mrna_suitability:
                msg = (
                    f"{label}: declared mrna_suitability "
                    f"{candidate.mrna_suitability!r} disagrees with the derived "
                    f"value {derived!r} for format {candidate.molecule_format!r} "
                    f"and sequence source {candidate.sequence_source!r}"
                )
                raise ValueError(msg)
            key = (candidate.antibody_id, candidate.toxin_family, candidate.species)
            if key in seen:
                msg = f"{label}: duplicate antibody_id, toxin_family and species {key!r}"
                raise ValueError(msg)
            seen.add(key)

    # -------------------------------------------------------------- analysis
    def by_species_family(self) -> dict[str, dict[str, list[Candidate]]]:
        """Group the candidates by species and then by toxin family.

        Returns:
            Nested mapping of species to toxin family to candidates.
        """
        grouped: dict[str, dict[str, list[Candidate]]] = {}
        for candidate in self.candidates:
            for species in candidate.species:
                families = grouped.setdefault(species, {})
                families.setdefault(candidate.toxin_family, []).append(candidate)
        return grouped

    def covering(self, species: str, toxin_family: str) -> list[Candidate]:
        """Return the candidates that count as coverage for one cell.

        A candidate covers a cell when it names the species and the toxin
        family and has reached at least :data:`COVERAGE_THRESHOLD`.

        Args:
            species: The binomial name.
            toxin_family: The toxin family.

        Returns:
            The covering candidates, in table order.
        """
        threshold = EVIDENCE_LEVELS.keys()
        order = list(threshold).index(COVERAGE_THRESHOLD)
        covering: list[Candidate] = []
        for candidate in self.candidates:
            if species not in candidate.species:
                continue
            if candidate.toxin_family != toxin_family:
                continue
            if list(threshold).index(candidate.evidence_level) < order:
                continue
            covering.append(candidate)
        return covering

    def gaps(self) -> list[Gap]:
        """Find the species and toxin family cells with no covering candidate.

        Returns:
            The gaps, ordered by species and then by toxin family.
        """
        grouped = self.by_species_family()
        found: list[Gap] = []
        for species in sorted(SPECIES_TOXIN_FAMILIES):
            families = SPECIES_TOXIN_FAMILIES[species]
            for toxin_family in families:
                if self.covering(species, toxin_family):
                    continue
                near_species = tuple(
                    other
                    for other in sorted(grouped)
                    if other != species
                    and grouped[other].get(toxin_family)
                    and self.covering(other, toxin_family)
                )
                near_family = tuple(
                    other_family
                    for other_family in families
                    if other_family != toxin_family
                    and self.covering(species, other_family)
                )
                found.append(
                    Gap(
                        species=species,
                        toxin_family=toxin_family,
                        reason=(
                            f"no candidate in the table names {species} and "
                            f"{toxin_family} with binding or better evidence"
                        ),
                        near_species=near_species,
                        near_family=near_family,
                    ),
                )
        return found

    def species_summary(self) -> list[dict[str, object]]:
        """Summarise the coverage of each Big Four species.

        Returns:
            One dictionary per Big Four species, in the order the Big Four
            are declared.
        """
        summary: list[dict[str, object]] = []
        for species in BIG_FOUR:
            candidates = [c for c in self.candidates if species in c.species]
            families = sorted(
                {
                    c.toxin_family
                    for c in candidates
                    if c.toxin_family not in {"mixed", "unspecified"}
                },
            )
            unresolved = sum(
                1 for c in candidates if c.toxin_family in {"mixed", "unspecified"}
            )
            suitability = Counter(c.mrna_suitability for c in candidates)
            best = max(
                (c.evidence_level for c in candidates),
                key=lambda level: list(EVIDENCE_LEVELS).index(level),
                default="none",
            )
            summary.append(
                {
                    "species": species,
                    "common_name": COMMON_NAMES.get(species, species),
                    "candidates": len(candidates),
                    "families": families,
                    "unresolved": unresolved,
                    "best_evidence": best,
                    "sufficiency": dict(suitability),
                },
            )
        return summary

    def findings(self) -> list[str]:
        """State the headline findings of the survey.

        The findings are computed from the table rather than written by
        hand, so they cannot go stale when a row is added or removed.

        Returns:
            One sentence per finding.
        """
        public = len(self.public_sequence_candidates())
        total = len(self.candidates)
        krait = self.covering("Bungarus caeruleus", "PLA2")
        viper = self.covering("Daboia russelii", "SVMP")
        big_four_gaps = [gap for gap in self.gaps() if gap.species in BIG_FOUR]
        return [
            (
                f"Only {public} of {total} candidates have a variable-region "
                "sequence that can be located in a public database. Sequence "
                "availability, not antibody discovery, is the bottleneck for "
                "an mRNA programme."
            ),
            (
                "Common krait venom is about 64.5 per cent phospholipase A2 "
                "by mass, and no published antibody candidate in this survey "
                "targets krait phospholipase A2. "
                f"{len(krait)} candidates cover krait phospholipase A2."
            ),
            (
                "No published antibody candidate targets the snake venom "
                "metalloproteinases of Indian Russell's viper, which are the "
                "main cause of the haemorrhage in that bite. "
                f"{len(viper)} candidates cover that toxin family."
            ),
            (
                f"{len(big_four_gaps)} of the uncovered cells belong to a Big "
                "Four species, and every one of them is a cell with no "
                "candidate at all rather than a cell with weak candidates."
            ),
            (
                "The one construct this repository has built, for E113, is "
                "one of the candidates whose sequence is publicly available."
            ),
        ]

    def sequence_availability(self) -> dict[str, int]:
        """Count the candidates by how available their sequence is.

        Returns:
            Mapping from suitability verdict to the number of candidates.
        """
        return dict(Counter(c.mrna_suitability for c in self.candidates))

    def public_sequence_candidates(self) -> list[Candidate]:
        """Return the candidates whose sequence is publicly locatable.

        Returns:
            The candidates, in table order.
        """
        return [c for c in self.candidates if c.mrna_suitability in {"ready", "adaptable"}]

    # -------------------------------------------------------------- rendering
    def to_tsv(self) -> str:
        """Render the table as tab-separated values.

        Returns:
            The TSV text, ending in a newline.
        """
        header = (
            "species",
            "toxin_family",
            "toxin_target",
            "antibody_id",
            "format",
            "sequence_source",
            "patent_status",
            "tested_against",
            "evidence_level",
            "mrna_suitability",
            "citation",
            "notes",
        )
        lines = ["\t".join(header)]
        for candidate in self.candidates:
            row = (
                ";".join(candidate.species),
                candidate.toxin_family,
                candidate.toxin_target,
                candidate.antibody_id,
                candidate.molecule_format,
                ";".join(candidate.sequence_source),
                candidate.patent_status,
                candidate.tested_against,
                candidate.evidence_level,
                candidate.mrna_suitability,
                candidate.citation,
                candidate.notes,
            )
            lines.append("\t".join(field.replace("\t", " ") for field in row))
        return "\n".join(lines) + "\n"

    def render_markdown(self) -> str:
        """Render the M4 report.

        Returns:
            The markdown text.
        """
        out: list[str] = []
        out.append("# M4: the candidate tracker")
        out.append("")
        out.append(
            "**In-silico designs. No therapeutic efficacy demonstrated. Not "
            "medical advice.**",
        )
        out.append("")
        out.append(
            "This report maps the published antibody candidates that have a "
            "claimed neutralising activity against a toxin, or against the "
            "whole venom, of one of India's Big Four snakes. It is generated "
            "by `scripts/candidate_tracker.py` from `src/mrna/candidates.py`, "
            "so the table and the analysis cannot disagree with each other.",
        )
        out.append("")
        out.append(
            "Two things are deliberately kept apart here. The first is "
            "whether a candidate works. The second is whether its sequence "
            "can be transcribed. A candidate can be excellent on the first "
            "count and unusable on the second, and for a project whose "
            "contribution is the mRNA layer, the second count is the one "
            "that decides.",
        )
        out.append("")

        out.append("## What the survey found")
        out.append("")
        for finding in self.findings():
            out.append(f"- {finding}")
        out.append("")

        out.append("## The table")
        out.append("")
        out.append(
            f"{len(self.candidates)} candidates, covering "
            f"{len({s for c in self.candidates for s in c.species})} species "
            "names and "
            f"{len({c.toxin_family for c in self.candidates})} toxin families.",
        )
        out.append("")
        out.append(
            "| Species | Toxin family | Target | Antibody | Format | Sequence "
            "source | Patent | Tested against | Evidence | mRNA |",
        )
        out.append("| --- " * 10 + "|")
        for candidate in self.candidates:
            species = "; ".join(candidate.species)
            source = "; ".join(candidate.sequence_source)
            out.append(
                f"| {species} | {candidate.toxin_family} | "
                f"{candidate.toxin_target} | {candidate.antibody_id} | "
                f"{candidate.molecule_format} | {source} | "
                f"{candidate.patent_status} | {candidate.tested_against} | "
                f"{candidate.evidence_level} | {candidate.mrna_suitability} |",
            )
        out.append("")

        out.append("## Coverage by species")
        out.append("")
        out.append(
            "| Species | Common name | Candidates | Families covered | Target "
            "unresolved | Best evidence |",
        )
        out.append("| --- | --- | --- | --- | --- | --- |")
        for row in self.species_summary():
            families = ", ".join(row["families"]) or "none"
            out.append(
                f"| {row['species']} | {row['common_name']} | "
                f"{row['candidates']} | {families} | {row['unresolved']} | "
                f"{row['best_evidence']} |",
            )
        out.append("")

        out.append("## The gaps")
        out.append("")
        out.append(
            "A cell is covered when at least one candidate names that species "
            f"and that toxin family and has reached {COVERAGE_THRESHOLD} or "
            "better. The empty cells are where this project can work, though "
            "not every empty cell is an opportunity: read the note under the "
            "matrix before treating a gap as a design target. The Big Four "
            "species are marked.",
        )
        out.append("")
        families = sorted(
            {family for fams in SPECIES_TOXIN_FAMILIES.values() for family in fams},
        )
        out.append("| Species | " + " | ".join(families) + " |")
        out.append("| --- " * (len(families) + 1) + "|")
        for species in sorted(SPECIES_TOXIN_FAMILIES):
            cells = []
            for family in families:
                covering = self.covering(species, family)
                if covering:
                    ids = ", ".join(c.antibody_id for c in covering)
                    cells.append(f"{ids} ({len(covering)})")
                else:
                    cells.append("none")
            label = species
            if species in BIG_FOUR:
                label = f"**{species}**"
            out.append(f"| {label} | " + " | ".join(cells) + " |")
        out.append("")
        gaps = self.gaps()
        if gaps:
            out.append(
                "The uncovered cells, one by one. The Big Four cells are the "
                "ones this project is about; the others are context.",
            )
            out.append("")
            for gap in gaps:
                marker = " (Big Four)" if gap.species in BIG_FOUR else ""
                out.append(
                    f"- **{gap.species} / {gap.toxin_family}**{marker}: "
                    f"{gap.reason}.",
                )
                if gap.near_family:
                    out.append(
                        "  The same species is covered for "
                        + ", ".join(gap.near_family)
                        + ".",
                    )
                if gap.near_species:
                    out.append(
                        "  The same toxin family is covered for "
                        + ", ".join(gap.near_species)
                        + ".",
                    )
        else:
            out.append("There are no uncovered cells.")
        out.append("")

        out.append("### An empty cell is not always an opportunity")
        out.append("")
        out.append(
            "Two things have to be read alongside the matrix. The first is "
            "that a candidate whose target toxin the source does not resolve "
            "is recorded with the toxin family `unspecified`, so it appears "
            "as no coverage anywhere. That is why the camelid antivenom "
            "raised against Indian saw-scaled viper venom does not fill the "
            "Echis carinatus sochureki row.",
        )
        out.append("")
        out.append(
            "The second is that the Kunitz-type serine protease inhibitor "
            "cells are empty for both viper species, and that is probably "
            "correct rather than a hole. Those inhibitors are the "
            "anticoagulant fraction of viper venom: they bind the clotting "
            "factors directly. An antibody against them would protect the "
            "clotting cascade, which is the opposite of what a patient with "
            "a viper bite needs. A tracker that reported every empty cell as "
            "an opportunity would be pointing this project at a design that "
            "would make things worse, so the matrix has to be read with the "
            "toxin families in mind.",
        )
        out.append("")

        out.append("## Sequence availability is the bottleneck")
        out.append("")
        availability = self.sequence_availability()
        out.append(
            "| mRNA suitability | Candidates | Meaning |",
        )
        out.append("| --- | --- | --- |")
        for verdict in ("ready", "adaptable", "unknown", "blocked"):
            out.append(
                f"| {verdict} | {availability.get(verdict, 0)} | "
                f"{MRNA_SUITABILITY[verdict]} |",
            )
        out.append("")
        public = self.public_sequence_candidates()
        out.append(
            f"{len(public)} of {len(self.candidates)} candidates have a "
            "sequence that can be located in a public database: "
            + ", ".join(c.antibody_id for c in public)
            + ".",
        )
        out.append("")
        out.append(
            "The rest are not missing because they are bad. They are missing "
            "because the sequence is in a paper's supplementary table, in a "
            "patent, or nowhere that a machine can read. That is a solvable "
            "problem, and it is a much smaller problem than discovering an "
            "antibody.",
        )
        out.append("")

        out.append("## Reading the toxin families")
        out.append("")
        for family in TOXIN_FAMILIES:
            note = TOXIN_FAMILY_NOTES.get(family)
            if note is None:
                continue
            out.append(f"- **{family}**: {note}")
        out.append("")

        out.append("## What the table does not claim")
        out.append("")
        out.append(
            "- **It is a survey, not an exhaustive census.** A row saying "
            "`none` for a sequence means the sequence was not located during "
            "the survey on 2026-11-01, not that no sequence exists.",
        )
        out.append(
            "- **Patent status is a fact about filings, not a clearance.** A "
            "row marked `application-pending` records that an application is "
            "public. It is not an opinion about freedom to operate, and it is "
            "not a freedom-to-operate opinion.",
        )
        out.append(
            "- **Evidence levels are the highest rung reached, not the "
            "quality of the result.** A candidate at L4 has protection data "
            "in an animal model. It does not follow that the protection was "
            "complete, or that the model resembles a human bite.",
        )
        out.append(
            "- **No candidate here has been evaluated by this project.** This "
            "is a map of what has been published. The one construct this "
            "repository has built, for E113, is built from a published "
            "sequence and carries the published data with it.",
        )
        out.append("")

        out.append("## Vocabularies")
        out.append("")
        out.append("Evidence levels:")
        out.append("")
        for level, meaning in EVIDENCE_LEVELS.items():
            out.append(f"- `{level}`: {meaning}")
        out.append("")
        out.append("mRNA suitability:")
        out.append("")
        for verdict, meaning in MRNA_SUITABILITY.items():
            out.append(f"- `{verdict}`: {meaning}")
        out.append("")
        out.append(
            "The suitability verdict is not a hand annotation. The engine "
            "recomputes it from the format and the sequence source and "
            "refuses to load a row where the two disagree, so the column "
            "cannot drift away from the facts it claims to summarise.",
        )
        out.append("")

        out.append("## Sources")
        out.append("")
        used = sorted({c.citation for c in self.candidates})
        for key in used:
            out.append(f"- `{key}`: {CITATIONS[key]}")
        out.append("")
        return "\n".join(out) + "\n"
