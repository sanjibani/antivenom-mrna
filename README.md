# antivenom-mrna

**Open, reproducible, full-construct mRNA designs for India's Big Four venomous snakes, and the auditable design engine that produced them.**

> **In-silico designs. No therapeutic efficacy demonstrated. Not medical advice.**

India records roughly 58,000 snakebite deaths a year, the highest of any
country, and the antivenom that exists does not fit the problem. Equine-derived
antivenom is made against a small number of venom pools, so it neutralises
poorly against the regional variation in Indian venoms, and it causes serum
sickness. Snakebite has been notifiable in India since Karnataka in February
2024 and nationally since November 2024, so the policy tailwind is real.

This repository is the attempt to build the thing nobody has published: an
open, reproducible set of full mRNA constructs for the Big Four venomous
snakes, plus the engine that produced them, auditable line by line.

## What is in here

Two artifacts, in this order of importance:

1. **The designs.** Full mRNA constructs for India-relevant venom targets:
   5' cap chemistry, 5' UTR, codon-optimised coding sequence, 3' UTR, poly(A)
   tail, and a provenance table that says where every antibody sequence came
   from. These arrive from milestone M3 onward.
2. **The engine.** A small, readable Python library that measures and builds
   mRNA constructs. It is useful on its own, and it is how the designs can be
   checked by anyone who does not trust the author.

## The honest part, stated up front

Intravenous antivenom works within minutes. An mRNA-encoded antibody is not
detectable for one to two hours and peaks somewhere between 12 and 24 hours.
So **mRNA antivenom is a poor answer to acute snakebite rescue**, and this
project does not claim otherwise.

What mRNA could plausibly be good for:

- **Prophylaxis for high-risk groups**, such as snake rescuers and
  agricultural workers in the highest-burden districts.
- **Regional manufacturing**, because the construct is synthetic and can be
  re-sequenced against a district's own venom profile instead of relying on a
  horse farm somewhere else.
- **Antibodies that are hard to manufacture biologically**, where expression
  from an mRNA template sidesteps cell-line or fermentation constraints.

Every design document in this repository leads with this caveat. It is the
thing that makes the rest of the work worth reading.

## What this project is not

- It does not discover new neutralising antibodies. That problem is
  computationally unsolved and is being worked on by better-equipped groups.
  This project reuses published sequences and designs the mRNA layer.
- It is not a clinical programme. There is no wet lab anywhere in the plan.

## Prior art, and what is left open

| Work | What it did | What is still open |
| --- | --- | --- |
| Almeida et al. 2025 | mRNA-LNP encoding toxin-specific antibodies, expressed and neutralised venom in animals | Not India-focused, not open |
| Baker lab, *Nature* 2025 | De novo AI-designed miniproteins neutralising 3FTx toxins, protecting mice | "AI-designed antivenom" is taken; this project makes no antibody-design claim |
| Glanville et al., *Cell* 2025 | Varespladib plus a designed protein cocktail, complete protection in animals | Different modality, proprietary |
| Sunagar (IISc) with Laustsen (DTU), 2026 | Camelid nanobodies cross-neutralising Indian cobra toxins | Nanobody format, not mRNA |
| Merck `mRNAid`, 2024 | Open-source mRNA design tool, CDS optimisation only | No UTR selection, no cap, no poly(A), no full-construct assembly; dormant since Jan 2024; non-commercial licence |
| Eterna OpenVaccine, Stanford Das lab, 2020 | Precedent that open mRNA designs are a legitimate genre | COVID-19 antigens, not venom |
| McLellan lab 2P construct, 2020 | A single openly published design artifact that four vaccine makers built on | Precedent that design can be the contribution |

**The open lane:** nobody has published open, reproducible, full-construct mRNA
designs for Indian snakebite. That is what this repository occupies. Regional
venom variation, mapped for the Indian spectacled cobra by Laxme and
colleagues in 2021, is treated as a first-class design input rather than a
footnote.

## Status

Milestones are sequenced, not scheduled. Every one ships an artifact, a public
post, and a validation signal, so that there are no silent stretches.

| Milestone | Artifact | Signal | Status |
| --- | --- | --- | --- |
| M0 | `vax_anatomy`: reverse-engineer the anatomy of a licensed mRNA vaccine | First public post; computed features match published values | Done, see [the M0 report](docs/m0_vax_anatomy.md) |
| M1 | Codon optimiser (CAI, GC, motif avoidance) with tests | Second artifact, first repository views | Not started |
| M2 | Full-construct assembler, curated UTR library, secondary-structure checks | First tool someone else can use | Not started |
| M3 | One toxin, one published neutralising antibody, one full construct, one design document | First attention from the venom research community | Not started |
| M4 | Literature map: Big Four toxins to antibodies to sequences to patents | First chance of being cited | Not started |
| M5 | Full Big Four coverage and a polyvalent cocktail design | Preprint and launch thread | Not started |
| M6+ | Regional variants, partner outreach, design iteration | Sustained cadence | Not started |

## What M0 established

M0 is done. The engine reads the WHO International Nonproprietary Names
Programme document for Comirnaty, reconstructs the full construct from it, and
checks every coordinate the document publishes. The report is
[`docs/m0_vax_anatomy.md`](docs/m0_vax_anatomy.md), and the headlines are:

- **Every published coordinate is reproduced.** Cap at 1 to 2, 5' UTR at 3 to
  54, signal peptide at 55 to 102, coding region to 3879, 3' UTR at 3880 to
  4174, poly(A) at 4175 to 4284. Total length 4,284 nt.
- **Three independent sources agree.** The regulatory document, an
  independently produced sequencing assembly of the commercial product
  (4,175 nt compared, zero mismatches), and the authoritative UniProt
  reference protein all match.
- **The antigen difference is exactly two residues**, K986P and V987P, the
  two prolines that lock the spike in its prefusion shape. The document, the
  TGA and MHRA regulatory assessments, and the engine's own diff all agree.
- **The design targets are calibrated on licensed products**, not on a
  textbook: GC 0.5699 and CAI 0.9515 for one vaccine, GC 0.6230 and CAI
  0.9786 for the other, which encode the *same* protein with different codon
  choices at 294 codons.

## Quickstart

```bash
uv sync                                  # create the environment
uv run python scripts/prepare_reference_data.py   # regenerate derived reference data
uv run pytest                            # run the test suite
uv run python scripts/vax_anatomy.py --write      # refresh docs/m0_vax_anatomy.md
```

Optional tools are installed on demand rather than declared as extras, so that
`uv sync` stays fast on every platform:

```bash
uv pip install ViennaRNA   # secondary-structure predictions, needed from M2
```

## Repository layout

```
src/mrna/sequence.py     primitive layer: codon table, translate, back-translate, GC
src/mrna/anatomy.py      measurement layer: ORFs, Kozak context, GC, CAI, protein diffs
src/mrna/params.py       design policy, every value carrying a citation
src/mrna/fasta.py        minimal FASTA reader and writer
scripts/                 prepare_reference_data.py, vax_anatomy.py
data/reference/          committed reference sequences, each with provenance
data/PROVENANCE.md       where every sequence came from, with checksums
data/reference/CHECKSUMS.sha256   checksums of the raw downloads
docs/m0_vax_anatomy.md   the M0 report, generated by the script so it cannot drift
design-log.md            the running log of design decisions
tests/                   pytest and hypothesis tests
```

## Rules this repository follows

- **No magic constants.** Every number in `src/mrna/params.py` carries the
  citation that justifies it. Every measurement is reproduced by a test.
- **Every sequence has a provenance record.** Sequences are CC-BY 4.0, code is
  MIT. Non-commercial terms are avoided deliberately: the closest comparable
  open tool is dormant under exactly that restriction.
- **Every artifact carries the disclaimer** at the top of this file.
- **Design is the contribution.** Validation runs on a ladder: L0 in silico
  (done here, rigorously and openly), then L1 expression, L2 binding, L3 ex
  vivo neutralisation, L4 animal ED50, L5 chemistry and control. Wet-lab
  partners take it from L1 upward. The McLellan 2P construct is the precedent
  that a single published design can change the field without a single
  experiment run by its authors.
- **The biological rationale is written by a human.** Code, tests and
  citations are machine-written. The reasoning about why a construct is
  designed the way it is, in the pull request body, is written by the project
  owner. That is the learning footprint, and it is the point.

## Contributing

Small, focused pull requests are welcome. Please read `design-log.md` first,
and cite every biological claim in the pull request body. If you find a
sequence or a number in this repository that does not match its source, that
is a bug worth reporting: the whole value of this project is that its numbers
can be checked.

## Citation

If you use the designs or the engine, please cite this repository and say
which version you used. A `CITATION.cff` file will be added once the designs
land.

## Licence

Code: MIT. Sequences and documentation: CC-BY 4.0. See `LICENSE`.
