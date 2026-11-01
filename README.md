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
   from. The first one has landed in [M3](docs/m3_ecarpholin_s_construct.md):
   an E113 scFv against ecarpholin S from the saw-scaled viper.
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
| M1 | Codon optimiser (CAI, GC, motif avoidance) with tests | Second artifact, first repository views | Done, see [the M1 report](docs/m1_codon_optimizer.md) |
| M2 | Full-construct assembler, curated UTR library, secondary-structure checks | First tool someone else can use | Done, see [the M2 report](docs/m2_construct_assembler.md) |
| M3 | One toxin, one published neutralising antibody, one full construct, one design document | First attention from the venom research community | Done, see [the M3 report](docs/m3_ecarpholin_s_construct.md) |
| M4 | Literature map: Big Four toxins to antibodies to sequences to patents | First chance of being cited | Done, see [the M4 report](docs/m4_candidate_tracker.md) |
| M5 | The Big Four cocktail, planned from the literature map and built where it can be | First design that names its own blockers | Done, see [the M5 report](docs/m5_big_four_cocktail.md) |
| M6+ | Regional variants, partner outreach, design iteration | Sustained cadence | Not started |

## What M5 established

The Big Four cocktail, and why the published record cannot supply one.

M5 asks the question M4's map was built to answer: if you had to build one
antibody cocktail covering India's Big Four snakes, which payloads would it
contain, and what would stop you? The planner in
[`src/mrna/cocktail.py`](src/mrna/cocktail.py) treats it as a set cover over a
matrix of species and toxin families, searches it exhaustively rather than
greedily, and then tightens the question one constraint at a time. The
attrition down those tiers is the result.

- **The published record runs out long before the design does.** Four tiers, and
  the cells the plan can cover fall from 9 of 14 to 6, then 2, then 1. The
  questions are: is there a published candidate at all, can its sequence be
  located in a public database, is it one polypeptide a single transcript can
  carry, and has a patent search found nothing covering it.
- **Indian Russell's viper leaves the plan entirely.** The one recombinant
  candidate against the dominant toxin of that venom, the camelid VHH MKSS1,
  has no publicly locatable sequence, so the whole species drops out at the
  second tier. That is the single most actionable line in the milestone, and
  it is a data-availability problem rather than a discovery problem.
- **Three of the five cells no tier covers are viper toxin families that
  nobody has taken on.** The two viper species of the Big Four between them
  carry every toxin class that makes a viper bite dangerous, and not one of
  those classes has a published antibody candidate in the table.
- **Common krait phospholipase A2 has no candidate at any tier.** Krait venom
  is about 64.5 per cent phospholipase A2 by mass, so the most under-covered
  species of the Big Four is the one where a single good antibody would buy
  the most.
- **Two constructs can be built from public material.** E113, the scFv already
  published in M3, reproduced byte for byte by the M5 pipeline, and
  TPL1158_01_C09, a 123-residue camelid VHH against alpha-cobratoxin from
  PDB 9GCN. The other four payloads are Fabs or two-chain antibodies:
  converting them into one transcript means inventing a single-chain format
  that the published neutralisation data does not describe, so the data would
  no longer belong to the molecule being transcribed.
- **The licence is the real blocker, and it is stated as a fact.** A public
  patent search finds an application covering TPL1158_01_C09: it is named 71
  times in WO2025078529A1, published by Danmarks Tekniske Universitet. For
  E113 no application was located. That is a record of what the public patent
  record said on a date about specific applications, not a freedom-to-operate
  opinion, and it means an open programme built on these sequences is not free
  to exist.
- **Three cells are empty on purpose.** The Kunitz-type serine protease
  inhibitors are the anticoagulant fraction of viper venom, so an antibody
  against them would make the coagulopathy worse. The planner says why rather
  than reporting them as opportunities.

The milestone ships two full constructs, an element table, and a report whose
every number is computed from the M4 table rather than written by hand.

## What M4 established

The literature map, and the surprising part of it.

M4 records every published antibody candidate with a claimed neutralising
activity against a toxin, or the whole venom, of a Big Four species: 22
candidates, with the species, toxin family, target, antibody, format, sequence
source, patent status, what each was tested against, its evidence level and
its mRNA suitability. The table is generated from
[`src/mrna/candidates.py`](src/mrna/candidates.py) through
[`src/mrna/tracker.py`](src/mrna/tracker.py), so it cannot drift away from
the code that validates it.

- **Sequence availability, not antibody discovery, is the bottleneck.** Only
  8 of the 22 candidates have a variable-region sequence that can be located
  in a public database. The rest are in a paper's supplementary table, in a
  patent, or nowhere a machine can read. For a project whose contribution is
  the mRNA layer, that is the number that matters, and it is much easier to
  move than discovering an antibody.
- **Common krait is almost entirely unprotected.** Krait venom is about
  64.5 per cent phospholipase A2 by mass, and no published antibody
  candidate in the survey targets krait phospholipase A2. The one krait
  candidate is a three-finger toxin antibody, and it only doubled survival
  time against whole krait venom.
- **Russell's viper haemorrhage has no antibody candidate.** No published
  antibody targets the snake venom metalloproteinases of Indian *Daboia
  russelii*, which are the main cause of the bleeding in that bite. The one
  recombinant candidate in the survey, a camelid VHH, targets secretory
  phospholipase A2 instead.
- **Seven of the uncovered cells belong to a Big Four species.** The report
  marks them, and it also says which empty cells are *not* an opportunity:
  the Kunitz-type serine protease inhibitors are the anticoagulant fraction
  of viper venom, so an antibody against them would make the coagulopathy
  worse.
- **The mRNA suitability column is computed, not annotated.** The engine
  recomputes it from the format and the sequence source and refuses to load a
  row where the two disagree, so the column cannot drift away from the facts
  it claims to summarise.

## What M3 established

The first venom construct in the repository. It encodes a published human
single-chain variable fragment against a toxin of the saw-scaled viper, on
the cap, UTR and poly(A) scaffold that M2 built and measured against a
licensed vaccine.

- **The toxin is ecarpholin S** (*Echis carinatus*, sochurek), the myotoxic
  Ser49 phospholipase A2, UniProt P48650, structure PDB 2QHE. The UniProt
  sequence and the crystal structure are the same molecule, checked
  base for base.
- **The antibody is E113**, a human scFv, PDB 8IA6, from Kumar, Madni,
  Chaturvedi and Salunke, *Molecular Immunology* 2024. The deposited entry
  holds a crystallisation construct, so the vector residue and the
  hexahistidine tag are removed by the reference-data script, not by hand.
- **The construct is 1,233 nt** and reads line by line against Comirnaty:
  cap at 1 to 2, 5' UTR at 3 to 54, signal peptide at 55 to 108, payload at
  109 to 828, 3' UTR at 829 to 1123, poly(A) at 1124 to 1233.
- **The payload is a signal peptide plus a 239-residue scFv**, 257 residues
  in all, with the coding region at GC 0.6381 and CAI 0.9854 against the
  human codon table.
- **Two independent sources put the epitope on the part of the toxin that
  does the damage.** UniProt annotates residues 105 to 117 as "important
  for membrane-damaging activities", and the source paper's structural model
  places the antibody's CDR regions against that same C-terminal myotoxic
  loop.

The design document spends a section on what the construct is **not**, which
is the part worth reading: an scFv is monovalent and clears in hours, the
published data for this clone is modest on its own, and an mRNA-encoded
antibody is the wrong answer to an acute bite.

## What M2 established

The assembler takes the parts of a construct in order and produces two things:
the sequence, and the element table that describes it, in the format a
regulatory document uses.

- **It rebuilds Comirnaty from published parts.** 4,284 of 4,284 nucleotides,
  **zero mismatches**, and all six published coordinates reproduced exactly.
  Three nucleotides of the rebuild, 0.07 per cent, come from the document
  being reproduced: its second of two consecutive stop codons, which the
  published coding-sequence transcription stops short of. Everything else
  comes from the published 5' UTR, the published coding sequence and the
  published 3' UTR and poly(A) specification.
- **The cap is a chemical footprint, not two nucleotides of chemistry.** The
  document prints `cap G1A2 = m7G+m3'-5'-ppp-5'-Am2'-3'-p-`, so the same 5'
  UTR is quoted at 52, 53 or 54 nt across the literature. The engine fixes
  the document's convention and records the others.
- **The licensed product carries the forbidden polyadenylation signal twice,
  and never in its coding region.** Once in the 5' UTR, once in the 3' UTR.
  That is why motif scanning is reported per element: a single verdict for
  the whole construct would have failed a product two regulators approved.
- **The UTR library is borrowed, not invented.** The human alpha-globin 5'
  UTR with an optimised Kozak, and the AES plus mitochondrial 12S rRNA 3'
  UTR, both from Comirnaty, both cited, both confirmed against independent
  publications.

## What M1 established

The optimiser takes a protein and a human codon usage table and returns a
coding sequence that holds its GC content inside a band, contains no
forbidden motif, and is optimal for the Codon Adaptation Index subject to
those constraints. The benchmark is the SARS-CoV-2 spike protein:

- **78.2 per cent codon agreement** with Comirnaty's published design, from
  the protein alone. Our optimiser independently picks the same codon as
  BioNTech at roughly four positions in five.
- The band is derived from the approved products, not from a paper.
  Comirnaty's blocks span GC 0.400 to 0.700 and Spikevax's 0.467 to 0.767, so
  the enforced band is 0.40 to 0.70. The 0.50 to 0.60 band that appears in
  the literature would have excluded Spikevax.
- The viral gene itself scores CAI 0.6734 and GC 0.3734. That gap is the
  reason codon optimisation exists.

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
uv run python scripts/codon_optimizer.py --write   # refresh docs/m1_codon_optimizer.md
uv run python scripts/construct_assembler.py --write  # refresh docs/m2_construct_assembler.md
uv run python scripts/venom_construct.py --write      # refresh docs/m3_ecarpholin_s_construct.md
uv run python scripts/candidate_tracker.py --write    # refresh docs/m4_candidate_tracker.md
uv run python scripts/cocktail_design.py --write     # refresh docs/m5_big_four_cocktail.md
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
src/mrna/params.py       design policy constants, every value carrying a citation
src/mrna/policy.py       the design policy dataclass the optimiser obeys
src/mrna/motifs.py       forbidden motifs, each with a citation and a verified flag
src/mrna/optimize.py     the codon optimiser: block dynamic program plus repairs
src/mrna/utr.py          the UTR library: borrowed elements, each with its source
src/mrna/leader.py       the signal peptide library, for secreted payloads
src/mrna/assemble.py     the construct assembler and the regulatory-format element table
src/mrna/candidates.py   the curated candidate table for the Big Four, with its sources
src/mrna/tracker.py      the tracker engine: validation, coverage and gap analysis
src/mrna/payload.py      the antibody payloads a construct can carry, and their chains
src/mrna/cocktail.py     the cocktail planner: exact set cover over species and toxin families
src/mrna/reference.py    loaders for the committed reference data
src/mrna/fasta.py        minimal FASTA reader and writer
scripts/                 prepare_reference_data.py, vax_anatomy.py, codon_optimizer.py,
                         construct_assembler.py, venom_construct.py, candidate_tracker.py,
                         cocktail_design.py
data/reference/          committed reference sequences, each with provenance
data/designs/            the designed constructs, with their element tables
data/candidates.tsv      the published candidate table, generated from the engine
data/PROVENANCE.md       where every sequence came from, with checksums
data/reference/CHECKSUMS.sha256   checksums of the raw downloads
docs/m0_vax_anatomy.md   the M0 report, generated by the script so it cannot drift
docs/m1_codon_optimizer.md  the M1 report, same deal
docs/m2_construct_assembler.md  the M2 report, same deal
docs/m3_ecarpholin_s_construct.md  the M3 design document, same deal
docs/m4_candidate_tracker.md  the M4 candidate map, same deal
docs/m5_big_four_cocktail.md  the M5 cocktail plan and its two constructs, same deal
docs/m1_plan.md          the M1 plan, written before the code
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
