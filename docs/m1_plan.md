# M1 plan: the codon optimiser

Status: ready to execute. Decisions that belong to the project owner are
listed at the bottom under "Open decisions".

## Goal

Take a protein, a host codon usage table and a design policy, and emit a
coding sequence that:

- translates back to exactly the input protein,
- sits inside the GC band, measured over fixed blocks and reported over
  sliding windows,
- contains no forbidden motif (homopolymer runs, cryptic polyadenylation
  signals, cloning sites),
- has a CAI at or above the target,
- is byte-identical on every run.

Ship signal: "I designed my first codon-optimised mRNA from scratch",
with before and after numbers and a benchmark against a licensed product.

## Why M1 is the optimiser and not something else

M0 gave us the measurements. Without a trustworthy scorer, an optimiser is
just a random walk with good manners. The scorer is already validated
against two approved products, so M1 can be scored honestly.

## Components

### 1. `src/mrna/policy.py`, the design policy dataclass

One object holding every knob, every knob cited:

| Field | Meaning | Source |
| --- | --- | --- |
| `gc_min`, `gc_max` | GC band for the coding region | measured from the two licensed products, see `docs/m0_vax_anatomy.md` |
| `gc_block_nt` | length of the block the band is enforced over | design choice, to be fixed |
| `cai_target` | minimum CAI | measured from the licensed products |
| `max_homopolymer` | longest tolerated single-base run | citation to be verified |
| `forbidden_motifs` | regex set | citations to be verified, see below |
| `avoid_upstream_aug` | reject internal AUG before the intended start | Kozak 1986, leaky scanning |

### 2. `src/mrna/optimize.py`, the optimiser

Proposed pipeline, three stages:

1. **Score.** For each amino acid, the candidate codons and their relative
   adaptiveness `w` from the host table, exactly as in `mrna.anatomy`.
2. **Choose.** A dynamic program over codons where the state is the position
   and the GC count inside the current block, and the objective is to
   maximise the sum of `log(w)`, which is the CAI. This gives a provable
   statement for the design document: *the CAI is optimal subject to the GC
   and motif constraints*, rather than *a greedy pass found something good*.
   State space is roughly 51 per codon position, so a 1,273-residue protein
   solves in milliseconds.
3. **Repair.** Scan for forbidden motifs and swap synonymous codons locally
   until none remain, then re-verify the GC band, because a swap can push a
   block out of band. Deterministic tie-breaking, no randomness.

Alternative if the DP proves slow or brittle in practice: greedy selection
with a GC steering controller plus the same repair pass. The DP is preferred
because it upgrades the claim we can make in a design document.

### 3. `src/mrna/motifs.py`, the forbidden motif set

Every entry carries a citation and a status flag, because some of these are
well established and some need verification before they go in a public
document.

| Motif | Why | Status |
| --- | --- | --- |
| `A{9,}` and other runs of 9 or more identical bases | synthesis and stability problems, and A-tracts are known to stall ribosomes | citation to verify |
| `AATAAA` | canonical polyadenylation signal; an internal copy causes premature cleavage | citation to verify |
| `GCGGCCGC` (NotI), `GGGCCC` (ApaI), `GCTCTTC` (SapI) | cloning sites in the linear DNA template | verify which the chosen template actually uses |
| `TAATACGACTCACTATAGGG` (T7 promoter) | must appear once, at the start of the template | standard |
| `GG` at the first transcribed base | efficient T7 initiation prefers a G at +1 | standard |

Nothing enters `motifs.py` without a citation. Where a citation is not yet
verified, the motif stays out of the default policy and lives in a clearly
labelled candidate list until it is checked.

## Benchmark and validation

The benchmark is free and it is strong: **the SARS-CoV-2 spike protein**.

BioNTech and Moderna both published their codon choices for the same 1,273
residue protein, differing at 294 codons. So we can run our optimiser on that
protein and ask a sharp question: *does it independently rediscover the codon
choices that regulators approved?*

| Benchmark | What it proves |
| --- | --- |
| Optimised spike translates to the spike | correctness |
| GC of every block inside the band | constraint satisfaction |
| CAI of our sequence against human Kazusa, compared with 0.9515 and 0.9786 | competitive with two approved products |
| Codon agreement with Comirnaty's published CDS | the headline number: how many codons we independently chose the same way |
| No forbidden motif remains | manufacturing readiness |

Second benchmark, for the public artifact: **GFP**, the classic 238-residue
protein, reported as a before and after table. The before is a
single-codon-per-amino-acid naive sequence, the after is our optimiser.

Third benchmark, on the measurement side: reproduce the Codon Adaptation
Index of a set of genes with published CAI values, to show the scorer agrees
with the literature rather than only with itself.

## Test plan

- output translates back to the input protein, for a set of real proteins
- output is byte-identical across runs
- every GC block inside the band
- no forbidden motif present, checked by an independent scan, not by the
  function that removed them
- CAI at or above target
- no internal AUG in the region before the intended start
- the repair pass cannot push a block out of band
- hypothesis tests: for random protein sequences, the invariants hold
- edge cases: proteins with long low-complexity stretches, proteins containing
  Met and Trp only options, very short proteins, and a protein whose only
  possible codons cannot satisfy the GC band, where the optimiser must fail
  loudly rather than silently emit a bad sequence

## Files to add

```
src/mrna/policy.py     the cited design policy dataclass
src/mrna/optimize.py   the optimiser
src/mrna/motifs.py     the forbidden motif set
tests/test_optimize.py
tests/test_motifs.py
docs/m1_codon_optimizer.md   the design report, generated so it cannot drift
```

## Explicitly out of scope for M1

- UTR selection and the full construct assembler, that is M2
- secondary structure prediction with ViennaRNA, also M2
- any venom construct, that is M3
- any claim about expression or efficacy, that is L1 and above and is not ours

## Open decisions for the project owner

1. **Algorithm.** Windowed dynamic program with optimality claim, or greedy
   with repair and no optimality claim. Recommendation: the dynamic program,
   because it changes what we can write in a design document.
2. **GC constraint granularity.** Enforced per 50-nt block, per 30-nt block,
   or per sliding window of 50 nt. The last is the strictest and the most
   expensive. Recommendation: enforce per 30-nt block and report the sliding
   window, because the published guidance is about the region near the start
   codon, not about every window.
3. **CAI target.** 0.90 minimum, or 0.95 minimum to match the licensed
   products. Recommendation: 0.95 as the target and 0.90 as the floor, so a
   construct that cannot reach it is reported rather than quietly shipped.
4. **Which host table.** Human, since the protein is made in human muscle
   and dendritic cells after injection. But if a partner manufactures in CHO
   or HEK cells, a CHO table changes the answer. Recommendation: human by
   default, host selectable, and say which host was used in every design
   document.
5. **Which protein for the public artifact.** GFP, because it is the classic
   teaching example and it is small. Recommendation: GFP, with the spike
   benchmark as the credibility half of the post.

## The post

> I wrote a codon optimiser for mRNA and immediately pointed it at the spike
> protein from a licensed vaccine. It translates back to the same 1,273
> amino acids, holds GC in band, and scores a CAI of [X] against the human
> codon table.
> Independently, it chose the same codon as the licensed product at [Y] per
> cent of positions.
> Same protein, different words.
> In-silico designs. No therapeutic efficacy demonstrated. Not medical advice.
