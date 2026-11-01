# Design log

One entry per decision that shaped the engine. Each entry is small on purpose:
a design log where nothing can be written is a design log that dies.

Entries are written in the order the decisions were made, with the source
that justifies them. The `Owner's note` lines are left for the project owner:
the biological reasoning behind a decision is the learning footprint of this
project and is written by a human, not by the tool that implements it.

---

## M0: reading a licensed vaccine backwards

### 2026-10-31, anchor the engine on a regulatory primary source

**Decision.** The calibration target for M0 is the WHO International
Nonproprietary Names Programme document 11889 for Comirnaty, which publishes
the full construct: a feature table with coordinates for every element, and
the complete 4,284-nt mRNA sequence with N1-methylpseudouridine marked.

**Why.** A regulatory document is a different kind of object from a GitHub
file. It has an issuing body, a document number, a date, and a table of
coordinates that can be checked. Anchoring the engine there means the first
thing the engine ever did was read the most authoritative public description
of an mRNA drug, not a blog post.

**Consequence.** The engine's job became concrete: parse the document,
reconstruct the construct, and reproduce every published coordinate. All six
elements match, and the poly(A) tail resolves to 30 adenosines, the 10-nt
linker `GCAUAUGACU`, and 70 more adenosines.

**Sources.** WHO INN document 11889, February 2021.

**Owner's note.** _to be written by the project owner_

---

### 2026-10-31, a second source is not optional

**Decision.** Every sequence used in M0 is checked against a second,
independently produced source, and the number of mismatches is reported.

**Why.** The single most common failure in computational biology is a
plausible sequence that is wrong. The cure is not more careful reading, it is
a second pair of eyes that was not looking at the first pair's assumptions.

**Consequence.** Three checks are in the M0 report and all of them pass:
the primary document against a read-derived sequencing assembly, 4,175 nt
compared with zero mismatches; the document's viral S coding sequence against
UniProt P0DTC2, an identical protein; and the antigen against the wild-type
protein, exactly two substitutions.

**Sources.** NAalytics et al. 2021 assemblies; UniProt P0DTC2.

**Owner's note.** _to be written by the project owner_

---

### 2026-10-31, validate the reference before using it

**Decision.** Before any comparison, the engine checks that the viral S coding
sequence taken from the primary document translates to exactly the UniProt
reference protein.

**Why.** If the reference were wrong, every diff downstream would be wrong in
a way that still looks like a result. This check costs one function call and
it turns the reference from an assumption into a verified input.

**Consequence.** The two agree across all 1,273 residues, so the K986P and
V987P differences that the engine reports are differences between the vaccine
and the virus, not artefacts.

**Sources.** UniProt P0DTC2; WHO INN document 11889.

**Owner's note.** _to be written by the project owner_

---

### 2026-10-31, read the "2P" construct correctly

**Decision.** The antigen difference is reported as K986P and V987P, located
at the junction between the HR1 helix and the central helix of the S2
subunit.

**Why.** The construct is universally called "2P" or "S2P", and it is easy to
read that as "serine at position 2 changed to proline". It is not. The two
prolines are at 986 and 987, and the substitution abolishes a key tryptic
cleavage site, which is what locks the spike in its prefusion shape. Getting
this naming right matters here: the project will talk about stabilised
antigens constantly, and a wrong shorthand becomes a wrong design.

**Consequence.** The engine reports positions and amino acids from the
sequences, never from a name, and the report cites the TGA and MHRA
regulatory assessments that use the same numbering.

**Sources.** TGA Nonclinical Evaluation Report for BNT162b2; MHRA Public
Assessment Report for BNT162b2.

**Owner's note.** _to be written by the project owner_

---

### 2026-10-31, implement the CAI rather than depend on it

**Decision.** The Codon Adaptation Index is implemented in
`src/mrna/anatomy.py` from the published definition, and no third-party CAI
package is a dependency.

**Why.** The obvious choice was the `CAI` package from PyPI. It does not
build on modern Python: version 1.0.3 uses the removed `use_2to3` setuptools
feature and fails during resolution. Beyond that, an engine that publishes
numbers should be able to show the arithmetic behind them. Fifteen lines
with a citation is worth more than a dependency that may stop installing.

**Consequence.** The implementation follows Sharp and Li 1987 exactly:
relative adaptiveness is the usage fraction divided by the best usage for
that amino acid, and the CAI is the geometric mean over the sense codons,
excluding the initiator ATG and the stop codon. The test suite pins the
identity case, where a sequence built only from preferred codons must score
exactly 1.0.

**Sources.** Sharp PM, Li WH. 1987. Mol Evol Biol 26(3):151-161.

**Owner's note.** _to be written by the project owner_

---

### 2026-10-31, calibrate the design targets on licensed products

**Decision.** The GC band and the CAI target in `src/mrna/params.py` are
justified against two licensed mRNA vaccines rather than against a
textbook median.

**Why.** The GC target of 0.55 was going to be cited to a literature
statement that has not been checked against the primary source. Measuring the
two products that people actually use gives a better answer than a citation
that might be about something else: BNT162b2 lands at GC 0.5699 and CAI
0.9515, and the other vaccine lands at GC 0.6230 and CAI 0.9786.

**Consequence.** Both sit inside the accepted band of 0.50 to 0.60 for
BNT162b2 and just above it for the other, which is a useful warning: a band
tuned on one product can exclude a second, equally good one. The band is
therefore recorded as a target with a cited median, and every construct that
uses it also reports its measured value.

**Sources.** WHO INN document 11889; Kazusa codon usage table; the M0
report in `docs/m0_vax_anatomy.md`.

**Owner's note.** _to be written by the project owner_

---

### 2026-10-31, the ORF closes on the first stop codon

**Decision.** `find_orf` returns the coding sequence up to and including the
first in-frame stop codon.

**Why.** Comirnaty has two consecutive stop codons at 3874 to 3879, and a
tool that counted to the end of both would report a coding sequence three
bases longer than the one that is actually translated. The first stop is the
one that terminates translation.

**Consequence.** The reported coding sequence is 3,822 nt, which is 1,273
residues plus one stop codon, and the published signal peptide start at
nucleotide 55 is reproduced exactly.

**Sources.** WHO INN document 11889, feature table.

**Owner's note.** _to be written by the project owner_

---

## M1: the codon optimiser

### 2026-11-01, the GC band comes from the approved products, not from a paper

**Decision.** The enforced GC band is 0.40 to 0.70 per 30-nucleotide block.

**Why.** The plan assumed the 0.50 to 0.60 band that appears in the
literature. Measuring the two approved products showed their blocks span
0.400 to 0.700 for Comirnaty and 0.467 to 0.767 for Spikevax. A band of
0.50 to 0.60 would have excluded one of the two products a regulator
approved. A constraint that rejects an approved product is a literature
number wearing a constraint's clothes.

**Consequence.** The enforced band is the band both approved products
satisfy, so a design inside it is never more extreme than something already
in clinical use. The narrower figure stays in `params.py` as a documented
target, not a constraint.

**Sources.** `docs/m0_vax_anatomy.md`, `docs/m1_codon_optimizer.md`.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, maximising the CAI alone is trivial, and that is not the point

**Decision.** The optimiser maximises the CAI *subject to* the GC band and
the motif constraints, and the report leads with the constraints rather than
with the CAI.

**Why.** Choosing the most favoured codon for every amino acid already gives
a CAI of 1.0000. There is no design content in that number. The content is in
what has to give way to satisfy a GC profile and a manufacturing constraint,
and in being able to say afterwards that nothing gave way unnecessarily.

**Consequence.** The claim in the design document is "the CAI is optimal
subject to the GC band", which is a different and stronger sentence than "a
greedy pass found something good". It is backed by a block dynamic program
rather than by a heuristic.

**Sources.** Sharp PM, Li WH. 1987. Mol Evol Biol 26(3):151-161.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, verify a motif or keep it out

**Decision.** Only the AAUAAA polyadenylation signal is enforced as a motif.
The homopolymer limit of 8 bases and the poly-lysine limit of 4 codons are
enforced as engineering defaults. The restriction sites are documented as
candidates and are not enforced.

**Why.** The rule is that nothing enters the engine without a source. The
polyadenylation signal has two primary papers behind it. The homopolymer and
poly-lysine limits are justified by mechanism, and the thresholds are chosen
from what the approved products actually do: neither contains a homopolymer
run of nine or more, nor even two consecutive AAA codons, so the limits sit
comfortably outside both. The restriction sites depend on which DNA template
a manufacturing partner will use, and this project does not know that yet, so
enforcing a guess would be exactly the sort of plausible-sounding fiction
this repository exists to avoid.

**Consequence.** `mrna.motifs` carries a verified flag, and a test asserts
that every default motif is verified and every candidate is not. The repair
passes are exercised against a sequence that contains the signal directly.

**Sources.** Proudfoot NJ, Brownlee BG. 1976. Nucleic Acids Res 3(7):1183-1190.
Fitzgerald PC, Shenk T. 1981. Proc Natl Acad Sci USA 78(7):4182-4185.
Mechanism of ribosome stalling during translation of a poly(A) tail,
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6900289/

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, the benchmark is a published design

**Decision.** The headline measurement is the codon agreement between our
optimiser's output for the spike protein and Comirnaty's published coding
sequence.

**Why.** Every other number in the report is a number we produced, which a
reader has to take on trust. This one is a comparison against a sequence that
was published by a regulator and that anyone can download. Given only the
protein and a codon table, the optimiser independently chooses the same codon
as BioNTech at 78.2 per cent of 1,273 positions.

**Consequence.** The viral gene scores CAI 0.6734 and GC 0.3734 against the
optimised 0.9981 and 0.6140. That gap is the entire reason codon optimisation
exists, and it is now a table in a public document rather than a claim in a
conversation.

**Sources.** WHO INN document 11889; UniProt P0DTC2; the NAalytics
assembly.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, a repair may not make an infeasible block worse

**Decision.** A repair swap must keep its GC block inside the band. If the
block was already outside the band because the protein admits no other
encoding, the swap is allowed provided it does not push the block further
out.

**Why.** The first version of the rule was absolute, and it had a bad
consequence: in a protein whose GC content cannot reach the band, no repair
was ever applied, so a forbidden motif survived. A motif left in place is a
worse defect than a band that was already out of reach.

**Consequence.** The rule is now relative, and the test suite covers both
halves: a motif in a repairable sequence is removed, and an infeasible
block is reported rather than hidden.

**Sources.** `src/mrna/optimize.py`, `tests/test_optimize.py`.

**Owner's note.** _to be written by the project owner_

---

## Next

The open items for milestone M2 are the full-construct assembler and the UTR
library, whose first entries will be the human alpha-globin 5' UTR and the
AES plus 12S rRNA 3' UTR that Comirnaty uses, both with their citations. The
assembler should emit the same element table that the WHO document uses, so
that a venom construct can be compared to a licensed product line by line.
