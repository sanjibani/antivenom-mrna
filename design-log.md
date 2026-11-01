# Design log

One entry per decision that shaped the engine. Each entry is small on purpose:
a design log where nothing can be written is a design log that dies.

Entries are written in the order the decisions were made, with the source
that justifies them. The `Owner's note` lines are left for the project owner:
the biological reasoning behind a decision is the learning footprint of this
project and is written by a human, not by the tool that implements it.

## Writing an owner's note

Every entry in this log ends with an `Owner's note` line, and every one of
those lines is still empty. They are empty on purpose, and they stay empty
until the project owner fills them in.

The reason is specific. The `Why` section above each entry is the reasoning
of the thing that implemented the decision: why a block dynamic program beats
a greedy pass, why a rebuild is a better test than a self-consistency check.
That reasoning belongs to a machine and it is already recorded. The `Owner's
note` is the part nobody else can write, which is what the decision assumes
about the biology and what would follow if the assumption turned out to be
wrong. A note written by the implementer would be a second summary of the
implementer's own work, and would teach the owner nothing.

A note is short. Three questions, answered in a few sentences each:

1. **What does this decision assume about the biology?**
2. **What would break if that assumption is wrong?**
3. **What would you do differently if you knew more?**

The template, to be copied for any entry:

> **Assumption.** ...
>
> **If wrong.** ...
>
> **With more information.** ...

To fill one in, replace the line

```
**Owner's note.** _to be written by the project owner_
```

in the relevant entry with your note under the heading. The note does not need
to be a biology essay, and it does not need to agree with the `Why` section.
If the reasoning disagrees with the implementation, that disagreement is the
most valuable sentence in this log, and it belongs in the log before it
belongs in a commit.

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

## M2: the assembler and the UTR library

### 2026-11-01, the cap is a chemical footprint, not a stretch of nucleotides

**Decision.** The assembler models the cap as a two-nucleotide footprint that
carries a chemical structure, and the 5' UTR element is numbered from
position 3, exactly as the WHO document numbers it.

**Why.** The document prints `cap G1A2 = m7G+m3'-5'-ppp-5'-Am2'-3'-p-`, which
says the 7-methylguanosine at position 1 is joined to the first backbone
nucleotide at position 2 by a 5'-to-5' triphosphate. So the cap occupies two
printed positions even though it is a modification rather than a sequence.
The consequence is that the same 5' UTR is quoted at 52 nt by the WHO
document, at 53 nt when the first backbone nucleotide is counted, and at 54 nt
when the cap guanosine is counted too. All three describe the same molecule,
and the engine had to pick one convention before it could emit a table.

**Consequence.** The engine uses the document's numbering, because the
document's coordinates are the ones a design has to match, and because the
rebuild test would fail against any other convention. The UTR library
records the alternative lengths in the note on the 5' UTR entry, so a reader
who counts by hand does not conclude the engine is wrong.

**Sources.** WHO INN document 11889, February 2021, section '5'-capping
structure'; Kariko K et al. 2005, Immunity 23(2):165-175; Kariko K et al.
2008, Mol Ther 16(11):1833-1840.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, the acceptance test is a rebuild, not a self-consistency check

**Decision.** The assembler's acceptance criterion is that it reproduces
Comirnaty from published parts: the same 4,284 nucleotides and the same six
rows of the element table, with the mismatch count reported.

**Why.** Every other number this repository publishes is a number we
produced, which a reader has to take on trust. A rebuild is different. The
inputs come from published sources, the expected output comes from a
regulatory document, and the two are compared without either side being
adjusted afterwards. A tool that can put a licensed product back together
from its parts is a tool whose output for an antivenom construct is worth
reading.

**Consequence.** The rebuild lands on 4,284 of 4,284 nucleotides with zero
mismatches and reproduces all six published coordinates exactly. Three of the
4,284 nucleotides, 0.07 per cent, come from the document being reproduced:
its second of two consecutive stop codons, which the published coding
sequence transcription stops short of. The report states that openly rather
than letting a zero look cleaner than it is.

**Sources.** `src/mrna/assemble.py`, `tests/test_assemble.py`,
`docs/m2_construct_assembler.md`.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, motif scope is per element, and a licensed product proves why

**Decision.** The assembler reports forbidden-motif occurrences per element,
never as a single verdict for the whole construct.

**Why.** The M1 motif policy was scoped to the coding region, and that scope
was correct. A polyadenylation signal inside a coding region can cause
premature cleavage and a transcript with no in-frame stop codon, which is a
real defect. The same signal in a 5' UTR has no cleavage site upstream of it
and is harmless. Applying the coding-region rule to the whole construct
would have failed a product that two regulators approved, which is exactly
the failure mode that makes a rule wrong rather than the rule being too
strict.

**Consequence.** Comirnaty carries the canonical AAUAAA signal exactly
twice, once in the 5' UTR at 4 to 9 and once in the 3' UTR at 4099 to 4104,
and exactly zero times in its coding region. The engine reports all three
facts, and the M2 report explains why the two in the UTRs are not defects.

**Sources.** Proudfoot NJ, Brownlee BG. 1976. Nucleic Acids Res
3(7):1183-1190. Fitzgerald PC, Shenk T. 1981. Proc Natl Acad Sci USA
78(7):4182-4185.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, upstream start codons belong to the UTR, not to the coding region

**Decision.** The check for ATG codons upstream of the intended start codon
lives in the assembler, where the real 5' UTR and the real Kozak context are
known, and it reports the context of each one rather than just counting them.

**Why.** The M1 amendment moved this check out of the optimiser on the
ground that it is a 5' UTR constraint. That was right about the location and
incomplete about the question. A start codon in the 5' UTR competes with the
intended one, and how well it competes depends on the Kozak context around
it: Kozak 1986 showed that the purine at position -3 and the G at position
+4 are the two positions that matter most. An upstream ATG in a weak
context is noise. An upstream ATG in a strong one is a design error.

**Consequence.** `Construct.upstream_starts()` returns every ATG in the 5'
UTR with its -3 and +4 bases and whether it matches the consensus at both.
Comirnaty's 5' UTR has none, so the new check does not retroactively
condemn a licensed product.

**Sources.** Kozak M. 1986. Cell 44(2):283-292.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, the UTR library holds borrowed elements, not invented ones

**Decision.** The UTR library ships two entries: the human alpha-globin 5'
UTR with an optimised Kozak sequence, and the AES plus mitochondrial 12S
rRNA 3' UTR. Both are the elements Comirnaty uses, and both carry the
citation that justifies them.

**Why.** The project's whole argument is that the design space is mostly
arithmetic and that the inputs should be published ones. Borrowing a UTR
that a regulator has already accepted inside a licensed product is a
smaller claim than inventing a better one, and it is the only claim this
project can support without a wet lab. It is also the claim that keeps the
antivenom construct comparable to a licensed product line by line, which is
the format the whole repository is written in.

**Consequence.** The 5' UTR is the 52 nucleotides the WHO document prints at
positions 3 to 54, independently confirmed by a peer-reviewed review. The
295-nucleotide 3' UTR decomposes into the 136-nt AES element and the
142-nt mtRNR1 element that Orlandini von Niessen and colleagues found by
cell-based screening, in the order the document names them, with the
document's own 6-nt prefix and 11-nt suffix either side.

**Sources.** WHO INN document 11889, February 2021; Orlandini von Niessen AG
et al. 2019, Mol Ther 27(4):824-836; Fang E et al. 2022, Signal Transduct
Target Ther 7(1):94.

**Owner's note.** _to be written by the project owner_

---

## M3: the first venom construct

### 2026-11-01, the payload is a published antibody, and the project does not discover antibodies

**Decision.** The first venom construct encodes E113, a human single-chain
variable fragment against *Echis carinatus* phospholipase A2, published by
Kumar, Madni, Chaturvedi and Salunke in *Molecular Immunology* in 2024 and
deposited as PDB entry 8IA6.

**Why.** Two things had to be true at once. The species had to be one of
India's Big Four, and the antibody sequence had to be retrievable from a
source with a permanent identifier, the way the Comirnaty construct was
retrievable from a WHO document in M0. E113 satisfies both. The alternative
candidates were rejected for concrete reasons: the broadly neutralising
Centivax antibodies are protected by patent families and come as Fabs,
which would need two open reading frames, and the DTU scFvs are published
without a deposited structure.

**Consequence.** The construct is a single polypeptide, one open reading
frame, 1,233 nt. The project claims no part in the antibody's discovery, and
the design document says so in the same words as the README.

**Sources.** Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2024. Mol Immunol
165:55-67. PDB entry 8IA6.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, the crystallisation tag is removed by script, not by hand

**Decision.** The payload is the deposited construct of PDB 8IA6 minus its
N-terminal vector residue and its C-terminal hexahistidine tag, and the
removal is performed by `scripts/prepare_reference_data.py` with an
assertion on both ends.

**Why.** The deposited entry holds the molecule that was crystallised, which
is the antibody with an expression handle on it. A therapeutic construct
should not encode a hexahistidine tag. But a trim done by hand is a trim
that nobody can check, and the crystallised construct is the sequence whose
structure was solved, so the reasoning about which residues are real has to
be recorded rather than assumed. The crystal locates residues 2 to 239 of
the deposited chain: residue 1, the vector alanine, was never located, and
so was the tag.

**Consequence.** `write_e113_payload()` refuses to run if the deposited
construct does not start with the vector residue, does not end with the tag,
or if the trimmed result does not start with the IGHV3-23 germline sequence
and end with the IGKC constant region. The payload is therefore pinned to
the deposited sequence by a test rather than by a copy.

**Sources.** PDB entry 8IA6; `scripts/prepare_reference_data.py`;
`tests/test_venom_construct.py`.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, the signal peptide is an element of its own

**Decision.** The N-terminal leader is a separate element of the construct
and a separate entry in the library, `src/mrna/leader.py`, rather than a
length argument passed to the assembler.

**Why.** A secreted antibody has to leave the cell, and a signal peptide is
what takes it there. Two properties made it worth its own module rather than
a string constant. First, the leader is a stretch of protein, not of DNA:
signal peptidase reads the protein and cuts it, and the DNA that encodes it
belongs to the codon optimiser. Second, the leader must be verified against
the coding region, not merely counted. Passing a `Leader` object means the
assembler translates the first codons of the coding region and checks that
they encode the leader it was given, so a mismatch between the design and
the DNA is an error rather than a silent off-by-one.

**Consequence.** The construct's element table has a `sig` row at 55 to 108,
exactly as Comirnaty's does at 55 to 102, and the engine raises an error if
the coding element does not encode the leader it was handed. The chosen
leader is the 18-residue human serum albumin signal peptide, which is one
of the two most commonly used in mRNA vaccine designs.

**Sources.** UniProt P02768; Attallah C et al. 2017, Protein Expr Purif
132:27-33; doi:10.1101/2025.11.25.690495.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, the design document has to say what the construct is not

**Decision.** The M3 document contains a section headed "What this construct
is not", and it is written before the section about what it is.

**Why.** This is the first document in the repository that a venom
researcher might read as a proposal, and the honest content of it is mostly
bad news: an scFv is monovalent, it clears in hours, the published data for
this clone is modest on its own, and an mRNA-encoded antibody is the wrong
answer to an acute bite. A reader who finds those points buried at the end
has been sold something. A reader who finds them in the fourth section
knows exactly what they are looking at.

**Consequence.** The document states that the binding and neutralisation data
belongs to the antibody made in bacteria, and that an scFv made in a
mammalian cell from an mRNA template is a different molecule until someone
shows it is not, because glycosylation, disulfide pairing and folding all
differ between the two expression systems. It also records that the patent
position is a negative search result rather than a freedom-to-operate
opinion.

**Sources.** `docs/m3_ecarpholin_s_construct.md`.

**Owner's note.** _to be written by the project owner_

---

### 2026-11-01, the toxin is checked against its own structure

**Decision.** The target of the antibody is taken from UniProt, and the
crystal structure of the same toxin is used as an independent check that the
two describe the same molecule.

**Why.** This is the rule M0 established, applied to the first venom design
rather than only to the calibration set. Ecarpholin S is available from two
places that were produced by different people for different reasons, and a
design document that names a toxin has to be able to say which toxin, down
to the base.

**Consequence.** The UniProt sequence and the PDB 2QHE sequence are
identical over all 122 residues, and the test suite checks it. The two
records also converge on the mechanism from different directions: UniProt
annotates residues 105 to 117 as important for membrane-damaging
activities, and the source paper places the antibody's CDR regions against
that same C-terminal myotoxic loop.

**Sources.** UniProt P48650; PDB 2QHE; Zhou X et al. 2008, Biophys J
95:3366-3380.

**Owner's note.** _to be written by the project owner_

---

## Next

The open items for milestone M4 are the tracker: every published candidate
in the field, recorded in one table with species, toxin family, target,
antibody, format, sequence source, patent status, what it was tested
against, evidence level and mRNA suitability. M5 is the cocktail: the same
engine run over the Big Four, with the combinations the evidence supports.


