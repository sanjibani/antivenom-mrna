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

## Next

The open items for milestone M1 are the codon optimiser, which must hold the
GC band while raising CAI, and the UTR library, whose first entries will be
the human alpha-globin 5' UTR and the AES plus 12S rRNA 3' UTR that Comirnaty
uses, both with their citations.
