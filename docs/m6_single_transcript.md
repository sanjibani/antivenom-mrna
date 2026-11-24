# M6: two payloads, one transcript

> In-silico designs. No therapeutic efficacy demonstrated. Not medical advice.

## 1. The design question M5 left open

M5 built one construct per payload and said why: a single transcript carrying several antibodies needs a 2A peptide or an internal ribosome entry site, and neither had been sourced to the standard this repository holds sequences to. The attempt to source a 2A peptide failed on every endpoint, and the decision that killed the route became a rule: no sequence enters this repository unless it has been located in the database that publishes it and recorded with a checksum.

The rule held, and the route is open now. The porcine teschovirus-1 2A peptide appears, identical over all 19 residues, in the expression constructs of two deposited structures from two different groups: PDB entry 31FW, entity 4, residues 51 to 69, and PDB entry 8E99, entity 2, residues 4 to 22, in both cases preceded by the glycine-serine-glycine spacer the characterisation recommends. Both entries are committed under `data/reference/raw/` with checksums, and the test suite extracts the peptide from each and checks the two agree with the engine's constant. The characterisation itself is Kim and colleagues, 2011.

## 2. What the skip costs, and which payload pays

A 2A element is not a proteolytic cleavage site. The ribosome simply does not make the peptide bond between the last glycine of the element and the proline that follows it: translation pauses, the upstream product is released with the element still attached, and the downstream product begins with that proline. Nothing is cut, so nothing is clean.

The upstream product carries a 21-residue C-terminal scar: the spacer and the 2A peptide up to the skip glycine. The published binding data for TPL1158_01_C09 belongs to the 123-residue nanobody, not to the 144-residue scarred version this transcript encodes, and the report says so rather than quietly assuming the scar is benign.

The downstream product starts with one proline, and when the downstream payload is secreted that proline sits in front of a signal peptide, which signal peptidase removes together with the proline. The mature downstream product is therefore the published molecule exactly.

That asymmetry is the design rule, and it decides the order: the payload whose published data matters more goes downstream. E113 has reached L4 on the validation ladder, in vivo neutralisation, against L2, binding only, for the nanobody, so E113 goes downstream and the scar goes on the nanobody whose data was binding-only to begin with.

## 3. The construct

The construct is 1722 nt:

| Element | Description | Position |
| --- | --- | --- |
| cap | A modified 5'-cap1 structure (m7G+m3'-5'-ppp-5'-Am). The two nucleotides printed for the cap are the 7-methylguanosine at position 1 and the first backbone nucleotide at position 2, joined by a 5'-to-5' triphosphate. | 1-2 |
| 5'-UTR | 5'-untranslated region derived from human alpha-globin RNA with an optimized Kozak sequence. As used in BNT162b2. | 3-54 |
| signal peptide | Human serum albumin signal peptide, which directs the nascent polypeptide into the endoplasmic reticulum for secretion. | 55-108 |
| TPL1158_01_C09 VHH | Codon-optimised camelid VHH TPL1158_01_C09 against alpha-cobratoxin, carrying the GSG-P2A scar at its C-terminus: 21 residues the published binding data does not describe | 109-477 |
| P2A | The porcine teschovirus-1 2A ribosome-skipping peptide, preceded by the glycine-serine-glycine spacer both deposited constructs use | 478-540 |
| signal peptide 2 | The proline the ribosome carries over from the skip, followed by a second human serum albumin signal peptide. Signal peptidase removes both together, so the mature downstream product starts exactly at the published E113 sequence | 541-597 |
| E113 scFv | Codon-optimised human scFv E113 against Echis carinatus phospholipase A2, reproduced residue for residue from the published molecule, followed by a stop codon | 598-1317 |
| 3'-UTR | A 3' untranslated region comprising two sequence elements derived from the amino-terminal enhancer of split (AES) mRNA and the mitochondrial encoded 12S ribosomal RNA, to confer RNA stability and high total protein expression. As used in BNT162b2. | 1318-1612 |
| poly(A) | A 110-nucleotide poly(A)-tail consisting of a stretch of 30 adenosine residues, followed by a 10-nucleotide linker sequence and another 70 adenosine residues. | 1613-1722 |

| Measurement | Value |
| --- | --- |
| Coding region length | 1260 nt |
| GC content of the coding region | 0.6381 |
| CAI against human codon usage | 0.9882 |
| GC range across 30-nucleotide blocks | 0.533 to 0.700 |
| Forbidden motifs in the coding elements | 0 |
| Forbidden motifs in the UTR elements | 2 |
| CAI is optimal subject to the GC band | yes |

## 4. What the engine checked

- **The P2A peptide is the same molecule in two deposited constructs.** The engine reads both raw files and checks each against the constant in `src/mrna/p2a.py`; the test suite does the same against the committed raw downloads.
- **The open reading frame translates to the designed protein.** 420 residues, checked by translation, not by counting codons.
- **The mature upstream product is the published nanobody plus its scar, 144 residues**, checked by translation of the corresponding element slices.
- **The mature downstream product is the published E113 sequence, 239 residues, reproduced exactly.** This is the check that makes the design claim true rather than asserted: the carry-over proline and the second signal peptide are gone, and what remains is the molecule the 2024 paper measured.
- **The coding elements carry no forbidden motif.** The policy enforces the polyadenylation signal out of every coding element, and this transcript carries 0 occurrences there; the 2 the construct does carry are both in the UTRs, where the licensed product carries them too, and the M2 report explains why they are not defects.

## 5. What this transcript is not

- **It is not two published antibodies.** It is one scarred nanobody and one exact scFv, on one message. The scar is stated in the element table, in the FASTA header of the upstream product, and in this paragraph, three times, because a design whose payload differs from the published molecule has to say so where a reader will look.
- **It is not a free cocktail.** Both payloads bind in vitro and the scFv has neutralisation data on its own; nothing published measures the two molecules together, on one message or on two. The M5 rule that a plan is not a molecule applies here with more force: co-expression is not co-formulation, and the ratio of the two products depends on skip efficiency, which the characterisation paper reports as high for P2A but which this project cannot measure.
- **It is not licensed.** TPL1158_01_C09 is named 71 times in WO2025078529A1, published by Danmarks Tekniske Universitet, and the position is unchanged from M5: an open programme built on that sequence is not free to exist, and this construct is a design artifact rather than a freedom-to-operate opinion.

## 6. The third payload, and what it waits for

The refresh added a third transcribable payload, D2, the 122-residue designed nanobody against alpha-cobratoxin whose sequences sit in an openly downloadable preprint supplement. It is committed reference data now, with the extraction pinned to the deposited row by the 9GCN framework prefix. It does not join this transcript, for two reasons the report states rather than resolves: a three-payload transcript has two scars and two ratios where this one has one, and D2 is a designed sequence whose liabilities the preprint reports but whose molecule nobody has neutralised in vivo against whole venom at the level E113 has. The second P2A site a three-payload version needs is one optimiser run away.

## 7. What comes next

The regional-variant strand: the venom variation inside one species, mapped for the Indian spectacled cobra by Laxme and colleagues in 2021, treated as a design input rather than a footnote. The tracker names the two gaps a variant would fill if a sequence ever lands: krait phospholipase A2, and the metalloproteinases of Indian Russell's viper.

## References

- Kim JH, Lee SR, Li LH, et al. 2011. High cleavage efficiency of a 2A peptide derived from porcine teschovirus-1 in human cell lines, zebrafish and mice. PLoS ONE 6(4):e18556. doi:10.1371/journal.pone.0018556.
- Donnelly ML, Hughes LE, Luke G, et al. 2001. The 2A/2B polyprotein 'cleavage' mechanism. J Gen Virol 82:1027-1041.
- Szymczak AL, Workman CJ, Wang Y, et al. 2004. Correction of multi-gene deficiency in vivo using a single 'self-cleaving' 2A peptide-based vector. Nat Biotechnol 22:589-594.
- PDB entries 31FW and 8E99, the two deposited constructs the peptide is read from.
- RCSB Protein Data Bank entries 9GCN and 8IA6, the payloads; Overath et al. 2026, bioRxiv 2026.09.01.748349, the third payload.
