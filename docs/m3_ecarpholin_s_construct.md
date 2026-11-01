# One toxin, one published antibody, one full construct

> In-silico designs. No therapeutic efficacy demonstrated. Not medical advice.

This is the M3 artifact, and it is the first venom construct in this repository. It encodes a published human single-chain variable fragment against a toxin of the saw-scaled viper, *Echis carinatus*, on the cap, UTR and poly(A) scaffold that M2 built and measured against a licensed vaccine.

The claim this document makes is narrow and checkable. Every sequence has a source, every number was computed by `scripts/venom_construct.py` from the committed reference data, and the construct is reproducible byte for byte. The document does not claim that the molecule works, because nobody has expressed it.

## 1. The toxin

| Field | Value |
| --- | --- |
| Species | *Echis carinatus*, the saw-scaled viper, sochurek |
| In India | One of the Big Four, the four species that account for the overwhelming majority of snakebite deaths in the country |
| Toxin | Ecarpholin S, the myotoxic Ser49 phospholipase A2 |
| Sequence | UniProt P48650, 122 residues |
| Structure | PDB 2QHE, 2.0 A resolution |
| Family | Snake venom phospholipase A2 |

Phospholipase A2 is the right place to start, and the source paper says why in numbers rather than adjectives: a proteomic survey of the venom found PLA2 to be the dominant toxin, at 28 per cent of *Echis carinatus* venom by mass. PLA2 is also the toxin family that does the most damage in viper bites. It is myotoxic because the C-terminal hydrophobic loop inserts into the muscle membrane, it is hemolytic because it hydrolyses the lipid bilayer of the red cell, and it is anticoagulant because it cleaves the phospholipids that clotting factors need.

Ecarpholin S sits in the Ser49 subgroup of the venom phospholipases, and that single residue is what makes it a myotoxin rather than an enzyme. Most Ca2-dependent PLA2 enzymes carry an aspartate at position 49, and that aspartate coordinates the calcium ion the catalytic cycle needs. Ecarpholin S has a serine there instead, and the structure shows a calcium-binding loop folded into a conformation that cannot coordinate the ion. The result, in the words of the structural paper, is 'low enzymatic but potent myotoxic activities'. The catalytic machinery is all still there; it simply cannot run.

The molecule is a problem for an mRNA product in a specific way. A snake venom PLA2 is a small, disulfide-rich protein: UniProt annotates seven disulfide bonds for this sequence, which means it has to fold correctly, with the right disulfide pairing, in an oxidising environment. The cytosol is reducing and the wrong place for it. That is the reason the design below carries a signal peptide at all, and it is the reason the payload is an antibody rather than the toxin itself: the antibody is made in the endoplasmic reticulum where the disulfides can form.

Two things in the UniProt record connect the toxin to the antibody. The first is the region annotated as 'important for membrane-damaging activities in eukaryotes and bacteria', residues 105 to 117, the C-terminal hydrophobic loop that does the actual damage. The second is the FUNCTION line, which says this is a 'snake venom phospholipase A2 homolog that lacks enzymatic activity. Shows high myotoxin activities and displays edema-inducing activities'. The source paper's structural model of the E113-PLA2 complex puts the antibody's CDR regions against exactly that C-terminal myotoxic loop and against the interfacial region of the toxin, which is the part of a PLA2 that docks into a membrane. Two independent sources, a sequence database and a paper, put the antibody's epitope on the part of the toxin that does the damage.

The target sequence, as committed:

```
SVVELGKMIIQETGKSPFPSYTSYGCFCGGGERGPPLDATDRCCLAHSCCYDTLPDCSPKTDRYKYKRENGEIICENSTSCKKRICECDKAVAVCLRKNLNTYNKKYTYYPNFWCKGDIEKC
```

## 2. The antibody

| Field | Value |
| --- | --- |
| Clone | E113 |
| Format | Human single-chain variable fragment, one polypeptide |
| Sequence | PDB 8IA6 entity 1, 2.75 A crystal structure |
| Source | Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2024. Mol Immunol 165:55-67 |
| Patent | No patent on this clone located in a public search on 2026-11-01 |

The payload is split into its three parts:

| Part | Length | Sequence |
| --- | --- | --- |
| VH | 116 aa | `EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSSISDNGQTTTYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKTPYLFDYWGQGTLVTVSS` |
| linker | 16 aa | `GGGGSGGGGSGGGGST` |
| VL | 107 aa | `DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYSASHLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQHRRPPPTFGQGTKVEIK` |

The heavy chain is IGHV3-23, the most used human immunoglobulin heavy chain variable gene. The light chain is an IGKV1 variable domain joined to the IGKC constant region, which is where the `FGQGTKVEIK` at the end comes from.

### What the deposited construct is, and what was removed

The Protein Data Bank entry holds a crystallisation construct, not a therapeutic molecule. It has a vector residue at the N-terminus and a hexahistidine expression tag at the C-terminus. The crystal locates residues 2 to 239 of the deposited chain: residue 1, the vector alanine, was never located, residues 118 to 133, the flexible linker, are disordered in every copy, and the expression tag was never located either.

The design encodes the 239-residue antibody and neither artefact. The removal is done by the reference-data script, not by hand, and the script refuses to run if the deposited construct does not have the expected ends.

## 3. The construct

| Element | Description | Position |
| --- | --- | --- |
| cap | A modified 5'-cap1 structure (m7G+m3'-5'-ppp-5'-Am). The two nucleotides printed for the cap are the 7-methylguanosine at position 1 and the first backbone nucleotide at position 2, joined by a 5'-to-5' triphosphate. | 1-2 |
| 5'-UTR | 5'-untranslated region derived from human alpha-globin RNA with an optimized Kozak sequence. As used in BNT162b2. | 3-54 |
| sig | Human serum albumin signal peptide, which directs the nascent polypeptide into the endoplasmic reticulum for secretion. | 55-108 |
| E113 scFv | Codon-optimised single-chain variable fragment E113 against Echis carinatus phospholipase A2, preceded by the human serum albumin signal peptide, followed by a stop codon at 826 to 828. | 109-828 |
| 3'-UTR | A 3' untranslated region comprising two sequence elements derived from the amino-terminal enhancer of split (AES) mRNA and the mitochondrial encoded 12S ribosomal RNA, to confer RNA stability and high total protein expression. As used in BNT162b2. | 829-1123 |
| poly(A) | A 110-nucleotide poly(A)-tail consisting of a stretch of 30 adenosine residues, followed by a 10-nucleotide linker sequence and another 70 adenosine residues. | 1124-1233 |

| Property | Value |
| --- | --- |
| Total length | 1233 nt |
| Translated coding sequence | 774 nt |
| Encoded protein | 257 aa |
| Payload | 239 aa |
| GC content of the coding region | 0.6381 |
| CAI against human codon usage | 0.9854 |
| GC per 30 nt block | 0.533 to 0.700 |
| Forbidden motifs in the coding region | 0 |
| Upstream ATG in the 5' UTR | 0 |
| AAUAAA in the whole construct | 2 (3'-UTR, 5'-UTR) |

The design decisions, each with the reason for it:

| Decision | Why |
| --- | --- |
| Human serum albumin signal peptide | One of the two most commonly used signal peptides in mRNA vaccine designs, and the address label that gets a disulfide-rich antibody out of the cell |
| Human alpha-globin 5' UTR | The element Comirnaty uses, so that the construct reads line by line against a licensed product |
| AES plus mitochondrial 12S rRNA 3' UTR | The same, and the 3' UTR that a screen of human UTR fragments found to be the most expression-augmenting |
| A110 template-encoded poly(A) | 30 adenosines, a 10-nucleotide linker, 70 adenosines. The linker is what a ribosome that misreads past the stop codon runs into instead of a longer poly-lysine tract |
| Single TGA stop codon | The minimal choice. Comirnaty prints two, which does not change what is translated |
| Codon optimisation subject to the GC band | The M1 result: the CAI is optimal subject to a band that both licensed products satisfy |

## 4. What the engine checked

Four checks, all of them mechanical, all of them reproduced by the test suite.

**The construct encodes what it claims.** The coding region translates back to 257 residues: the 18-residue signal peptide followed by the 239-residue antibody, with no substitutions.

**The initiation context is the right one.** The alpha-globin 5' UTR ends in GCCACC, so the start codon sits in the Kozak context that the same construct uses in a licensed product.

**There is no upstream start codon.** A start codon in the 5' UTR would compete with the intended one, and the ribosome usually takes the first it meets. There are none.

**The coding region holds every GC block inside the band and carries no forbidden motif.** The polyadenylation signal that Comirnaty carries twice in its untranslated regions is absent from the coding region, which is where it would matter.

## 5. What the source paper actually shows

The evidence table is quoted rather than paraphrased, because the difference between the two is the difference between an honest design document and a promotional one.

| | Reported |
| --- | --- |
| Binding | PLA2 was 'prominently captured by the phage-anchored scFv antibodies', and was the dominant toxin of both venoms by proteomics, at 24 per cent of Naja naja venom and 28 per cent of Echis carinatus venom. |
| Target | E113 and E10, the two anti-Echis carinatus clones, both interact with the basic PLA2 of the venom. |
| In vivo | E113 gave 'significant prolongation of survival time' in Echis carinatus venom-challenged mice. The equal-dose combination of E113 and E10 gave 33 per cent survival, against 16 per cent for the individual scFvs. |
| Mechanism | E113 inhibited venom-induced hemolysis by 50 per cent, and a structural model of the E113-PLA2 complex placed the CDR regions against the reported C-terminal myotoxic loop and interfacial region of PLA2. |
| Format | A single-chain variable fragment from a human antibody library, so one polypeptide, one open reading frame, one construct. |

## 6. What this construct is not

An scFv is the weakest possible starting point for an antivenom, and the reasons are worth stating before anyone builds it.

**It is monovalent.** One binding site means one toxin per antibody. The published data on this clone is about survival in mice challenged with venom, which is a whole-venom endpoint, and the combination of E113 with a second clone did better than either alone. A bivalent format would be the obvious improvement, and it would need two chains or a single-chain Fab.

**It is small and it is cleared.** An scFv of this size is roughly 26 kDa, has a serum half-life in the range of hours, and is cleared by the kidney. An IgG or an Fc fusion would last weeks.

**The timing is wrong for an acute bite.** Intravenous antivenom works within minutes. An mRNA-encoded antibody is not detectable for one to two hours and peaks somewhere between 12 and 24 hours. This construct is a prophylactic candidate, and the repository says so at the top of its README.

**Nothing has been expressed.** No cell has made this protein. The next rung of the validation ladder is L1, expression in cells, and that belongs to a wet-lab partner. This project does not have one and does not intend to.

**The antigen-binding claim is about the antibody, not about the mRNA.** The published binding and neutralisation data belongs to E113 as a protein made in bacteria. An scFv made in a mammalian cell from an mRNA template is the same sequence, but it is a different molecule until someone shows it is not: glycosylation, disulfide pairing and folding all differ between a bacterial and a mammalian expression system, and for a disulfide-rich antibody that is not a detail.

**The patent position is a negative search result, not a clearance.** No patent covering the E113 clone was located in a public patent search on 2026-11-01. That is not a freedom-to-operate opinion, and milestone M4 exists to record patent status properly across the whole candidate set.

**Venom varies regionally.** The same species does not have the same venom everywhere, and a construct aimed at India has to be aimed at Indian venom. The east-coast and the west-coast *Echis carinatus* are not the same animal.

## 7. What comes next

M4 is the tracker: species, toxin family, target, antibody, format, sequence source, patent status, what it was tested against, evidence level and mRNA suitability, for every candidate in the field. M5 is the cocktail: the same engine, run over the Big Four, with the combinations that the evidence supports.

## References

- Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2024. Recombinant human scFv antibody fragments against phospholipase A2 from Naja naja and Echis carinatus snake venoms: in vivo neutralization and mechanistic insights. Mol Immunol 165:55-67. doi:10.1016/j.molimm.2023.12.006.
- Protein Data Bank entry 8IA6, crystal structure of scFv antibody against phospholipase A2 of Echis carinatus venom, 2.75 A, deposited by Kumar A, Madni ZK, Salunke DM.
- Protein Data Bank entry 2QHE, crystal structure of the phospholipase A2 of Echis carinatus, the toxin E113 binds.
- UniProt P48650, ecarpholin S, basic phospholipase A2 homolog from Echis carinatus.
- UniProt P02768, human serum albumin precursor, the source of the signal peptide.
- Zhou X, Tan TC, Valiyaveettil S, Go ML, Kini RM, Velazquez-Campoy A, Sivaraman J. 2008. Structural characterization of myotoxic ecarpholin S from Echis carinatus venom. Biophys J 95:3366-3380. doi:10.1529/biophysj.107.117747.
- Lomonte B, Rangel J. 2012. Snake venom Lys49 myotoxins: from phospholipases A2 to non-enzymatic membrane disruptors. Toxicon 60:520-530.
- Orlandini von Niessen AG, et al. 2019. Improving mRNA-based therapeutic gene delivery by expression-augmenting 3' UTRs identified by cellular library screening. Mol Ther 27(4):824-836.
- WHO International Nonproprietary Names Programme, document 11889, Comirnaty (BNT162b2), February 2021.

_In-silico designs. No therapeutic efficacy demonstrated. Not medical advice._
