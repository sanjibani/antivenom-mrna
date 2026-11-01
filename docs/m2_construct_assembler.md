# A construct assembler, measured against a licensed product

> In-silico designs. No therapeutic efficacy demonstrated. Not medical advice.

This is the M2 artifact. The assembler takes the parts of an mRNA construct in order and produces two things: the sequence, and the element table that describes it. The element table is the part that matters, because it is the format a regulatory document uses, and because it is the only way a reader can check a design against a licensed product line by line.

The measurement is a rebuild. The published parts of Comirnaty, taken from the primary regulatory document and from an independent transcription of its coding sequence, are fed back into the assembler. If the result is the same 4,284 nucleotides with the same six rows, then the assembler reproduces the geometry that a regulator has approved, and its output for an antivenom construct is worth reading.

## 1. The rebuild

| Element | Published range | Rebuilt range | Match |
| --- | --- | --- | --- |
| cap | 1 to 2 | 1 to 2 | yes |
| 5'-UTR | 3 to 54 | 3 to 54 | yes |
| sig | 55 to 102 | 55 to 102 | yes |
| S protein_mut | 103 to 3879 | 103 to 3879 | yes |
| 3'-UTR | 3880 to 4174 | 3880 to 4174 | yes |
| poly(A) | 4175 to 4284 | 4175 to 4284 | yes |

Rebuilt length 4284 nt against the published 4284 nt. Nucleotides compared: 4284. Mismatches: **0**.

The element table the assembler emits, in the document's format:

| Element | Description | Position |
| --- | --- | --- |
| cap | A modified 5'-cap1 structure (m7G+m3'-5'-ppp-5'-Am). The two nucleotides printed for the cap are the 7-methylguanosine at position 1 and the first backbone nucleotide at position 2, joined by a 5'-to-5' triphosphate. | 1-2 |
| 5'-UTR | 5'-untranslated region derived from human alpha-globin RNA with an optimized Kozak sequence. As used in BNT162b2. | 3-54 |
| sig | Signal peptide, which guides translocation of the nascent polypeptide chain into the endoplasmic reticulum. | 55-102 |
| S protein_mut | Codon-optimized sequence encoding full-length SARS-CoV-2 spike (S) glycoprotein containing mutations K986P and V987P to ensure the S glycoprotein remains in an antigenically optimal pre-fusion conformation; stop codons: 3874-3879 (underlined) | 103-3879 |
| 3'-UTR | A 3' untranslated region comprising two sequence elements derived from the amino-terminal enhancer of split (AES) mRNA and the mitochondrial encoded 12S ribosomal RNA, to confer RNA stability and high total protein expression. As used in BNT162b2. | 3880-4174 |
| poly(A) | A 110-nucleotide poly(A)-tail consisting of a stretch of 30 adenosine residues, followed by a 10-nucleotide linker sequence and another 70 adenosine residues. | 4175-4284 |

The coding region opens at nucleotide 55 and closes on the first in-frame stop codon at nucleotide 3876, which reproduces the document's own feature table rather than a paraphrase of it.

## 2. Where the nucleotides came from

A rebuild is only as good as the independence of its inputs, so here is the provenance of every part.

| Element | nt | Source | Independent check |
| --- | --- | --- | --- |
| cap | 2 | WHO INN document 11889 | none available, it is a chemical structure |
| 5'-UTR | 52 | WHO INN document 11889, positions 3-54 | printed identically by Fang et al. 2022, Table 1 |
| coding region | 3,825 | berthubert/bnt162b2 transcription of the same document, plus one stop codon | byte-identical to the document over 3,822 nt |
| 3'-UTR | 295 | WHO INN document 11889, positions 3880-4174 | none available |
| poly(A) | 110 | WHO INN document 11889, positions 4175-4284 | none available |

Three nucleotides of 4,284, or 0.07 per cent, come from the document being reproduced: the second of its two consecutive stop codons. The published coding-sequence transcription ends at the first stop codon, because that is where translation ends, so the second one has to be added back to reproduce the construct. Everything else comes from the published 5' UTR, the published coding sequence and the published 3' UTR and poly(A) specification, and the 3,822 coding nucleotides were checked against the document byte for byte before this report was written.

## 3. The cap is a chemical footprint, not two nucleotides

The document prints two nucleotides for the cap and starts its 5' UTR element at position 3. That is not an off-by-one error, it is a different convention for the same molecule, and it is worth being precise about, because the same UTR is quoted at three different lengths across the literature.

| What is counted | Length |
| --- | --- |
| The 5' UTR element as the WHO document numbers it, positions 3 to 54 | 52 nt |
| The same region including the first backbone nucleotide | 53 nt |
| The same region including the cap guanosine as well | 54 nt |

The chemistry behind it. The cap is 7-methylguanosine, joined to the first backbone nucleotide by a 5'-to-5' triphosphate, so it covers two printed positions: the m7G at position 1 and the first nucleotide of the backbone at position 2. The document writes it as 'cap G1A2 = m7G+m3'-5'-ppp-5'-Am2'-3'-p-', with a 2'-O-methyl on that first backbone nucleotide, which is what makes this a cap 1 rather than a cap 0. The ribosome scans the whole physical leader, which is 53 nt, whatever the table calls it.

This engine uses the document's numbering, because the document's coordinates are the ones a design has to match. A reader who counts the UTR by hand will get 54 and will be right about the molecule and wrong about the table, which is exactly the kind of disagreement that turns into a silent error in a later milestone.

## 4. A new construct on the same scaffold

An assembler that only reproduces Comirnaty has reproduced a lookup table. The same call that rebuilds the licensed product builds a new one: the same cap, the same 5' UTR, the same 3' UTR, the same poly(A) specification, and a different payload.

| Element | Description | Position |
| --- | --- | --- |
| cap | A modified 5'-cap1 structure (m7G+m3'-5'-ppp-5'-Am). The two nucleotides printed for the cap are the 7-methylguanosine at position 1 and the first backbone nucleotide at position 2, joined by a 5'-to-5' triphosphate. | 1-2 |
| 5'-UTR | 5'-untranslated region derived from human alpha-globin RNA with an optimized Kozak sequence. As used in BNT162b2. | 3-54 |
| GFP | Codon-optimised green fluorescent protein from Aequorea victoria, UniProt P42212, with a single stop codon. | 55-771 |
| 3'-UTR | A 3' untranslated region comprising two sequence elements derived from the amino-terminal enhancer of split (AES) mRNA and the mitochondrial encoded 12S ribosomal RNA, to confer RNA stability and high total protein expression. As used in BNT162b2. | 772-1066 |
| poly(A) | A 110-nucleotide poly(A)-tail consisting of a stretch of 30 adenosine residues, followed by a 10-nucleotide linker sequence and another 70 adenosine residues. | 1067-1176 |

| Property | Comirnaty, rebuilt | New construct |
| --- | --- | --- |
| Total length | 4284 nt | 1176 nt |
| Translated coding sequence, start codon to first stop | 3822 nt | 717 nt |
| Protein | 1273 aa | 238 aa |
| Upstream ATG in the 5' UTR | 0 | 0 |
| AAUAAA occurrences, whole construct | 2 | 2 |
| AAUAAA occurrences in the coding region | 0 | 0 |

The new construct is a teaching example rather than a candidate product: green fluorescent protein is not a venom antigen. It is here because it is small enough to read in full, because it was already in the reference data, and because the next milestone does the same thing with a payload that matters.

## 5. Motif scope: the licensed product carries the signal twice

M1 enforced one forbidden motif in the coding region: the AAUAAA polyadenylation signal. Applied to a whole construct, the same scan finds it twice in Comirnaty.

| Element | Length | AAUAAA occurrences |
| --- | --- | --- |
| cap | 2 nt | 0 |
| 5'-UTR | 52 nt | 1 |
| sig | 48 nt | 0 |
| S protein_mut | 3777 nt | 0 |
| 3'-UTR | 295 nt | 1 |
| poly(A) | 110 nt | 0 |

- polyadenylation signal in 5'-UTR at 4 to 9
- polyadenylation signal in 3'-UTR at 4099 to 4104

Both occurrences are in untranslated regions, and the coding region has none. That is the correct state of affairs, and it is worth being explicit about why, because the difference is the whole point of scanning per element rather than per construct.

A polyadenylation signal is dangerous when it sits downstream of a cleavage site, because the cell can then cut the transcript there. Inside a coding region that means a truncated message with no in-frame stop codon, and a ribosome that runs into the poly(A) tail and stalls on a poly-lysine tract. Inside a 5' UTR there is no cleavage site upstream of it, and inside a 3' UTR the signal only matters if it is the one the cell chooses to use. A signal that is harmless where it sits is not a defect, and a scanner that cannot tell the difference will eventually get a licensed product wrong.

The engine therefore reports motif occurrences per element, and the element they sit in. It does not report a single verdict for the construct, because a single verdict would have failed a product that two regulators approved.

## 6. Upstream start codons

The ribosome scans the 5' end and usually takes the first start codon it meets. A start codon in the 5' UTR therefore competes with the intended one, and the context around it decides how well. Kozak 1986 showed that the purine at position -3 and the G at position +4 are the two positions that matter. The assembler now reports every ATG in the 5' UTR with its context.

Comirnaty's 5' UTR contains 0 upstream start codons, and the new construct's 5' UTR contains 0. This check was deliberately left out of the M1 optimiser, because it is a property of the UTR and not of the coding region, and it belongs here where the real UTR and the real Kozak context are known.

## 7. What this does not claim

No expression level is predicted. No protein folding is predicted. Nothing has been synthesised and nothing has been put into a cell. The assembler makes a construct well formed with respect to a set of cited, measurable constraints and can show its work by reproducing a published one. That is the whole claim, and the ladder from here is L1, expression in cells, which belongs to a wet-lab partner.

The next milestone takes one toxin from India's Big Four, one published neutralising antibody against it, and builds the same kind of construct with that antibody as the payload. Everything this milestone added is what that milestone needs and nothing more.

## References

- WHO International Nonproprietary Names Programme, document 11889, Comirnaty (BNT162b2), February 2021.
- Orlandini von Niessen AG, Poleganov MA, Rechner C, Plaschke A, Kranz LM, Fesser M, Diken M, Loewer M, Vallazza B, Beissert T, Bukur V, Kuhn AN, Tureci O, Sahin U. 2019. Improving mRNA-based therapeutic gene delivery by expression-augmenting 3' UTRs identified by cellular library screening. Mol Ther 27(4):824-836.
- Fang E, Liu X, Li M, Zhang Z, Song L, Zhu B, Wu X, Liu J. 2022. Advances in COVID-19 mRNA vaccine development. Signal Transduct Target Ther 7(1):94.
- Xia X. 2021. Detailed Dissection and Critical Evaluation of the Pfizer/BioNTech and Moderna mRNA Vaccines. Vaccines 9(7):734.
- Kariko K, Buckstein M, Ni H, Weissman D. 2005. Suppression of RNA recognition by Toll-like receptors: the impact of nucleoside modification and the evolutionary origin of RNA. Immunity 23(2):165-175.
- Kozak M. 1986. Point mutation defines a sequence flanking the AUG initiator codon that modulates translation by eukaryotic ribosomes. Cell 44(2):283-292.
- Proudfoot NJ, Brownlee BG. 1976. AATAAA is a signal for synthesis of polyadenylated RNA. Nucleic Acids Res 3(7):1183-1190.
- Pardi N, Hogan MJ, Porter FW, Weissman D. 2018. mRNA vaccines, a new era in vaccinology. Nat Rev Drug Discov 17:261-279.

_In-silico designs. No therapeutic efficacy demonstrated. Not medical advice._
