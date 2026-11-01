# M5: the Big Four cocktail, and what the published record cannot cover

> In-silico designs. No therapeutic efficacy demonstrated. Not medical advice.

This is the M5 artifact, and most of it is a negative result. The question was: if you had to build one antibody cocktail covering India's Big Four snakes, which payloads would it contain, and what would stop you? The answer is computed by `scripts/cocktail_design.py` from the M4 candidate table, not asserted, and the short version is that the published record cannot supply the cocktail.

## 1. What the planner was asked

A cell is a species and a toxin family. The planner treats cocktail design as a set cover: it finds the smallest set of antibodies that covers as much of the matrix as the published record allows, and it searches exhaustively rather than greedily, so the answer is exact. 4 tiers tighten the question one constraint at a time.

The species are the Big Four plus the monocled cobra. The monocled cobra, *Naja kaouthia*, is not one of the Big Four, but it is the other major Indian cobra, it is one of the species Indian polyvalent antivenom is raised against, and a cocktail aimed at Indian snakebite that ignored it would be incomplete in a way that matters.

| Tier | The question | Payloads | Cells covered |
| --- | --- | --- | --- |
| `published` | Is there a published candidate with binding evidence at all? | 8 | 9 of 14 |
| `sequence-available` | Can its sequence be located in a public database? | 5 | 6 of 14 |
| `single-polypeptide` | Is it one polypeptide, so that one transcript can carry it? | 2 | 2 of 14 |
| `licence-unblocked` | Has a public patent search found no application covering it? | 1 | 1 of 14 |

The attrition down that table is the finding. Four tiers, and the number of cells the plan can actually cover falls from 9 to 1.

## 2. The coverage matrix

| Species | Toxin family | `published` | `sequence-available` | `single-polypeptide` | `licence-unblocked` |
| --- | --- | --- | --- | --- | --- |
| **Naja naja** | 3FTx-L | yes | yes | no | no |
| **Naja naja** | 3FTx-S | yes | yes | no | no |
| **Naja naja** | PLA2 | yes | no | no | no |
| **Bungarus caeruleus** | 3FTx-L | yes | yes | no | no |
| **Bungarus caeruleus** | PLA2 | no | no | no | no |
| **Daboia russelii** | LAAO | no | no | no | no |
| **Daboia russelii** | PLA2 | yes | no | no | no |
| **Daboia russelii** | SVMP | no | no | no | no |
| **Echis carinatus** | LAAO | no | no | no | no |
| **Echis carinatus** | PLA2 | yes | yes | yes | yes |
| **Echis carinatus** | SVMP | yes | yes | no | no |
| Naja kaouthia | 3FTx-L | yes | yes | yes | no |
| Naja kaouthia | 3FTx-S | no | no | no | no |
| Naja kaouthia | PLA2 | yes | no | no | no |

The Big Four species are in bold. The `licence-unblocked` column is a floor rather than a clearance: a row reaches it only when a public patent search was actually performed and found nothing. Of the 22 rows in the M4 table, 9 have had a search performed at all, and 2 of those came back with nothing: E113, E10.

### What each tier drops

**`published`**: N194, Centi-LNX-D09, Centi-SNX-B03, 2554_01_D11, TPL0637_01_A07, MKSS1, E113, H11.

**`sequence-available`**: Centi-LNX-D09, Centi-SNX-B03, TPL1158_01_C09, E113, H11.

**`single-polypeptide`**: TPL1158_01_C09, E113.

**`licence-unblocked`**: E113.

Read down the list of payloads and the order is the story. The published record has candidates for most of the matrix. The public sequence databases do not. The transcripts that can be built from what is public are 2.

## 3. The cells that no tier covers

Even at the most permissive tier, with every published candidate allowed, 5 cells have no candidate at all:

- **Bungarus caeruleus / PLA2** (Big Four)
- **Daboia russelii / SVMP** (Big Four)
- **Daboia russelii / LAAO** (Big Four)
- **Echis carinatus / LAAO** (Big Four)
- **Naja kaouthia / 3FTx-S**

Of these, 3 are viper toxin families that nobody has taken on: Daboia russelii / SVMP, Daboia russelii / LAAO, Echis carinatus / LAAO. The 2 viper species of the Big Four between them carry every toxin class that makes a viper bite dangerous, and not one of those classes has a published antibody candidate in this table.

The other 2 are elapid cells: Bungarus caeruleus / PLA2, Naja kaouthia / 3FTx-S. The krait phospholipase A2 cell is the one to start from, because phospholipase A2 is the dominant fraction of common krait venom by mass, and no published candidate targets it. Krait venom is the most under-covered of the Big Four.

### The cells the planner will not target

3 cells are empty on purpose. The Kunitz-type serine protease inhibitors are the anticoagulant fraction of viper venom: they bind the clotting factors directly, so an antibody against them would protect the clotting cascade in a patient whose problem is that the cascade is being consumed. A planner that reported every empty cell as an opportunity would be pointing this project at a design that would make the patient worse off.

- **KSPI**: the anticoagulant fraction of viper venom; an antibody against it would make the coagulopathy worse rather than better.

## 4. The 2 constructs that can be built

At the `single-polypeptide` tier there are 2 payloads, and every one of them can be transcribed from a published sequence as it stands, with no engineering step in between. The others are Fab fragments or two-chain antibodies: converting them into one transcript needs a single-chain format that does not exist in the published data, which means the neutralisation data would no longer belong to the molecule being transcribed.

### TPL1158_01_C09

| Field | Value |
| --- | --- |
| Target | alpha-cobratoxin |
| Species | Naja kaouthia |
| Format | VHH, one polypeptide |
| Payload length | 123 aa |
| Total construct length | 885 nt |
| Translated coding sequence | 426 nt |
| GC content of the coding region | 0.6217 |
| CAI against human codon usage | 0.9910 |
| GC per 30 nt block | 0.533 to 0.700 |
| Forbidden motifs in the coding region | 0 |
| Upstream ATG in the 5' UTR | 0 |
| AAUAAA in the whole construct | 2 (3'-UTR, 5'-UTR) |
| Source | PDB entry 9GCN, entity 1: a co-crystal of the camelid single-domain antibody TPL1158_01_C09 with alpha-cobratoxin from Naja kaouthia. Released 2025-08-06. |

| Element | Description | Position |
| --- | --- | --- |
| cap | A modified 5'-cap1 structure (m7G+m3'-5'-ppp-5'-Am). The two nucleotides printed for the cap are the 7-methylguanosine at position 1 and the first backbone nucleotide at position 2, joined by a 5'-to-5' triphosphate. | 1-2 |
| 5'-UTR | 5'-untranslated region derived from human alpha-globin RNA with an optimized Kozak sequence. As used in BNT162b2. | 3-54 |
| sig | Human serum albumin signal peptide, which directs the nascent polypeptide into the endoplasmic reticulum for secretion. | 55-108 |
| TPL1158_01_C09 payload | Codon-optimised VHH TPL1158_01_C09 against alpha-cobratoxin, preceded by the human serum albumin signal peptide, followed by a stop codon at 478 to 480. | 109-480 |
| 3'-UTR | A 3' untranslated region comprising two sequence elements derived from the amino-terminal enhancer of split (AES) mRNA and the mitochondrial encoded 12S ribosomal RNA, to confer RNA stability and high total protein expression. As used in BNT162b2. | 481-775 |
| poly(A) | A 110-nucleotide poly(A)-tail consisting of a stretch of 30 adenosine residues, followed by a 10-nucleotide linker sequence and another 70 adenosine residues. | 776-885 |

### E113

| Field | Value |
| --- | --- |
| Target | ecarpholin S, the Ser49 myotoxic phospholipase A2 |
| Species | Echis carinatus |
| Format | scFv, one polypeptide |
| Payload length | 239 aa |
| Total construct length | 1233 nt |
| Translated coding sequence | 774 nt |
| GC content of the coding region | 0.6381 |
| CAI against human codon usage | 0.9854 |
| GC per 30 nt block | 0.533 to 0.700 |
| Forbidden motifs in the coding region | 0 |
| Upstream ATG in the 5' UTR | 0 |
| AAUAAA in the whole construct | 2 (3'-UTR, 5'-UTR) |
| Source | Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2024. Mol Immunol 165:55-67. doi:10.1016/j.molimm.2023.12.006. PDB entry 8IA6, entity 1, with the crystallisation vector residue and the hexahistidine tag removed by scripts/prepare_reference_data.py. |

| Element | Description | Position |
| --- | --- | --- |
| cap | A modified 5'-cap1 structure (m7G+m3'-5'-ppp-5'-Am). The two nucleotides printed for the cap are the 7-methylguanosine at position 1 and the first backbone nucleotide at position 2, joined by a 5'-to-5' triphosphate. | 1-2 |
| 5'-UTR | 5'-untranslated region derived from human alpha-globin RNA with an optimized Kozak sequence. As used in BNT162b2. | 3-54 |
| sig | Human serum albumin signal peptide, which directs the nascent polypeptide into the endoplasmic reticulum for secretion. | 55-108 |
| E113 payload | Codon-optimised scFv E113 against ecarpholin S, the Ser49 myotoxic phospholipase A2, preceded by the human serum albumin signal peptide, followed by a stop codon at 826 to 828. | 109-828 |
| 3'-UTR | A 3' untranslated region comprising two sequence elements derived from the amino-terminal enhancer of split (AES) mRNA and the mitochondrial encoded 12S ribosomal RNA, to confer RNA stability and high total protein expression. As used in BNT162b2. | 829-1123 |
| poly(A) | A 110-nucleotide poly(A)-tail consisting of a stretch of 30 adenosine residues, followed by a 10-nucleotide linker sequence and another 70 adenosine residues. | 1124-1233 |

## 5. The licence position, which is the real blocker

These 2 payloads are the only ones that can be transcribed, and a public patent search finds an application covering 1 of them.

- **TPL1158_01_C09** is named 71 times in WO2025078529A1, *Nanobodies capable of cross-binding to several different elapid three-finger toxins*, published by Danmarks Tekniske Universitet.
- **E113** is the clone this repository already built in M3, and no application covering it was located.

This is not a clearance, and it is not a freedom-to-operate opinion. It is a record of what the public patent record says, on a date, about the applications this project searched. It does mean that an open, reproducible mRNA programme built on the published sequences is not free to exist, and that the honest obstacle to an open Big Four cocktail is as much legal as scientific.

## 6. What this cocktail is not

**It covers 1 Big Four species out of 4.** 2 of 14 cells, and the only Big Four species in it is Echis carinatus. It is the first honest attempt, not a finished product, and the report above says exactly which cells are missing.

**Indian Russell's viper is not in it at all.** The one recombinant candidate against the dominant toxin of that venom, the camelid VHH MKSS1, has no publicly locatable sequence, so at the `sequence-available` tier the whole species drops out of the plan. That is the single most actionable line in this document.

**Nothing has been expressed.** No cell has made any of these proteins. The next rung of the validation ladder is L1, and it belongs to a wet-lab partner this project does not have.

**The timing is wrong for an acute bite, exactly as in M3.** An mRNA-encoded antibody is not detectable for one to two hours and peaks somewhere between 12 and 24 hours. Intravenous antivenom works within minutes. This is a prophylactic and a regional-manufacturing proposition.

**A cocktail is not a single molecule.** Each payload is its own transcript here, because a single transcript carrying several antibodies needs a 2A peptide or an internal ribosome entry site, and neither has been sourced for this repository to the standard the rest of the sequences meet. That is the first thing M6 should do, and it is a real engineering step rather than a formatting one.

**The payloads are monovalent.** An scFv and a VHH each have one binding site, so each covers one toxin per antibody molecule. The published combination data for the E113 clone shows that a second clone helps, which is an argument for a cocktail and not for a single construct.

## 7. What comes next

M6 is regional variants, partner outreach and design iteration. The tracker has now said what the work should be: the krait phospholipase A2 cell, the Russell's viper metalloproteinase cell, and a licence conversation with the two institutions whose applications cover the sequences that exist. None of those is a design problem, which is why this milestone ends with a map rather than with more constructs.

## References

- Samanta A, Rudresha GV, Ahmadi S, et al. 2026. Oligoclonal nanobody-based recombinant antivenom protects mice challenged with venom from cobras and king cobras from India. Sci Transl Med 18:866. doi:10.1126/scitranslmed.aed4290.
- Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2024. Mol Immunol 165:55-67. doi:10.1016/j.molimm.2023.12.006.
- Glanville J, Kwong K, et al. 2025. Broad-spectrum neutralisation of elapid venoms by a minimal antibody combination. Cell. Patent US20240294618A1.
- WO2025078529A1. Nanobodies capable of cross-binding to several different elapid three-finger toxins. Danmarks Tekniske Universitet.
- WO2023209212A1. Neutralisation of alpha-neurotoxins using human recombinant IgG antibodies. Danmarks Tekniske Universitet.
- Szymczak AL, Workman CJ, Wang Y, et al. 2004. Correction of multi-gene deficiency in vivo using a single 'self-cleaving' 2A peptide-based retroviral vector. Nat Biotechnol 22:589-594.
- Protein Data Bank entries 9GCN, 8DA1, 8D9Z, 8V13, 9CLP, 9HUO, 9FYT and 9HXO.

_In-silico designs. No therapeutic efficacy demonstrated. Not medical advice._
