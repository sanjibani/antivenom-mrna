# A codon optimiser, measured against two approved products

> In-silico designs. No therapeutic efficacy demonstrated. Not medical advice.

This is the M1 artifact. The optimiser takes a protein and a human codon usage table and returns a coding sequence that holds its GC content inside a band, contains no forbidden motif, and is optimal for the Codon Adaptation Index subject to those constraints.

The interesting measurement is not what the optimiser produces. It is whether the optimiser, given only the protein, arrives anywhere near the codon choices that BioNTech and Moderna published. They chose those codons with more information than we have: budgets, manufacturing constraints, immunogenicity data and years of iteration. If our optimiser lands near their answers from nothing but the protein sequence and a codon table, then the design space is mostly a matter of arithmetic, which is good news for a project that has to make many designs.

## 1. The benchmark

The SARS-CoV-2 spike protein, in its 2P form, as encoded by Comirnaty. The input to the optimiser is the protein only. The comparisons are the published designs.

| Sequence | GC | CAI | Forbidden motifs |
| --- | --- | --- | --- |
| Wild-type viral spike gene | 0.3734 | 0.6734 | 0 |
| Naive: best codon everywhere | 0.6169 | 1.0000 | 0 |
| This optimiser | **0.6140** | **0.9981** | 0 |
| Comirnaty, as published | 0.5700 | 0.9515 | 0 |

**Codon agreement with Comirnaty's published design: 78.2 per cent** of 1273 positions, derived from the protein alone.

Three things are worth reading off that table.

The virus codes this protein for human cells at a CAI of 0.6734 and GC 0.3734. That is why codon optimisation exists: the gene is not bad, it is written for a different host. The gap between the viral coding sequence and the optimised one is the entire reason this project is worth doing.
Choosing the most used codon everywhere gives a CAI of 1.0000, which sounds perfect and is not. It ignores GC content, which is a constraint the design has to live with, and it ignores motifs, which are a manufacturing failure. The optimiser reaches 0.9981 while holding every GC block inside the band and leaving no forbidden motif behind. A lower CAI bought with a controllable GC profile is a better construct, not a worse one.
Comirnaty sits at CAI 0.9515 and GC 0.5700. Our design lands at a higher CAI and a slightly higher GC, and agrees with theirs on roughly four codon choices in five. The places where we differ are the places where they deliberately took a less favoured codon, which is exactly where the interesting engineering is.

## 2. What the optimiser guarantees

| Property | This design |
| --- | --- |
| Length | 3819 nt, encoding 1273 aa |
| GC per block | 0.500 to 0.700 (band 0.40 to 0.70) |
| GC per 50-nt window | 0.440 to 0.800 |
| CAI | 0.9981 (target 0.95) |
| Forbidden motifs | 0 |
| Homopolymer runs over 8 | 0 |
| Poly-lysine runs over 4 codons | 0 |
| CAI optimal subject to the GC band | yes |
| Design is clean | yes |

The GC band is not a literature number. It is the band that both approved products satisfy: Comirnaty's 30-nucleotide blocks span 0.400 to 0.700 and Spikevax's span 0.467 to 0.767. A design that stays inside it is never more extreme than something a regulator has already approved. The narrower 0.50 to 0.60 band that appears in the literature would have excluded Spikevax, whose coding sequence sits at 0.6230.

One honest limit: when the band is tightened to within one nucleotide of a target, some blocks become infeasible, because a stretch of residues with only one or two possible codons cannot reach the target GC. Tightened to within one nucleotide of the target GC, this design reports 10 infeasible blocks and reports its CAI as not optimal, rather than hiding either. A tool that quietly returns a sequence that violates its own policy is worse than no tool.

## 3. The forbidden motifs

Enforced:

- **polyadenylation signal** (`AATAAA`). AAUAAA is the canonical polyadenylation signal of vertebrate pre-mRNAs, sitting 10 to 30 nucleotides upstream of the cleavage site. A copy inside a coding region can be used as a cryptic polyadenylation site, which yields a truncated transcript with no in-frame stop codon, and a ribosome that then reads into the poly(A) tail and stalls on a poly-lysine tract.
  Source: Proudfoot NJ, Brownlee BG. 1976. AATAAA is a signal for synthesis of polyadenylated RNA. Nucleic Acids Res 3(7):1183-1190. Fitzgerald PC, Shenk T. 1981. The sequence AAUAAA is required for accurate and efficient polyadenylation of the 5' long control region of SV40. Proc Natl Acad Sci USA 78(7):4182-4185. For the observation that polyadenylation inside a coding region produces a transcript without an in-frame stop codon, see 'Signals for pre-mRNA cleavage and polyadenylation', https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4451228/

Documented but not enforced, because the constraint depends on information this project does not have yet:

- **NotI site** (`GCGGCCGC`). A cloning site used in many linear DNA templates. A copy inside the coding region has to be removed by changing a codon, not the protein. Template dependent. The restriction sites actually present in the chosen IVT template must be confirmed with the manufacturing partner before this motif can be enforced.
- **ApaI site** (`GGGCCC`). A cloning site used in some linear DNA templates. Template dependent, see the NotI entry above.
- **SapI site** (`GCTCTTC`). A type IIS cloning site used in some linear DNA templates. It is also the target of the guide-directed nuclease used in some template workflows. Template dependent, see the NotI entry above.

Neither approved product contains a single AATAAA, a homopolymer run of nine or more, or even two consecutive AAA codons. That is not a coincidence to copy blindly, it is a set of defaults that two regulators accepted.

## 4. The teaching artifact: GFP

Green fluorescent protein from Aequorea victoria, UniProt P42212, 238 residues. The classic example, and small enough to read in full.

| Sequence | GC | CAI | Forbidden motifs |
| --- | --- | --- | --- |
| Naive: best codon everywhere | 0.6092 | 1.0000 | 0 |
| This optimiser | 0.6050 | 0.9964 | 0 |

## 5. What this does not claim

No expression level is predicted. No structure is predicted. Nothing here has been synthesised, and nothing here has been administered to a cell. The optimiser makes a sequence well formed with respect to a set of cited, measurable constraints, and that is the whole claim. The ladder from here is L1, expression in cells, which belongs to a wet-lab partner.

## References

- Sharp PM, Li WH. 1987. The codon adaptation index, a measure of directional synonymous codon usage bias. Mol Evol Biol 26(3):151-161.
- Kudla G, Wozniak M, Hossmann T, Krainer AR. 1999. High guanine and cytosine content increases mRNA levels in mammalian cells. Nucleic Acids Res 27(9):2085-2088.
- Proudfoot NJ, Brownlee BG. 1976. AATAAA is a signal for synthesis of polyadenylated RNA. Nucleic Acids Res 3(7):1183-1190.
- Fitzgerald PC, Shenk T. 1981. The sequence AAUAAA is required for accurate and efficient polyadenylation of the 5' long control region of SV40. Proc Natl Acad Sci USA 78(7):4182-4185.
- Mechanism of ribosome stalling during translation of a poly(A) tail, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6900289/
- Engineering the poly(A) tail for therapeutic mRNA, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC13212533/

_In-silico designs. No therapeutic efficacy demonstrated. Not medical advice._
