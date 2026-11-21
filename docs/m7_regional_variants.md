# M7: regional venom variation as a design input

> In-silico designs. No therapeutic efficacy demonstrated. Not medical advice.

## 1. Why the species-level matrix is not enough

The coverage matrix of M4 and the planner of M5 ask their questions per species. A species is not one venom. The spectacled cobra of Punjab carries about 80 per cent neurotoxic three-finger toxins by venom mass; the desert population of Rajasthan spreads its venom across cytotoxins, about 42 per cent, and phospholipase A2, about 20 per cent, and the commercial antivenom fails against it completely at five times the lethal dose. The Russell's viper of the Gangetic Plain carries about 39 per cent Kunitz-type inhibitors, the Western Ghats population 59 per cent phospholipase A2. The same species, different venoms, and a cocktail sized for the average is sized for neither.

This report does one thing with the printed regional data: it weights the coverage matrix by venom mass, per population. The species-level question is which cells a plan covers. The regional question is which cells a plan covers *in the venom a patient in a given state is actually carrying*.

## 2. What the surveys measured, and what they did not

Both surveys are open access and both print their mass spectrometry for three populations per species, chosen for their distinct chromatographic and toxicity profiles. The other populations of both surveys were characterised by chromatography alone, and no composition is recorded for them here, because a fraction nobody printed is not a fraction anybody can check.

The surveys' toxin categories are not the matrix's toxin families, and the correspondence this report uses is an approximation it names. The cobra survey's neurotoxic three-finger toxins include both long-chain and short-chain alpha-neurotoxins; the matrix separates 3FTx-L, the long-chain neurotoxins, from 3FTx-S, the cytotoxins. This report maps the survey's neurotoxic fraction onto 3FTx-L and its cytotoxic fraction onto 3FTx-S. A survey category with no matrix family at all, such as the snake venom serine proteases and C-type lectins of the viper survey, is recorded as unmapped rather than squeezed into the nearest family. The unmapped fraction is itself a finding: it is the part of the venom the matrix does not model.

## 3. The populations, weighted tier by tier

Every mapped family below is checked against the plan at all four tiers, and the table records the strictest tier at which the family is still covered. The tiers are the planner's: published, sequence-available, single-polypeptide, licence-unblocked. The single-polypeptide tier is the one where a transcript can carry the payload, and its plan is TPL1158_01_C09, E113. A family the planner excludes on purpose, the anticoagulant Kunitz-type fraction, is recorded as excluded, not as uncovered: the planner refuses it because an antibody against it would make the coagulopathy worse, which is the M5 rule and it applies region by region.

### PB (Punjab, North, semi-arid) - Naja naja

| Survey category | Venom mass | Matrix family | Status |
| --- | --- | --- | --- |
| neurotoxic 3FTx | 0.80 | 3FTx-L | covered through `sequence-available` |
| cytotoxic 3FTx | 0.10 | 3FTx-S | covered through `sequence-available` |

Note: The most neurotoxic cobra venom of the survey. The PLA2 fraction is inside the survey's 0.04 to 20 per cent range but the per-population value is printed only in the supplementary tables, so it is not recorded here.

### WB (West Bengal, East, Gangetic Plain) - Naja naja

| Survey category | Venom mass | Matrix family | Status |
| --- | --- | --- | --- |
| neurotoxic 3FTx | 0.73 | 3FTx-L | covered through `sequence-available` |
| cytotoxic 3FTx | 0.24 | 3FTx-S | covered through `sequence-available` |

### RJ (Rajasthan, Desert) - Naja naja

| Survey category | Venom mass | Matrix family | Status |
| --- | --- | --- | --- |
| cytotoxic 3FTx | 0.42 | 3FTx-S | covered through `sequence-available` |
| PLA2 | 0.20 | PLA2 | covered through `published` |

Note: The venom the commercial antivenom failed against completely at five times the lethal dose. Neurotoxic 3FTx was reported only as about 30 per cent, which is not a printed fraction and is therefore not recorded; the mapped total of this profile is accordingly an undercount, not an overcount.

### WB (West Bengal, East, Gangetic Plain) - Daboia russelii

| Survey category | Venom mass | Matrix family | Status |
| --- | --- | --- | --- |
| Kunitz-type inhibitors | 0.39 | KSPI | excluded by design |
| PLA2 | 0.35 | PLA2 | covered through `published` |
| SVSP | 0.08 | none | unmapped: the matrix does not model this family |
| snaclec | 0.05 | none | unmapped: the matrix does not model this family |

Note: The Kunitz-type fraction is the anticoagulant fraction the planner excludes on purpose. The serine protease and C-type lectin fractions have no matrix family and stay unmapped.

### MH (Maharashtra, Southwest, Western Ghats) - Daboia russelii

| Survey category | Venom mass | Matrix family | Status |
| --- | --- | --- | --- |
| PLA2 | 0.59 | PLA2 | covered through `published` |
| SVSP | 0.27 | none | unmapped: the matrix does not model this family |
| Kunitz-type inhibitors | 0.04 | KSPI | excluded by design |
| snaclec | 0.02 | none | unmapped: the matrix does not model this family |

### MP (Madhya Pradesh, Central, Deccan Plateau) - Daboia russelii

| Survey category | Venom mass | Matrix family | Status |
| --- | --- | --- | --- |
| PLA2 | 0.45 | PLA2 | covered through `published` |
| snaclec | 0.20 | none | unmapped: the matrix does not model this family |
| SVSP | 0.13 | none | unmapped: the matrix does not model this family |
| Kunitz-type inhibitors | 0.07 | KSPI | excluded by design |
| VEGF | 0.08 | none | unmapped: the matrix does not model this family |

Note: The snaclec and VEGF fractions have no matrix family. The printed snaclec value is 'nearly about 20 per cent', recorded here as 0.20.

## 4. What the weighting changes

At the species level, the spectacled cobra's three cells all have candidates at the published tier, so the species reads as covered. Weighted by population and read down the tiers, the picture is different, and the difference is the point:

- **Every cobra cell is covered only by candidates that cannot be transcribed as they stand, so the whole species goes uncovered at the tier where a transcript can carry it.** The two neurotoxin cells are held by Fabs, which have sequences but two-chain formats; the phospholipase A2 cell is held by an scFv whose sequence was never located. In Punjab, about 80 per cent of the venom mass sits behind the 3FTx-L cell when it drops. In Rajasthan it is about 42 per cent behind the cytotoxin cell and 20 per cent behind phospholipase A2, and the antivenom fails against this venom completely at five times the lethal dose. The gap is not a row in a matrix; it is the venom a patient carries.
- **For the viper, the biggest mapped fraction in two of three populations is the family the planner excludes on purpose.** The Gangetic Plain venom is about 39 per cent Kunitz-type inhibitors, and no antibody against them belongs in a design, because they are the anticoagulant fraction of a venom whose patient is already bleeding. The biggest actionable mapped fraction in every viper population is phospholipase A2, covered only at the published tier, because its one recombinant candidate has no publicly locatable sequence: the same data-availability problem M4 and M5 named, now with a mass attached, 59 per cent of the Western Ghats venom, 45 per cent of the Deccan Plateau's, 35 per cent of the Gangetic Plain's.
- **The viper matrix under-models the venom.** Snake venom serine proteases are 8 to 27 per cent of the mapped populations' venom by mass and C-type lectins up to 20 per cent, and neither family has a column in the gap matrix. Adding them is a decision for the next milestone, and this report records the need rather than silently widening the matrix underneath M4's published artifacts.

## 5. What this report is not

- **It is not a regional cocktail design.** Six populations, three per species, with printed composition, is the survey's record, not India's map. The bioclimatic modelling of Sarangi and colleagues in 2025 extends the coverage; the composition data behind this report is deliberately limited to what was printed.
- **It is not an antivenom evaluation.** The regional failures the surveys report, the below-claim potency in Punjab, the complete failure against the Rajasthan cobra at five times the lethal dose, the unneutralised clotting against the Gangetic Plain and coastal vipers, are the surveys' findings and this report's motivation. An mRNA construct is not the answer to an acute bite regardless of region, and the README's caveat applies here with full force.

## References

- Senji Laxme RR, Attarde S, Khochare S, et al. 2021. Biogeographical venom variation in the Indian spectacled cobra (Naja naja) underscores the pressing need for pan-India efficacious snakebite therapy. PLoS Negl Trop Dis 15(2):e0009150. doi:10.1371/journal.pntd.0009150.
- Senji Laxme RR, Khochare S, Attarde S, et al. 2021. Biogeographic venom variation in Russell's viper (Daboia russelii) and the preclinical inefficacy of antivenom therapy in snakebite hotspots. PLoS Negl Trop Dis 15(3):e0009247. doi:10.1371/journal.pntd.0009247.
- Sarangi MK, et al. 2025. Predictive modelling of bioclimatic venom variation in Russell's viper. PLoS Negl Trop Dis. doi:10.1371/journal.pntd.0012949.
