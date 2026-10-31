# Anatomy of Comirnaty, reverse-engineered in Python

> In-silico designs. No therapeutic efficacy demonstrated. Not medical advice.

This is the M0 artifact of the antivenom-mrna engine. It takes the primary regulatory document for a licensed mRNA vaccine, reconstructs the construct from it, checks every coordinate the document publishes, cross-checks the result against an independent sequencing assembly and against the authoritative reference protein, and measures the coding region. Every number below was computed by `scripts/vax_anatomy.py` from the committed reference data in `data/reference/`.

## Sources

| Data | Source | Nature |
| --- | --- | --- |
| Comirnaty full mRNA, 4,284 nt | WHO International Nonproprietary Names Programme document 11889 (Comirnaty), archived at https://berthub.eu/articles/11889.doc | primary regulatory document |
| Comirnaty coding sequence | transcription of the same document in [berthubert/bnt162b2](https://github.com/berthubert/bnt162b2) | human transcription |
| Read-derived mRNA assemblies | [NAalytics, Assemblies of putative SARS-CoV-2 spike-encoding mRNA sequences](https://github.com/NAalytics/Assemblies-of-putative-SARS-CoV2-spike-encoding-mRNA-sequences-for-vaccines-BNT-162b2-and-mRNA-1273) | independent sequencing assembly |
| SARS-CoV-2 spike reference | UniProt [P0DTC2](https://rest.uniprot.org/uniprotkb/P0DTC2.fasta) | authoritative reference protein |
| Codon usage | Kazusa Codon Usage Database, *Homo sapiens* (Nakamura et al. 2000), 93,487 coding sequences | reference table |
| Antigen mutations | TGA Nonclinical Evaluation Report for BNT162b2; MHRA Public Assessment Report for COVID-19 mRNA Vaccine BNT162b2 | regulatory assessment |

## 1. The published architecture, reproduced

The document gives a name and a coordinate range for every part of the construct. Slicing the sequence at those coordinates turns a drawing into a test.

| Element | Published range | Published length | Engine measurement | Match |
| --- | --- | --- | --- | --- |
| cap | 1 to 2 | 2 nt | 2 nt | yes |
| 5'-UTR | 3 to 54 | 52 nt | 52 nt | yes |
| signal peptide | 55 to 102 | 48 nt | 48 nt | yes |
| S protein | 103 to 3879 | 3777 nt | 3777 nt | yes |
| 3'-UTR | 3880 to 4174 | 295 nt | 295 nt | yes |
| poly(A) | 4175 to 4284 | 110 nt | 110 nt | yes |

Total length: 4284 nt. Every published coordinate is reproduced. The poly(A) tail is 30 adenosines, the 10-nucleotide linker `GCAUAUGACU`, then 70 more adenosines, and the two stop codons sit at 3874 to 3879.

## 2. Three independent cross-checks

| Check | Result |
| --- | --- |
| Primary document against independent sequencing assembly | 4175 nt compared, 0 mismatches |
| Document's viral S coding sequence against UniProt P0DTC2 | identical protein, 1,273 aa |
| Document's viral S against the vaccine antigen | 2 amino-acid changes: K986P, V987P |
| Engine's ORF start against the published signal peptide start | nt 55 against nt 55 |

Three sources that were produced by different people for different reasons agree with each other. That is the standard this engine will be held to for every antivenom construct it produces: one source is not enough, and the sequence must survive being checked against a second one.

## 3. Measurements

### BNT162b2 (Comirnaty)

| Feature | BNT162b2 |
| --- | --- |
| Coding sequence | 3822 nt |
| Encoded protein | 1273 aa |
| 5' UTR in this source | 54 nt |
| 3' UTR in this source | 299 nt |
| poly(A) in this source | 0 nt |
| GC content, whole CDS | 0.5699 |
| GC content, first 50 nt | 0.5600 |
| GC content, sliding-window min / max | 0.3800 / 0.6600 |
| CAI (Homo sapiens, Kazusa) | 0.9515 |
| Kozak context (-3 to +4) | `ACCATGTTC` |
| Kozak consensus | deviates: +4 is not G |
| Differences vs wild-type spike | 2 |
| Which | K986P, V987P |

### mRNA-1273 (Spikevax)

Read-derived assembly. The authors of the assembly state that the Moderna sequence could not be checked against a published reference, so this row is reported as an observation, not as a specification.

| Feature | mRNA-1273 |
| --- | --- |
| Coding sequence | 3822 nt |
| Encoded protein | 1273 aa |
| 5' UTR in this source | 57 nt |
| 3' UTR in this source | 125 nt |
| poly(A) in this source | 0 nt |
| GC content, whole CDS | 0.6230 |
| GC content, first 50 nt | 0.6000 |
| GC content, sliding-window min / max | 0.4400 / 0.7400 |
| CAI (Homo sapiens, Kazusa) | 0.9786 |
| Kozak context (-3 to +4) | `ACCATGTTC` |
| Kozak consensus | deviates: +4 is not G |
| Differences vs wild-type spike | 2 |
| Which | K986P, V987P |

## 4. Calibration for the antivenom constructs

The point of measuring licensed products is to set the engine's design targets against reality instead of against a textbook. Two licensed mRNA vaccines, both in clinical use, agree on the following:

| Parameter | Engine target | BNT162b2 | mRNA-1273 |
| --- | --- | --- | --- |
| GC content of CDS | 0.50 to 0.60 | 0.5699 | 0.6230 |
| CAI against human codon usage | 0.70 | 0.9515 | 0.9786 |
| 5' UTR | 50 to 150 nt | 54 nt | 57 nt |
| Initiation context | Kozak gccRccATGG | `ACCATGTTC` | `ACCATGTTC` |

Two observations worth carrying forward. First, both vaccines use a purine at position -3 of the start codon and both deviate at position +4, so the consensus is a guide rather than a rule. Second, the two products encode the *same* protein with different codon choices, 363 synonymous positions apart, and the GC-richer of the two is the one with the higher CAI. Codon choice is a free variable in this design space, and it is the variable this engine controls.

## 5. What the engine can now do

- Read a regulatory-grade primary document and reconstruct the construct it describes.
- Locate a coding region and prove that it translates to the intended protein, with the differences from the wild-type protein stated as positions and amino acids.
- Verify a sequence against a second, independently produced source and report the number of mismatches.
- Measure GC content overall and locally, CAI against a host codon usage table, and the initiation context of the start codon.

And what it cannot do yet: it cannot yet design anything. That is milestone M1 and M2, and the targets above are what those milestones will be measured against.

## References

- WHO International Nonproprietary Names Programme, document 11889, Comirnaty (BNT162b2), February 2021.
- Australian TGA, Nonclinical Evaluation Report for BNT162b2 (K986P and V987P, expressed protein P2 S).
- MHRA, COVID-19 mRNA Vaccine BNT162b2 Public Assessment Report (pre-fusion stabilised variant, K986P and V987P).
- Nakamura Y, Gojobori T, Ikemura T. 2000. Codon usage tabulated from international DNA sequence databases: status for the year 2000. Nucleic Acids Res 28(1):292.
- Sharp PM, Li WH. 1987. The codon adaptation index. Mol Evol Biol 26(3):151-161.
- Kozak M. 1986. Point mutation defines a sequence flanking the AUG initiator codon. Cell 44(2):283-292.
- Kudla G, Wozniak M, Hossmann T, Krainer AR. 1999. High guanine and cytosine content increases mRNA levels in mammalian cells. Nucleic Acids Res 27(9):2085-2088.

_In-silico designs. No therapeutic efficacy demonstrated. Not medical advice._
