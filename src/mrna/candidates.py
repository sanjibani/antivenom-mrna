"""The curated candidate table for India's Big Four venomous snakes.

This module is the source of truth for the M4 literature map. The tab-separated
table under ``data/`` and the markdown report under ``docs/`` are both generated
from it by ``scripts/candidate_tracker.py``, so the published table cannot drift
away from the code that validates it.

The table records **published antibody candidates** that have a claimed
neutralising activity against a toxin, or against the whole venom, of a Big
Four species. It is not a catalogue of everything ever made against snake
venom. The scope rule is deliberate: a row exists when the source makes a
neutralisation or protection claim against a Big Four species or a toxin from
one. Small molecules, aptamers, plant extracts and antivenom pools are out of
scope, because this project designs the mRNA layer for an antibody and a
small molecule is not an antibody.

Every row carries the source it came from. Where a sequence could not be
located in a public database, the row says ``none`` rather than guessing,
because "the sequence is not publicly available" is a finding, and it is the
finding that decides whether a candidate can be transcribed at all.

Note that the table is a snapshot of a survey performed on 2026-11-01. It is
not exhaustive, and a row saying ``none`` means "not located during this
survey", not "does not exist".
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = [
    "BIG_FOUR",
    "CANDIDATES",
    "CITATIONS",
    "CLINICAL_SET",
    "COMMON_NAMES",
    "EVIDENCE_LEVELS",
    "FORMATS",
    "MRNA_SUITABILITY",
    "PATENT_STATUS",
    "SPECIES_FAMILIES",
    "SPECIES_TOXIN_FAMILIES",
    "TOXIN_FAMILIES",
    "TOXIN_FAMILY_NOTES",
    "Candidate",
]


@dataclass(frozen=True)
class Candidate:
    """One published antibody candidate against a Big Four toxin or venom.

    Attributes:
        species: Binomial names of the species the candidate targets. A
            semicolon separates several species, which is how a candidate that
            cross-neutralises two of them is recorded.
        toxin_family: The toxin family the candidate targets.
        toxin_target: The specific toxin or venom component.
        antibody_id: The clone or candidate name as the source publishes it.
        molecule_format: The antibody format, one of :data:`FORMATS`.
        sequence_source: Where the variable-region sequence can be found.
            Each entry is a typed token, ``kind:value``, with kind one of
            ``PDB``, ``patent``, ``paper`` or ``none``.
        patent_status: What the public patent record shows, one of
            :data:`PATENT_STATUS`. This is a fact about filings, never a
            freedom-to-operate opinion.
        tested_against: The species, toxin and assay the source reports.
        evidence_level: The highest rung of the validation ladder the
            candidate has reached, one of :data:`EVIDENCE_LEVELS`.
        mrna_suitability: Whether the candidate can be transcribed from
            public material, one of :data:`MRNA_SUITABILITY`. The engine
            recomputes this from the format and the sequence source and
            refuses to load a table where the two disagree.
        citation: Key into :data:`CITATIONS`.
        notes: Anything a reader needs in order to not be misled.
    """

    species: tuple[str, ...]
    toxin_family: str
    toxin_target: str
    antibody_id: str
    molecule_format: str
    sequence_source: tuple[str, ...]
    patent_status: str
    tested_against: str
    evidence_level: str
    mrna_suitability: str
    citation: str
    notes: str = ""
    tags: tuple[str, ...] = field(default=())


#: The four snakes that account for almost all snakebite death in India.
#:
#: Sources: Suraj et al. 2019, *Lancet* 393:2273-2280 (burden of snakebite
#: mortality in India); and the National Commission for the Control of
#: Snakebite Envenoming reporting, which names the same four.
BIG_FOUR: tuple[str, ...] = (
    "Naja naja",
    "Bungarus caeruleus",
    "Daboia russelii",
    "Echis carinatus",
)

#: The Big Four plus the monocled cobra.
#:
#: The monocled cobra, *Naja kaouthia*, is not one of the Big Four, but it is
#: the other major Indian cobra, it is one of the species Indian polyvalent
#: antivenom is raised against, and it accounts for a large share of cobra
#: bites in the country. A cocktail aimed at Indian snakebite that ignored it
#: would be incomplete in a way that matters, so the planner can be asked for
#: the Big Four alone or for this wider set.
#:
#: Sources: Suraj et al. 2019, *Lancet* 393:2273-2280; and the Indian
#: polyvalent antivenom literature, which raises pools against both Indian
#: cobras.
CLINICAL_SET: tuple[str, ...] = (*BIG_FOUR, "Naja kaouthia")

#: English names, used only in the generated report.
COMMON_NAMES: dict[str, str] = {
    "Naja naja": "spectacled cobra",
    "Naja kaouthia": "monocled cobra",
    "Bungarus caeruleus": "common krait",
    "Bungarus multicinctus": "many-banded krait",
    "Daboia russelii": "Russell's viper",
    "Daboia russelii formosensis": "Russell's viper, Taiwan subspecies",
    "Echis carinatus": "saw-scaled viper",
    "Echis carinatus sochureki": "saw-scaled viper, sochurek",
}

#: The venom family each species belongs to, used only in the generated report.
#:
#: This is a fact about the animals rather than a judgement about them. The
#: three elapid genera carry the neurotoxic and cytotoxic toxins; the two
#: viper genera carry the haemostatic and haemorrhagic ones. The report
#: groups its findings by this, because a reader who sees "six cells
#: uncovered" learns less than a reader who sees which venom and which
#: toxin class went untouched.
SPECIES_FAMILIES: dict[str, str] = {
    "Naja naja": "elapid",
    "Naja kaouthia": "elapid",
    "Bungarus caeruleus": "elapid",
    "Bungarus multicinctus": "elapid",
    "Daboia russelii": "viper",
    "Daboia russelii formosensis": "viper",
    "Echis carinatus": "viper",
    "Echis carinatus sochureki": "viper",
}

#: Toxin families that the tracker knows about.
TOXIN_FAMILIES: tuple[str, ...] = (
    "3FTx-L",
    "3FTx-S",
    "PLA2",
    "SVMP",
    "KSPI",
    "LAAO",
    "mixed",
    "unspecified",
)

#: Antibody formats, grouped by how many open reading frames they need.
FORMATS: tuple[str, ...] = (
    "scFv",
    "VHH",
    "VHH-Fc",
    "Fab",
    "mAb-IgG",
    "polyclonal-IgG",
    "polyclonal-Fab",
    "IgY",
)

#: Formats that are one polypeptide and therefore one open reading frame.
SINGLE_ORF_FORMATS: frozenset[str] = frozenset({"scFv", "VHH", "VHH-Fc"})

#: Formats with no defined sequence, so there is nothing to transcribe.
POLYCLONAL_FORMATS: frozenset[str] = frozenset(
    {"polyclonal-IgG", "polyclonal-Fab", "IgY"},
)

#: The validation ladder, reused from the repository README so that a rung
#: means the same thing here as it does for the designs in this repository.
EVIDENCE_LEVELS: dict[str, str] = {
    "L0": "in silico only; no experimental data",
    "L1": "expressed and purified",
    "L2": "binding measured, in vitro",
    "L3": "neutralisation measured, ex vivo or in cells",
    "L4": "protection measured in an animal model",
    "L5": "manufacturing, formulation or chemistry demonstrated",
}

#: Whether a candidate can be transcribed from public material.
MRNA_SUITABILITY: dict[str, str] = {
    "ready": "a defined sequence is public and the format is a single open reading frame",
    "adaptable": (
        "a defined sequence is public, but the format needs more than one ORF "
        "or an engineering step"
    ),
    "blocked": "no defined sequence exists to transcribe, because the candidate is polyclonal",
    "unknown": "the sequence was not located, so nothing can be said about buildability",
}

#: What the public patent record shows. This is never a clearance.
PATENT_STATUS: tuple[str, ...] = (
    "none-located",
    "application-pending",
    "granted",
    "not-assessed",
)

#: Which toxin families matter for which species.
#:
#: The four elapid families are the three-finger toxins and phospholipase A2.
#: The two viper families are dominated by phospholipase A2, snake venom
#: metalloproteinases and Kunitz-type serine protease inhibitors, with
#: L-amino acid oxidase contributing to the coagulopathy of *Daboia* bites.
#:
#: Sources: Lomonte and Rangel 2012, *Toxicon* 60:520-530 (PLA2 myotoxins);
#: Patra, Chanda and Mukherjee 2019, *Expert Rev Proteomics* 16:457-469 (krait
#: venom composition); Chanda et al. 2024 and the proteomics literature for
#: *Daboia russelii*; and Kumar et al. 2024, *Mol Immunol* 165:55-67, which
#: reports phospholipase A2 as the most abundant toxin in both *Naja naja* and
#: *Echis carinatus* venom.
SPECIES_TOXIN_FAMILIES: dict[str, tuple[str, ...]] = {
    "Naja naja": ("3FTx-L", "3FTx-S", "PLA2"),
    "Naja kaouthia": ("3FTx-L", "3FTx-S", "PLA2"),
    "Bungarus caeruleus": ("3FTx-L", "PLA2", "KSPI"),
    "Bungarus multicinctus": ("3FTx-L", "PLA2", "KSPI"),
    "Daboia russelii": ("PLA2", "SVMP", "KSPI", "LAAO"),
    "Daboia russelii formosensis": ("PLA2", "SVMP", "KSPI", "LAAO"),
    "Echis carinatus": ("PLA2", "SVMP", "KSPI", "LAAO"),
    "Echis carinatus sochureki": ("PLA2", "SVMP", "KSPI", "LAAO"),
}

#: What each toxin family is, and what an antibody against it would or would
#: not achieve.
#:
#: This matters for reading the gap matrix honestly. An empty cell is not
#: always an opportunity. The Kunitz-type serine protease inhibitors are the
#: anticoagulant fraction of viper venom, because they inhibit the clotting
#: factors directly, so neutralising them would make the coagulopathy worse
#: rather than better, and an empty cell there may be a correct decision by
#: the field rather than a hole in it.
#:
#: Sources: Lomonte and Rangel 2012, *Toxicon* 60:520-530 (PLA2 myotoxins);
#: Gutiérrez and Lomonte 2013, *Toxicon* 42:915-931 (PLA2 myotoxicity);
#: Escalante, Rucavado and Fox 2011, *J Proteomics* 74:1781-1796 (SVMP
#: microvascular damage); and the standard venom nomenclature reviews for
#: the Kunitz-type and L-amino-acid-oxidase fractions.
TOXIN_FAMILY_NOTES: dict[str, str] = {
    "3FTx-L": (
        "Long-chain three-finger toxins are the neurotoxins of the elapids. "
        "They block the nicotinic acetylcholine receptor and cause flaccid "
        "paralysis."
    ),
    "3FTx-S": (
        "Short-chain three-finger toxins are the cytotoxins of the elapids. "
        "They damage cell membranes and cause the local tissue necrosis of a "
        "cobra bite."
    ),
    "PLA2": (
        "Phospholipase A2 is the dominant toxin of Indian cobra and saw-scaled "
        "viper venom, and the dominant fraction of common krait venom by "
        "mass. It is myotoxic, hemolytic and anticoagulant at the same time."
    ),
    "SVMP": (
        "Snake venom metalloproteinases are the haemorrhagic toxins of the "
        "viper bites. They degrade the basement membrane of small blood "
        "vessels, which is what produces the bleeding."
    ),
    "KSPI": (
        "Kunitz-type serine protease inhibitors are the anticoagulant "
        "fraction of viper venom, because they inhibit the clotting factors "
        "directly. Neutralising them would make the coagulopathy worse "
        "rather than better, so an empty cell here may be a correct decision "
        "by the field rather than a hole in it."
    ),
    "LAAO": (
        "L-amino acid oxidase contributes to the coagulopathy of viper bites "
        "by damaging the vessel wall rather than by acting on the clotting "
        "cascade itself."
    ),
    "mixed": "A candidate or a cocktail that covers more than one toxin family.",
    "unspecified": (
        "The source does not resolve the target toxin, so the row cannot be "
        "placed in the gap matrix."
    ),
}

#: Primary sources, keyed by the short citation used in each row.
CITATIONS: dict[str, str] = {
    "Kumar2024": (
        "Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2024. Recombinant human "
        "scFv antibody fragments against phospholipase A2 from Naja naja and "
        "Echis carinatus snake venoms: In vivo neutralization and mechanistic "
        "insights. Mol Immunol 165:55-67. "
        "doi:10.1016/j.molimm.2023.12.006. PMID 38154407."
    ),
    "Kumar2023preprint": (
        "Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2023. Antibodies against "
        "Elapidae and Viperidae snake venoms: in vivo neutralization and "
        "mechanistic insights. bioRxiv 2023.06.16.545250. "
        "doi:10.1101/2023.06.16.545250."
    ),
    "Glanville2025": (
        "Glanville J, Kwong K, et al. 2025. Broad-spectrum neutralisation of "
        "elapid venoms by a minimal antibody combination plus a small-molecule "
        "inhibitor. Cell. Centivax anti-venom antibodies, "
        "US20240294618A1 and EP4319816A4."
    ),
    "Ledsgaard2022": (
        "Ledsgaard L, Laustsen AH, Pus U, et al. 2022. In vitro discovery of a "
        "human monoclonal antibody that neutralizes lethality of cobra snake "
        "venom. MAbs 14:2085536. doi:10.1080/19420862.2022.2085536."
    ),
    "Ledsgaard2023": (
        "Ledsgaard L, Wade J, Jenkins TP, et al. 2023. Discovery and "
        "optimization of a broadly-neutralizing human monoclonal antibody "
        "against long-chain alpha-neurotoxins from snakes. Nat Commun 14:682. "
        "doi:10.1038/s41467-023-36393-4."
    ),
    "Khalek2024": (
        "Khalek A, et al. Synthetic development of a broadly neutralizing "
        "antibody against snake venom long-chain alpha-neurotoxins. "
        "Sci Transl Med. PDB 8SXP."
    ),
    "Laustsen2024": (
        "Laustsen AH, et al. 2024. Engineering of a nanobody-based "
        "recombinant antivenom. Nat Commun. "
        "doi:10.1038/s41467-024-48539-z."
    ),
    "Samanta2026": (
        "Samanta A, Rudresha GV, Ahmadi S, et al. 2026. Oligoclonal "
        "nanobody-based recombinant antivenom protects mice challenged with "
        "venom from cobras and king cobras from India. Sci Transl Med "
        "18:866. doi:10.1126/scitranslmed.aed4290. PMID 42715346."
    ),
    "Richard2013": (
        "Richard G, Meyers AJ, McLean MD, et al. 2013. In vivo neutralization "
        "of alpha-cobratoxin with high-affinity llama single-domain antibodies "
        "(VHHs) and a VHH-Fc antibody. PLoS ONE 8:e69495. "
        "doi:10.1371/journal.pone.0069495."
    ),
    "Lee2017": (
        "Lee CH, Lee YC, Lee YL, et al. 2017. Single chain antibody fragment "
        "against venom from the snake Daboia russelii formosensis. Toxins "
        "9:347. doi:10.3390/toxins9110347."
    ),
    "Li2019": (
        "Li W, et al. 2019. Secretory PLA2 specific single domain antibody "
        "neutralizes Russell viper venom induced cellular and organismal "
        "toxicity. Toxicon. PMID 31689425."
    ),
    "RCSB8IA6": (
        "RCSB Protein Data Bank entry 8IA6. Crystal structure of scFv "
        "antibody against Phospholipase A2 of Echis carinatus venom. Released "
        "2024-02-14."
    ),
    "RCSB9GCN": (
        "RCSB Protein Data Bank entry 9GCN. Co-crystal VHH "
        "(TPL1158_01_C09) with alpha-cobratoxin from Naja kaouthia. Released "
        "2025-08-06."
    ),
    "RCSB9CLP": (
        "RCSB Protein Data Bank entry 9CLP. Structure of ecarin from the "
        "venom of Kenyan saw-scaled viper in complex with the Fab of "
        "neutralizing antibody H11. Released 2024-09-11."
    ),
    "RCSB9HUO": (
        "RCSB Protein Data Bank entries 9HUO, 9FYT and 9HXO. A01 monoclonal "
        "antibody bound to alpha-cobratoxin at pH 5.5, 4.5 and 6.0."
    ),
    "DTU2022IgG": (
        "WO2023209212A1. Neutralisation of alpha-neurotoxins using human "
        "recombinant IgG antibodies. Danmarks Tekniske Universitet. Priority "
        "date 2022-04-29, published 2023-11-02. The application names "
        "2554_01_D11, 368_01_C05 and A01."
    ),
    "DTU2025Nanobody": (
        "WO2025078529A1. Nanobodies capable of cross-binding to several "
        "different elapid three-finger toxins. Danmarks Tekniske Universitet. "
        "Priority date 2023-10-10, published 2025-04-17. The application "
        "names TPL1158_01_C09."
    ),
    "Ledsgaard2024": (
        "Ledsgaard L, Wade J, Jenkins TP, et al. 2024. Phage display assisted "
        "discovery of a pH-dependent anti-alpha-cobratoxin antibody."
    ),
    "IMPACT1864": (
        "Investigational candidates for snakebite envenoming, IMPACT Global "
        "Health, PCR ID 1864. Human recombinant polyclonal Fab against "
        "Echis carinatus."
    ),
    "IMPACT1708": (
        "Investigational candidates for snakebite envenoming, IMPACT Global "
        "Health, PCR ID 1708. Camelid monospecific IgG against Echis "
        "sochureki venom."
    ),
}


#: The curated table.
CANDIDATES: tuple[Candidate, ...] = (
    # ------------------------------------------------------------------ cobra
    Candidate(
        species=("Naja naja",),
        toxin_family="PLA2",
        toxin_target="acidic phospholipase A2",
        antibody_id="N194",
        molecule_format="scFv",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against=(
            "in vivo survival in mice challenged with Naja naja venom, up to "
            "50 per cent protection"
        ),
        evidence_level="L4",
        mrna_suitability="unknown",
        citation="Kumar2024",
        notes=(
            "Selected by phage display against Naja naja venom. No PDB entry "
            "and no public sequence located in the M4 survey."
        ),
    ),
    Candidate(
        species=("Naja naja",),
        toxin_family="unspecified",
        toxin_target="not resolved in the public record",
        antibody_id="N248",
        molecule_format="scFv",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against=(
            "in vitro binding and in vivo venom inhibition against Naja naja "
            "venom"
        ),
        evidence_level="L4",
        mrna_suitability="unknown",
        citation="Kumar2024",
        notes=(
            "Second Naja naja clone from the same screen. Its target toxin is "
            "not named in the public record, so the toxin_family is recorded "
            "as unspecified rather than guessed."
        ),
    ),
    Candidate(
        species=("Naja naja",),
        toxin_family="3FTx-L",
        toxin_target="long neurotoxin",
        antibody_id="Centi-LNX-D09",
        molecule_format="mAb-IgG",
        sequence_source=(
            "PDB:8DA1",
            "PDB:8D9Z",
            "PDB:8DA0",
            "PDB:8D9Y",
            "patent:US20240294618A1",
        ),
        patent_status="application-pending",
        tested_against=(
            "in vivo whole-venom challenge across 19 elapid species, plus "
            "recombinant long neurotoxins from cobra, krait, taipan and mamba"
        ),
        evidence_level="L4",
        mrna_suitability="adaptable",
        citation="Glanville2025",
        notes=(
            "The same clone is published as Centi-3FTX-D09 in the bioRxiv "
            "record of the discovery. The four PDB entries are Fab complexes "
            "with neurotoxins from Bungarus multicinctus, Naja nivea, "
            "Dendroaspis polylepis and Oxyuranus scutellatus."
        ),
    ),
    Candidate(
        species=("Naja naja",),
        toxin_family="3FTx-S",
        toxin_target="short neurotoxin",
        antibody_id="Centi-SNX-B03",
        molecule_format="mAb-IgG",
        sequence_source=("PDB:8V13", "patent:US20240294618A1"),
        patent_status="application-pending",
        tested_against="in vivo whole-venom challenge",
        evidence_level="L4",
        mrna_suitability="adaptable",
        citation="Glanville2025",
        notes="PDB 8V13 is a Fab complex with a short neurotoxin from Dendroaspis polylepis.",
    ),
    Candidate(
        species=("Naja kaouthia",),
        toxin_family="3FTx-L",
        toxin_target="alpha-cobratoxin",
        antibody_id="2554_01_D11",
        molecule_format="mAb-IgG",
        sequence_source=("patent:WO2023209212A1",),
        patent_status="application-pending",
        tested_against=(
            "in vitro patch-clamp neutralisation, and in vivo protection "
            "against Naja kaouthia whole venom"
        ),
        evidence_level="L4",
        mrna_suitability="unknown",
        citation="Ledsgaard2023",
        notes=(
            "Derived from 368_01_C05 by light-chain shuffling. The variable "
            "regions are disclosed in the DTU application WO2023209212A1, "
            "which names the clone 89 times, but the sequences have not been "
            "transcribed for this repository and the licence position is "
            "unresolved."
        ),
    ),
    Candidate(
        species=("Naja kaouthia",),
        toxin_family="3FTx-L",
        toxin_target="alpha-cobratoxin",
        antibody_id="368_01_C05",
        molecule_format="scFv",
        sequence_source=("patent:WO2023209212A1",),
        patent_status="application-pending",
        tested_against="in vitro neutralisation of alpha-cobratoxin",
        evidence_level="L2",
        mrna_suitability="unknown",
        citation="Ledsgaard2022",
        notes=(
            "The parental clone of the affinity-maturation campaign. It did "
            "not neutralise alpha-bungarotoxin, which is the reason its "
            "light chain was shuffled. Named in the DTU application "
            "WO2023209212A1."
        ),
    ),
    Candidate(
        species=("Naja kaouthia",),
        toxin_family="3FTx-L",
        toxin_target="alpha-cobratoxin",
        antibody_id="TPL1158_01_C09",
        molecule_format="VHH",
        sequence_source=("PDB:9GCN", "patent:WO2025078529A1"),
        patent_status="application-pending",
        tested_against="in vitro binding, by co-crystal structure",
        evidence_level="L2",
        mrna_suitability="ready",
        citation="RCSB9GCN",
        notes=(
            "The only single-domain candidate for an Indian cobra toxin whose "
            "sequence is directly available from a structural deposit. The DTU "
            "application WO2025078529A1 names the clone 71 times, so the "
            "sequence being public does not make the candidate free to build."
        ),
    ),
    Candidate(
        species=("Naja kaouthia",),
        toxin_family="3FTx-L",
        toxin_target="alpha-cobratoxin",
        antibody_id="A01",
        molecule_format="mAb-IgG",
        sequence_source=("PDB:9HUO", "PDB:9FYT", "PDB:9HXO", "patent:WO2023209212A1"),
        patent_status="application-pending",
        tested_against=(
            "in vitro receptor-blocking at pH 4.5, 5.5 and 6.0"
        ),
        evidence_level="L2",
        mrna_suitability="adaptable",
        citation="RCSB9HUO",
        notes=(
            "Selected for pH-dependent binding, which is the point of the "
            "study: venom is acidic at the bite site and near neutral in the "
            "blood. Named in the DTU application WO2023209212A1."
        ),
    ),
    Candidate(
        species=("Naja kaouthia",),
        toxin_family="3FTx-L",
        toxin_target="alpha-cobratoxin",
        antibody_id="C2, C21",
        molecule_format="VHH",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against=(
            "in vitro and in vivo neutralisation of alpha-cobratoxin in mice"
        ),
        evidence_level="L4",
        mrna_suitability="unknown",
        citation="Richard2013",
        notes=(
            "Two llama single-domain clones reported in the same study, "
            "together with a VHH-Fc construct. No public sequence located."
        ),
    ),
    Candidate(
        species=("Naja kaouthia",),
        toxin_family="3FTx-L",
        toxin_target="alpha-neurotoxin",
        antibody_id="TPL0629_01_D11",
        molecule_format="VHH",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against="in vitro neutralisation and in vivo protection",
        evidence_level="L4",
        mrna_suitability="unknown",
        citation="Laustsen2024",
        notes=(
            "Anti-long-chain alpha-neurotoxin camelid VHH. Also a component "
            "of the five-nanobody cocktail in Samanta2026."
        ),
    ),
    Candidate(
        species=("Naja kaouthia",),
        toxin_family="PLA2",
        toxin_target="phospholipase A2",
        antibody_id="TPL0637_01_A07",
        molecule_format="VHH",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against="in vitro neutralisation and in vivo protection",
        evidence_level="L4",
        mrna_suitability="unknown",
        citation="Laustsen2024",
        notes=(
            "Anti-PLA2 camelid VHH. Also a component of the five-nanobody "
            "cocktail in Samanta2026."
        ),
    ),
    Candidate(
        species=("Naja naja", "Naja kaouthia"),
        toxin_family="mixed",
        toxin_target=(
            "alpha-neurotoxins, cytotoxins and phospholipases A2 across cobra "
            "and king cobra venoms"
        ),
        antibody_id="five-VHH cocktail",
        molecule_format="VHH",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against=(
            "in vivo preincubation and rescue models against Naja naja, Naja "
            "kaouthia and both Indian king cobra species"
        ),
        evidence_level="L4",
        mrna_suitability="unknown",
        citation="Samanta2026",
        notes=(
            "The closest published thing to an Indian cobra cocktail. Only "
            "two of the five nanobodies, TPL0629_01_D11 and TPL0637_01_A07, "
            "are named in sources that could be verified during the M4 "
            "survey."
        ),
    ),
    # ------------------------------------------------------------------ krait
    Candidate(
        species=("Bungarus multicinctus",),
        toxin_family="3FTx-L",
        toxin_target="alpha-bungarotoxin",
        antibody_id="95Mat5",
        molecule_format="mAb-IgG",
        sequence_source=("PDB:8SXP",),
        patent_status="not-assessed",
        tested_against=(
            "in vivo protection against alpha-bungarotoxin and against whole "
            "venoms of Naja kaouthia and Dendroaspis polylepis"
        ),
        evidence_level="L4",
        mrna_suitability="adaptable",
        citation="Khalek2024",
        notes=(
            "PDB 8SXP is a Fab complex with a king cobra long neurotoxin, "
            "3FTx-L15. The protection data against krait alpha-bungarotoxin "
            "is in the paper."
        ),
    ),
    Candidate(
        species=("Bungarus caeruleus",),
        toxin_family="3FTx-L",
        toxin_target="alpha-delta-bungarotoxin, 3L2A_BUNCE",
        antibody_id="Centi-LNX-D09",
        molecule_format="mAb-IgG",
        sequence_source=(
            "PDB:8DA1",
            "PDB:8D9Z",
            "PDB:8DA0",
            "PDB:8D9Y",
            "patent:US20240294618A1",
        ),
        patent_status="application-pending",
        tested_against=(
            "in vivo protection against recombinant common krait long "
            "neurotoxin, and against common krait whole venom in combination "
            "with varespladib"
        ),
        evidence_level="L4",
        mrna_suitability="adaptable",
        citation="Glanville2025",
        notes=(
            "PDB 8DA1 is the Fab complex with alpha-bungarotoxin from "
            "Bungarus multicinctus. Common krait whole venom is 64.5 per "
            "cent phospholipase A2, which is why the antibody alone only "
            "doubled survival time against it."
        ),
    ),
    Candidate(
        species=("Bungarus multicinctus",),
        toxin_family="3FTx-L",
        toxin_target="alpha-bungarotoxin",
        antibody_id="2554_01_D11",
        molecule_format="mAb-IgG",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against="in vitro partial neutralisation of alpha-bungarotoxin",
        evidence_level="L2",
        mrna_suitability="unknown",
        citation="Ledsgaard2023",
        notes=(
            "The light-chain shuffling that produced this clone was what "
            "extended its cross-neutralisation into the krait toxins."
        ),
    ),
    # -------------------------------------------------------- Russell's viper
    Candidate(
        species=("Daboia russelii",),
        toxin_family="PLA2",
        toxin_target="secretory phospholipase A2, gp12B",
        antibody_id="MKSS1",
        molecule_format="VHH",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against=(
            "in vitro phospholipase inhibition, cellular toxicity and in vivo "
            "neutralisation in zebrafish"
        ),
        evidence_level="L4",
        mrna_suitability="unknown",
        citation="Li2019",
        notes=(
            "The only recombinant antibody in the survey with a named target "
            "toxin in Indian Russell's viper venom. No public sequence "
            "located. The in vivo model is zebrafish, not a mammal."
        ),
    ),
    Candidate(
        species=("Daboia russelii formosensis",),
        toxin_family="PLA2",
        toxin_target="PLA2, RV-4 or RV-7 by mass spectrometry",
        antibody_id="DRFS3, DRFS6, DRFS9, DRFS13, DRFL1",
        molecule_format="scFv",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against=(
            "in vitro binding to Daboia russelii formosensis venom proteins, "
            "and increased survival time in mice"
        ),
        evidence_level="L3",
        mrna_suitability="unknown",
        citation="Lee2017",
        notes=(
            "Chicken scFvs from egg-yolk-immunised hens. The subspecies is "
            "Taiwanese, not Indian, and the clones did not cross-react with "
            "other snake venoms. Protection was survival prolongation, not "
            "rescue."
        ),
    ),
    # ------------------------------------------------------- saw-scaled viper
    Candidate(
        species=("Echis carinatus",),
        toxin_family="PLA2",
        toxin_target="ecarpholin S, the Ser49 myotoxic phospholipase A2",
        antibody_id="E113",
        molecule_format="scFv",
        sequence_source=("PDB:8IA6",),
        patent_status="none-located",
        tested_against=(
            "inhibition of myotoxicity and haemolysis, and in vivo survival "
            "in mice challenged with Echis carinatus venom"
        ),
        evidence_level="L4",
        mrna_suitability="ready",
        citation="Kumar2024",
        notes=(
            "The construct in this repository's M3 design. The deposited "
            "structure is a crystallisation construct, so the vector residue "
            "and the hexahistidine tag are removed by script."
        ),
    ),
    Candidate(
        species=("Echis carinatus",),
        toxin_family="PLA2",
        toxin_target="ecarpholin S, the Ser49 myotoxic phospholipase A2",
        antibody_id="E10",
        molecule_format="scFv",
        sequence_source=("none",),
        patent_status="none-located",
        tested_against=(
            "inhibition of myotoxicity and haemolysis, and in vivo survival "
            "in mice challenged with Echis carinatus venom"
        ),
        evidence_level="L4",
        mrna_suitability="unknown",
        citation="Kumar2024",
        notes=(
            "Together with E113 the combination gave 33 per cent survival, "
            "against 16 per cent for an equal amount of either clone alone. "
            "No public sequence located."
        ),
    ),
    Candidate(
        species=("Echis carinatus",),
        toxin_family="SVMP",
        toxin_target="ecarin, a zinc metalloproteinase-disintegrin",
        antibody_id="H11",
        molecule_format="Fab",
        sequence_source=("PDB:9CLP",),
        patent_status="not-assessed",
        tested_against="in vitro binding, by co-crystal structure",
        evidence_level="L2",
        mrna_suitability="adaptable",
        citation="RCSB9CLP",
        notes=(
            "The only structural candidate found against a snake venom "
            "metalloproteinase. The deposited toxin comes from the Kenyan "
            "form of the saw-scaled viper, so the epitope has to be checked "
            "against Indian E. carinatus before the antibody can be used "
            "here."
        ),
    ),
    Candidate(
        species=("Echis carinatus",),
        toxin_family="unspecified",
        toxin_target="unspecified venom components",
        antibody_id="human recombinant polyclonal Fab",
        molecule_format="polyclonal-Fab",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against="in vivo neutralisation in laboratory mice",
        evidence_level="L4",
        mrna_suitability="blocked",
        citation="IMPACT1864",
        notes=(
            "Produced in bacteria from a human gene library. There is no "
            "defined sequence to transcribe, so this candidate cannot become "
            "an mRNA construct without first being turned into a monoclonal."
        ),
    ),
    Candidate(
        species=("Echis carinatus sochureki",),
        toxin_family="unspecified",
        toxin_target="unspecified venom components",
        antibody_id="camelid monospecific IgG",
        molecule_format="polyclonal-IgG",
        sequence_source=("none",),
        patent_status="not-assessed",
        tested_against=(
            "preclinical in vivo neutralisation of lethal, haemorrhagic, "
            "coagulant and local necrotising activities"
        ),
        evidence_level="L4",
        mrna_suitability="blocked",
        citation="IMPACT1708",
        notes=(
            "Raised in dromedaries against Indian E. carinatus sochureki "
            "venom. Useful as a benchmark for what the existing polyclonal "
            "products achieve, and not a template for a defined mRNA "
            "construct."
        ),
    ),
)
