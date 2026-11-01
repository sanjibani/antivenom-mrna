"""Tests for the first venom construct.

These tests pin the things a reader would want to check about a design
document: that every sequence is the one its source says it is, that the
construct encodes what it claims to encode, and that the committed design
files match the design the engine produces.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from mrna.assemble import A110_POLYA, CAP1, Construct, full_construct
from mrna.fasta import read_fasta
from mrna.leader import HSA_SIGNAL_PEPTIDE
from mrna.motifs import DEFAULT_MOTIFS
from mrna.optimize import optimize_cds
from mrna.policy import HUMAN_POLICY
from mrna.reference import load_codon_usage, load_single_sequence
from mrna.sequence import translate
from mrna.utr import AES_MTRNR1_3UTR, ALPHA_GLOBIN_5UTR

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"
RAW = REFERENCE / "raw"
DESIGNS = ROOT / "data" / "designs"

STOP_CODON = "TGA"


def payload() -> str:
    """Return the committed E113 single-chain variable fragment."""
    return load_single_sequence(REFERENCE / "e113_scfv.fasta").strip()


def toxin() -> str:
    """Return the committed ecarpholin S sequence."""
    return load_single_sequence(REFERENCE / "ecarpholin_S_P48650.fasta").strip()


def build() -> Construct:
    """Build the venom construct the way the script does."""
    usage_table = load_codon_usage(REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv")
    protein = HSA_SIGNAL_PEPTIDE.sequence + payload()
    design = optimize_cds(protein, usage_table, HUMAN_POLICY)
    return full_construct(
        "E113 scFv against Echis carinatus ecarpholin S",
        utr5=ALPHA_GLOBIN_5UTR,
        coding=design.cds + STOP_CODON,
        utr3=AES_MTRNR1_3UTR,
        poly_a=A110_POLYA,
        leader=HSA_SIGNAL_PEPTIDE,
        coding_name="E113 scFv",
    )


def test_the_payload_is_the_deposited_construct_without_its_artefacts() -> None:
    """The payload must be the deposited antibody minus vector residue and tag."""
    deposited = read_fasta(RAW / "rcsb_8IA6_scfv_e113.fasta")[0][1]
    assert deposited.startswith("A")  # vector residue
    assert deposited.endswith("RAAALEHHHHHH")  # vector tag
    assert payload() == deposited[1:-12]
    assert len(payload()) == 239


def test_the_payload_starts_and_ends_at_the_antibody() -> None:
    """The trimmed payload must start at the germline and end at the constant region."""
    assert payload().startswith("EVQLLESGGG")
    assert payload().endswith("FGQGTKVEIK")


def test_the_toxin_matches_its_crystal_structure() -> None:
    """The UniProt sequence and the PDB structure must be the same molecule.

    This is the cross-check that M0 established as the project's standard:
    one source is not enough.
    """
    structure = read_fasta(RAW / "rcsb_2QHE_ecarpholin_S.fasta")[0][1]
    assert toxin() == structure
    assert len(toxin()) == 122


def test_the_toxin_has_the_seven_annotated_disulfide_bonds() -> None:
    """UniProt annotates seven disulfide bonds for this sequence."""
    assert toxin().count("C") == 14


def test_the_signal_peptide_is_the_albumin_leader() -> None:
    """The leader must be the first 18 residues of the albumin precursor."""
    albumin = load_single_sequence(REFERENCE / "albumin_P02768.fasta")
    assert albumin.startswith(HSA_SIGNAL_PEPTIDE.sequence)
    assert len(HSA_SIGNAL_PEPTIDE.sequence) == 18


def test_the_construct_encodes_the_leader_and_the_payload() -> None:
    """The coded protein must be the signal peptide followed by the antibody."""
    construct = build()
    orf = construct.open_reading_frame()
    protein = translate(orf.cds).rstrip("*")
    assert protein == HSA_SIGNAL_PEPTIDE.sequence + payload()
    assert len(protein) == 257


def test_the_element_table_has_the_expected_coordinates() -> None:
    """Every element must sit where the design document says it does."""
    construct = build()
    observed = tuple(
        (element.name, element.start, element.end) for element in construct.elements
    )
    assert observed == (
        ("cap", 1, 2),
        ("5'-UTR", 3, 54),
        ("sig", 55, 108),
        ("E113 scFv", 109, 828),
        ("3'-UTR", 829, 1123),
        ("poly(A)", 1124, 1233),
    )
    assert construct.length == 1233


def test_the_coding_region_carries_no_forbidden_motif() -> None:
    """The coding region must be free of the polyadenylation signal."""
    construct = build()
    hits = construct.find_motifs(DEFAULT_MOTIFS)
    assert [hit for hit in hits if hit.element in {"sig", "E113 scFv"}] == []


def test_there_is_no_upstream_start_codon() -> None:
    """A start codon in the 5' UTR would compete with the intended one."""
    assert build().upstream_starts() == ()


def test_the_design_is_deterministic() -> None:
    """The same inputs must produce the same construct every time."""
    assert build().sequence == build().sequence


def test_the_committed_design_files_match_the_engine() -> None:
    """The committed design must be the design the engine produces."""
    records = read_fasta(DESIGNS / "e113_scfv_ecarpholin_s_mrna.fasta")
    assert len(records) == 3
    construct = build()
    assert records[0][1] == construct.sequence
    assert records[1][1] == HSA_SIGNAL_PEPTIDE.sequence + payload()
    assert records[2][1] == toxin()


def test_the_committed_element_table_matches_the_engine() -> None:
    """The committed element table must be the one the engine computes."""
    rows = (DESIGNS / "e113_scfv_ecarpholin_s_mrna.elements.tsv").read_text().splitlines()
    assert rows[0].split("\t") == ["element", "start", "end", "length", "description"]
    construct = build()
    assert len(rows) == len(construct.elements) + 1
    for row, element in zip(rows[1:], construct.elements, strict=True):
        fields = row.split("\t")
        assert fields[0] == element.name
        assert int(fields[1]) == element.start
        assert int(fields[2]) == element.end
        assert int(fields[3]) == element.length


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("cap", 2),
        ("5'-UTR", 52),
        ("sig", 54),
        ("E113 scFv", 720),
        ("3'-UTR", 295),
        ("poly(A)", 110),
    ],
)
def test_each_element_has_its_documented_length(name: str, expected: int) -> None:
    """Each element must be as long as the design document says."""
    construct = build()
    element = next(item for item in construct.elements if item.name == name)
    assert element.length == expected


def test_the_cap_is_the_same_one_the_licensed_product_uses() -> None:
    """The construct must start with the same cap footprint as Comirnaty."""
    construct = build()
    assert construct.part("cap").sequence == CAP1.prefix
