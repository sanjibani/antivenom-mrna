"""Tests for the M6 single-transcript design.

The transcript carries two payloads on one open reading frame, joined by
the P2A ribosome-skipping element. The tests protect three things: the P2A
peptide is the same molecule in two deposited constructs, the two mature
products the cell makes are checked by translation against the published
molecules, and the committed design files cannot drift away from the code
that generates them.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from mrna.fasta import read_fasta
from mrna.p2a import GSG_P2A_LINKER, P2A, P2A_DOWNSTREAM_PROLINE
from mrna.payload import load_payloads

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"
RAW = REFERENCE / "raw"
DESIGNS = ROOT / "data" / "designs"

_SPEC = importlib.util.spec_from_file_location(
    "cocktail_transcript",
    ROOT / "scripts" / "cocktail_transcript.py",
)
cocktail_transcript = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(cocktail_transcript)


@pytest.fixture(scope="module")
def design() -> tuple[object, object, dict[str, str]]:
    """The built construct, its optimisation result, and its products."""
    from mrna.reference import load_codon_usage

    usage = load_codon_usage(REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv")
    return cocktail_transcript.build_construct(usage)


# ------------------------------------------------------------------- the 2A


def test_the_p2a_peptide_is_the_same_in_two_deposited_constructs() -> None:
    """The M0 two-source rule, applied to the new element.

    Both raw entries are read, and each must carry the peptide the engine
    constant holds, identical over all 19 residues.
    """
    found = cocktail_transcript.read_p2a_from_raw()
    assert len(found) == 2
    for filename, peptide in found.items():
        assert peptide == P2A.sequence
        assert Path(filename).name in {p.name for p in RAW.glob("*.fasta")}


def test_the_skip_site_is_between_the_last_glycine_and_the_proline() -> None:
    """The element's two products follow from where the bond is not made."""
    assert P2A.sequence.endswith("GP")
    assert P2A.upstream_retained == P2A.sequence[:-1]
    assert P2A.downstream_proline == "P"
    assert P2A_DOWNSTREAM_PROLINE == "P"


def test_the_linker_is_the_spacer_and_the_peptide_up_to_the_glycine() -> None:
    """The upstream product carries the spacer plus 18 residues, 21 in all."""
    assert len(GSG_P2A_LINKER) == 21
    assert GSG_P2A_LINKER == "GSG" + P2A.sequence[:-1]
    assert GSG_P2A_LINKER.endswith("ENPG")


# --------------------------------------------------------------- the design


def test_the_construct_length_and_orf(design) -> None:
    """One open reading frame of 420 residues on a 1,722-nt message."""
    construct, optimization, _products = design
    assert construct.length == 1722
    assert len(optimization.protein) == 420


def test_the_orf_translates_to_the_designed_protein(design) -> None:
    """Translation, not codon counting, decides the protein."""
    construct, optimization, _products = design
    from mrna.sequence import translate

    assert translate(construct.open_reading_frame().cds).rstrip("*") == (
        optimization.protein
    )


def test_the_upstream_mature_product_is_the_nanobody_plus_its_scar(design) -> None:
    """The published 123 residues, then the 21-residue GSG-P2A scar."""
    payloads = {p.name: p for p in load_payloads()}
    _construct, _optimization, products = design
    expected = payloads["TPL1158_01_C09"].sequence + GSG_P2A_LINKER
    assert products["upstream_mature"] == expected
    assert len(products["upstream_mature"]) == 144


def test_the_downstream_mature_product_is_the_published_e113(design) -> None:
    """The molecule the 2024 paper measured, residue for residue.

    The check is against the committed reference file rather than against
    the payload loader, so a regression in the loader would be caught by
    this test being wrong about the file rather than agreeing with it.
    """
    e113 = read_fasta(REFERENCE / "e113_scfv.fasta")[0][1]
    _construct, _optimization, products = design
    assert products["downstream_mature"] == e113
    assert len(products["downstream_mature"]) == 239


def test_no_forbidden_motif_in_any_coding_element(design) -> None:
    """The polyadenylation signal stays out of every coding element."""
    construct, _optimization, _products = design
    coding = set(cocktail_transcript.ORF_ELEMENTS)
    hits = construct.find_motifs([])
    from mrna.motifs import DEFAULT_MOTIFS

    hits = construct.find_motifs(DEFAULT_MOTIFS)
    assert not [hit for hit in hits if hit.element in coding]


def test_the_element_table_names_all_nine_elements(design) -> None:
    """The regulatory-format table, with the two leaders and the skip."""
    construct, _optimization, _products = design
    names = [element.name for element in construct.elements]
    assert names == [
        "cap",
        "5'-UTR",
        "signal peptide",
        "TPL1158_01_C09 VHH",
        "P2A",
        "signal peptide 2",
        "E113 scFv",
        "3'-UTR",
        "poly(A)",
    ]


# ---------------------------------------------------------------- artifacts


def test_the_committed_report_matches_the_generated_one() -> None:
    """The published report cannot drift away from the code that generates it."""
    report_path = ROOT / "docs" / "m6_single_transcript.md"
    if not report_path.exists():
        pytest.skip("docs/m6_single_transcript.md has not been generated yet")
    assert report_path.read_text(encoding="utf-8") == (
        cocktail_transcript.build_report()
    )


def test_the_committed_fasta_matches_the_generated_design(design) -> None:
    """The committed design sequence must be the one the engine builds."""
    fasta_path = cocktail_transcript.FASTA_PATH
    if not fasta_path.exists():
        pytest.skip("the design FASTA has not been generated yet")
    construct, _optimization, _products = design
    header, sequence = read_fasta(fasta_path)[0]
    assert sequence == construct.sequence
    assert "no therapeutic efficacy demonstrated" in header


def test_the_committed_element_table_matches_the_construct(design) -> None:
    """The element table cannot drift away from the construct it describes."""
    table_path = cocktail_transcript.TABLE_PATH
    if not table_path.exists():
        pytest.skip("the element table has not been generated yet")
    construct, _optimization, _products = design
    rows = table_path.read_text(encoding="utf-8").strip().split("\n")
    assert rows[0] == "element\tstart\tend\tlength\tdescription"
    assert len(rows) == len(construct.elements) + 1
    for row, element in zip(rows[1:], construct.elements, strict=True):
        fields = row.split("\t")
        assert fields[0] == element.name
        assert int(fields[1]) == element.start
        assert int(fields[2]) == element.end
