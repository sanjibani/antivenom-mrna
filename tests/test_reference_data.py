"""Tests for the committed reference data.

These tests protect the inputs to every other test. If a reference file is
regenerated wrongly, or a download silently changes upstream, these fail
first.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from mrna.anatomy import LayoutElement, verify_layout
from mrna.fasta import read_fasta

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "data" / "reference"
RAW = REFERENCE / "raw"

PUBLISHED_LAYOUT: tuple[LayoutElement, ...] = (
    LayoutElement("cap", 1, 2, "Modified 5'-cap1"),
    LayoutElement("5'-UTR", 3, 54, "Human alpha-globin derived"),
    LayoutElement("signal peptide", 55, 102, "Extended leader"),
    LayoutElement("S protein", 103, 3879, "Codon-optimised spike then stop codons"),
    LayoutElement("3'-UTR", 3880, 4174, "AES plus 12S rRNA"),
    LayoutElement("poly(A)", 4175, 4284, "A30 plus linker plus A70"),
)


def load_prepare_script():
    """Import the reference preparation script as a module.

    Returns:
        The imported module.
    """
    path = ROOT / "scripts" / "prepare_reference_data.py"
    spec = importlib.util.spec_from_file_location("prepare_reference_data", path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def who_mrna() -> str:
    """Return the full Comirnaty mRNA parsed from the WHO document."""
    return read_fasta(REFERENCE / "bnt162b2_mrna_who.fasta")[0][1]


def test_the_who_document_sequence_is_4284_nt() -> None:
    """The published construct length must be reproduced from the document."""
    assert len(who_mrna()) == 4284


def test_every_published_coordinate_is_reproduced() -> None:
    """Each element of the published feature table must slice cleanly."""
    checks = verify_layout(who_mrna(), PUBLISHED_LAYOUT)
    assert all(check.matches for check in checks)
    assert [check.observed_length for check in checks] == [2, 52, 48, 3777, 295, 110]


def test_the_poly_a_tail_has_its_published_structure() -> None:
    """30 adenosines, the 10-nt linker, then 70 adenosines."""
    poly_a = who_mrna()[4174:4284]
    assert poly_a[:30] == "A" * 30
    assert poly_a[30:40] == "GCATATGACT"
    assert poly_a[40:] == "A" * 70


def test_the_two_stop_codons_sit_at_3874_to_3879() -> None:
    """The document places the stop codons at 3874 to 3879."""
    assert who_mrna()[3873:3879] == "TGATGA"


def test_the_signal_peptide_starts_where_the_document_says() -> None:
    """The 48-nt signal peptide must translate to the extended leader."""
    from mrna.sequence import translate

    assert translate(who_mrna()[54:102]) == "MFVFLVLLPLVSSQCV"


def test_the_committed_mrna_file_matches_a_fresh_parse_of_the_document() -> None:
    """Re-parsing the raw document must reproduce the committed file exactly."""
    module = load_prepare_script()
    parsed = module.parse_who_sequence(RAW / "who_inn_11889_bnt162b2.txt")
    assert parsed.replace("U", "T") == who_mrna()


def test_the_read_derived_assembly_matches_the_document_base_for_base() -> None:
    """The independent sequencing assembly must agree with the primary source."""
    assembly = read_fasta(RAW / "naalytics_figure1figure2_032321.fasta")[0][1]
    document = who_mrna().replace("U", "T")
    overlap = min(len(document), len(assembly))
    mismatches = [
        index for index in range(overlap) if document[index] != assembly[index]
    ]
    assert mismatches == []
    assert len(assembly) == 4175


def test_the_codon_usage_table_is_complete() -> None:
    """The Kazusa table must cover all 61 sense codons."""
    rows = [
        line.split("\t")
        for line in (REFERENCE / "homo_sapiens_codon_usage_kazusa.tsv")
        .read_text()
        .splitlines()
        if line and not line.startswith("#")
    ]
    sense = [row for row in rows if row[1] != "*"]
    assert len(sense) == 61
    assert len({row[0] for row in sense}) == 61


def test_the_kazusa_page_still_parses() -> None:
    """The downloaded Kazusa page must still yield a usable table."""
    module = load_prepare_script()
    rows = module.parse_kazusa(RAW / "kazusa_homo_sapiens_2020.txt")
    assert len(rows) == 64
    assert sum(1 for row in rows if row[1] == "*") == 3


def test_raw_files_have_not_changed() -> None:
    """The raw downloads must match the checksums recorded at download time."""
    module = load_prepare_script()
    recorded = {}
    for line in (REFERENCE / "CHECKSUMS.sha256").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        digest, name = line.split()
        recorded[name] = digest
    assert recorded, "no checksums recorded"
    for name, digest in recorded.items():
        path = RAW / name
        assert path.exists(), f"missing raw file {name}"
        assert module.sha256(path) == digest, f"{name} changed since download"


def test_every_raw_file_has_a_recorded_checksum() -> None:
    """A raw download without a recorded checksum is unverified.

    Without this check a new raw file could be committed and never added to
    the checksum file, and the test above would pass without ever looking at
    it.
    """
    recorded = set()
    for line in (REFERENCE / "CHECKSUMS.sha256").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        _digest, name = line.split()
        recorded.add(name)
    present = {
        path.name
        for path in RAW.iterdir()
        if path.is_file() and path.name != ".DS_Store"
    }
    assert present == recorded, (
        f"raw files without a recorded checksum: {sorted(present - recorded)}; "
        f"checksums recorded for files that are absent: {sorted(recorded - present)}"
    )


def test_the_prepared_cds_file_is_reproducible() -> None:
    """Re-running the preparation must reproduce the committed CDS file."""
    load_prepare_script()
    records = read_fasta(RAW / "berthubert_vaccine-s.fasta")
    committed = read_fasta(REFERENCE / "bnt162b2_cds.fasta")[0][1]
    assert records[0][1] == committed


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("cap", 2),
        ("5'-UTR", 52),
        ("signal peptide", 48),
        ("S protein", 3777),
        ("3'-UTR", 295),
        ("poly(A)", 110),
    ],
)
def test_layout_elements_have_their_published_lengths(name: str, expected: int) -> None:
    """Each published element must have the length the document states."""
    element = next(item for item in PUBLISHED_LAYOUT if item.name == name)
    assert element.length == expected
