"""Loaders for the committed reference data.

The engine reads its inputs from plain files under ``data/reference/`` so that
every input carries its own provenance in its file header and its own
checksum. These helpers are the one place that knows the shape of those
files, so that a change in the file format breaks in one place rather than in
every script that reads one.

The reference data is part of the repository rather than part of the
installed package, because it is evidence: a downloaded file with a recorded
checksum and a stated origin can be checked by anyone who clones the
repository. An installed wheel carries the code, not the evidence.
"""

from __future__ import annotations

from pathlib import Path

from mrna.fasta import read_fasta

__all__ = ["load_codon_usage", "load_single_sequence", "reference_dir"]


def reference_dir() -> Path:
    """Return the path of the committed reference-data directory.

    Returns:
        The ``data/reference`` directory of the repository.

    Raises:
        FileNotFoundError: If the directory is not present, which happens
            when the package is installed without the repository.
    """
    path = Path(__file__).resolve().parent.parent.parent / "data" / "reference"
    if not path.is_dir():
        msg = (
            f"{path} is not available. The reference data lives in the "
            "repository, not in the installed package. Run from a clone of "
            "the repository."
        )
        raise FileNotFoundError(msg)
    return path


def load_codon_usage(path: str | Path) -> dict[str, float]:
    """Load a Kazusa-style codon usage table into a codon-to-fraction mapping.

    The table is written by ``scripts/prepare_reference_data.py`` with one
    row per codon: the codon, the amino acid it encodes, and the fraction of
    that amino acid's usage it accounts for. Stop codons are skipped, because
    the Codon Adaptation Index is defined over sense codons only.

    Args:
        path: Path to the tab-separated table.

    Returns:
        Mapping of codon to synonymous usage fraction.

    Raises:
        ValueError: If the file contains no codon rows.
    """
    table: dict[str, float] = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3 or fields[1] == "*":
            continue
        table[fields[0]] = float(fields[2])
    if not table:
        msg = f"{path}: no codon usage rows found"
        raise ValueError(msg)
    return table


def load_single_sequence(path: str | Path) -> str:
    """Load the first record of a FASTA file as a plain sequence.

    Args:
        path: Path to the FASTA file.

    Returns:
        The sequence string.

    Raises:
        ValueError: If the file has no records.
    """
    records = read_fasta(path)
    if not records:
        msg = f"{path}: no FASTA records found"
        raise ValueError(msg)
    return records[0][1]
