"""Minimal FASTA reading and writing.

The engine reads reference sequences from plain FASTA files so that every input
carries its own provenance in its header. This module is deliberately tiny and
dependency-free: the alternative, a full bioinformatics toolkit, would make the
engine harder to audit for no benefit at this scale.
"""

from __future__ import annotations

from pathlib import Path

__all__ = ["read_fasta", "write_fasta"]


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Read a FASTA file into ``(header, sequence)`` pairs.

    Sequence lines are concatenated and whitespace is removed. Header text is
    returned verbatim, minus the leading ``>``, so that provenance information
    written into the header by hand or by a download script survives.

    Args:
        path: Path to the FASTA file.

    Returns:
        A list of ``(header, sequence)`` tuples in file order.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If a sequence line appears before any header.
    """
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []

    text = Path(path).read_text(encoding="utf-8")
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header = line[1:].strip()
            chunks = []
        else:
            if header is None:
                msg = f"{path}: found sequence data before any '>' header line"
                raise ValueError(msg)
            chunks.append(line)

    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def write_fasta(path: str | Path, records: list[tuple[str, str]], *, width: int = 60) -> None:
    """Write ``(header, sequence)`` pairs to a FASTA file.

    Args:
        path: Destination path.
        records: The records to write.
        width: Number of sequence characters per line.

    Raises:
        ValueError: If a record has an empty sequence.
    """
    lines: list[str] = []
    for header, seq in records:
        if not seq:
            msg = f"refusing to write empty sequence for header {header!r}"
            raise ValueError(msg)
        lines.append(f">{header}")
        lines.extend(seq[i : i + width] for i in range(0, len(seq), width))

    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
