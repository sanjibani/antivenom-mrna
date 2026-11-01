"""Turn the raw downloads in ``data/reference/raw`` into committed reference files.

Every reference file in this repository is derived by this script from a raw
download, so that any reader can regenerate the derived file and diff it. The
script is idempotent: running it twice produces byte-identical output.

Usage::

    uv run python scripts/prepare_reference_data.py

Sources and their licences are recorded in ``data/PROVENANCE.md``.
"""

from __future__ import annotations

import hashlib
import re
from pathlib import Path

from mrna.fasta import read_fasta, write_fasta

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "reference" / "raw"
OUT = ROOT / "data" / "reference"

CODON_ROW = re.compile(
    r"(?P<codon>[ACGTU]{3})\s+(?P<aa>[A-Z*])\s+(?P<fraction>\d+\.\d+)\s+"
    r"(?P<per_thousand>\d+\.\d+)\s+\(\s*(?P<count>\d+)\)",
)

#: A line of the WHO document's sequence block: groups of nucleotides in groups
#: of ten, followed by the running position of the last base on that line.
WHO_SEQUENCE_LINE = re.compile(r"^(?P<bases>[ACGTUΨ\s]+?)\s+(?P<position>\d+)\s*$")

#: The feature table of the WHO INN document 11889 for Comirnaty, as
#: published. Positions are 1-based and inclusive.
WHO_BNT162B2_LAYOUT: tuple[tuple[str, int, int, str], ...] = (
    ("cap", 1, 2, "Modified 5'-cap1 structure (m7G plus triphosphate)"),
    ("5'-UTR", 3, 54, "Derived from human alpha-globin RNA, with an optimised Kozak sequence"),
    ("signal peptide", 55, 102, "Extended signal sequence of the S glycoprotein"),
    ("S protein", 103, 3879, "Codon-optimised spike with K986P and V987P, then stop codons"),
    ("3'-UTR", 3880, 4174, "AES mRNA and mitochondrial 12S rRNA elements"),
    ("poly(A)", 4175, 4284, "30 adenosines, a 10-nucleotide linker, then 70 adenosines"),
)


def parse_who_sequence(path: Path) -> str:
    """Extract the full mRNA sequence from the WHO INN document text.

    The document prints the sequence in groups of ten bases with a running
    position number at the end of each line. Every line is validated against
    that running position, so a line of prose can never be mistaken for
    sequence: the parser only accepts a line whose nucleotides bring the
    total length exactly to the number printed on it.

    N1-methylpseudouridine is written as ``Psi`` in the document and becomes
    ``T`` here, because the engine works in the DNA alphabet. The
    modification is recorded in the element table, not in the string.

    Args:
        path: Path to the plain-text conversion of the document.

    Returns:
        The full nucleotide sequence, 5' to 3'.

    Raises:
        ValueError: If the sequence could not be assembled consistently.
    """
    pieces: list[str] = []
    expected = 0
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = WHO_SEQUENCE_LINE.match(raw_line.strip())
        if match is None:
            continue
        bases = "".join(match.group("bases").split()).replace("Ψ", "U")
        position = int(match.group("position"))
        if len(bases) != position - expected:
            continue
        if set(bases) - set("ACGTU"):
            continue
        pieces.append(bases)
        expected = position
    if not pieces:
        msg = f"{path}: no sequence lines found in the document text"
        raise ValueError(msg)
    return "".join(pieces)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file, for provenance records."""
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def parse_kazusa(path: Path) -> list[tuple[str, str, float, float, int]]:
    """Parse the Kazusa Codon Usage Database HTML page into rows.

    The page embeds the table in a ``<PRE>`` block as fixed-width text, with
    the RNA alphabet. Rows look like::

        UUU F 0.46 17.6 (714298)

    Args:
        path: Path to the downloaded HTML file.

    Returns:
        A list of ``(codon, amino_acid, fraction, per_thousand, count)`` rows,
        with DNA codon letters.

    Raises:
        ValueError: If no codon rows could be parsed.
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    rows: list[tuple[str, str, float, float, int]] = []
    for match in CODON_ROW.finditer(text):
        codon = match.group("codon").replace("U", "T")
        rows.append(
            (
                codon,
                match.group("aa"),
                float(match.group("fraction")),
                float(match.group("per_thousand")),
                int(match.group("count")),
            ),
        )
    if not rows:
        msg = f"{path}: no codon usage rows found, the page layout must have changed"
        raise ValueError(msg)
    return rows


def write_codon_usage(rows: list[tuple[str, str, float, float, int]]) -> Path:
    """Write the parsed codon usage table to a TSV file.

    Args:
        rows: Rows as returned by :func:`parse_kazusa`.

    Returns:
        The path of the written file.
    """
    target = OUT / "homo_sapiens_codon_usage_kazusa.tsv"
    lines = ["# Homo sapiens codon usage table, Kazusa Codon Usage Database"]
    lines.append("# Source: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606&aa=1&style=N")
    lines.append("# Reference: Nakamura Y, Gojobori T, Ikemura T. 2000. Codon usage")
    lines.append("# tabulated from international DNA sequence databases: status for the")
    lines.append("# year 2000. Nucleic Acids Res 28(1):292.")
    lines.append("# Fields: codon(T) amino_acid fraction fraction_per_thousand count")
    for codon, aa, fraction, per_thousand, count in rows:
        lines.append(f"{codon}\t{aa}\t{fraction}\t{per_thousand}\t{count}")
    target.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return target


def write_bnt162b2_cds() -> Path:
    """Extract the BNT162b2 coding sequence from the WHO-document transcription."""
    records = read_fasta(RAW / "berthubert_vaccine-s.fasta")
    if not records:
        msg = "berthubert_vaccine-s.fasta contains no records"
        raise ValueError(msg)
    _, sequence = records[0]
    header = (
        "BNT162b2 (Comirnaty) coding sequence | source: WHO Comirnaty vaccine description "
        "document, transcribed at https://github.com/berthubert/bnt162b2 (file vaccine-s.fasta) "
        "| length 3822 nt = 1273 codons plus stop"
    )
    target = OUT / "bnt162b2_cds.fasta"
    write_fasta(target, [(header, sequence)])
    return target


def write_reference_gfp() -> Path:
    """Normalise the downloaded UniProt GFP record."""
    records = read_fasta(RAW / "uniprot_P42212_gfp.fasta")
    if not records:
        msg = "uniprot_P42212_gfp.fasta contains no records"
        raise ValueError(msg)
    header, sequence = records[0]
    target = OUT / "gfp_P42212.fasta"
    write_fasta(target, [(header, sequence)])
    return target


def write_reference_protein() -> Path:
    """Normalise the downloaded UniProt SARS-CoV-2 spike record."""
    records = read_fasta(RAW / "uniprot_P0DTC2_spike.fasta")
    if not records:
        msg = "uniprot_P0DTC2_spike.fasta contains no records"
        raise ValueError(msg)
    header, sequence = records[0]
    target = OUT / "sars_cov2_spike_P0DTC2.fasta"
    write_fasta(target, [(header, sequence)])
    return target


def write_bnt162b2_mrna(sequence: str) -> Path:
    """Write the full Comirnaty mRNA sequence taken from the WHO document.

    Args:
        sequence: The full nucleotide sequence.

    Returns:
        The path of the written file.
    """
    header = (
        "BNT162b2 (Comirnaty) full mRNA, 4284 nt | source: WHO International "
        "Nonproprietary Names Programme document 11889 (Comirnaty), archived at "
        "https://berthub.eu/articles/11889.doc | N1-methylpseudouridine written as T "
        "(DNA alphabet, as in GenBank mRNA records)"
    )
    target = OUT / "bnt162b2_mrna_who.fasta"
    write_fasta(target, [(header, sequence.replace("U", "T"))])
    return target


def main() -> None:
    """Regenerate every derived reference file."""
    OUT.mkdir(parents=True, exist_ok=True)

    usage_path = write_codon_usage(parse_kazusa(RAW / "kazusa_homo_sapiens_2020.txt"))
    print(f"wrote {usage_path.relative_to(ROOT)}")

    who_sequence = parse_who_sequence(RAW / "who_inn_11889_bnt162b2.txt")
    print(f"WHO document sequence: {len(who_sequence)} nt")
    mrna_path = write_bnt162b2_mrna(who_sequence)
    print(f"wrote {mrna_path.relative_to(ROOT)}")

    for path in (write_bnt162b2_cds(), write_reference_protein(), write_reference_gfp()):
        print(f"wrote {path.relative_to(ROOT)}")

    print("\nsha256 of raw downloads:")
    for raw in sorted(RAW.glob("*")):
        print(f"  {raw.name}: {sha256(raw)}")


if __name__ == "__main__":
    main()
