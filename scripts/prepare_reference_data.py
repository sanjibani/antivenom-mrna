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


#: The C-terminal expression tag of the crystallised E113 construct: the
#: vector's own residues followed by a hexahistidine tag. It is an
#: expression and purification handle, not part of the antibody, so it is
#: not encoded into the designed construct.
#:
#: Source: the deposited construct of PDB entry 8IA6, entity 1.
EXPRESSION_TAG = "RAAALEHHHHHH"

#: The first residue of the deposited E113 construct is an alanine that is
#: not part of the human IGHV3-23 open reading frame, and it was not located
#: in the crystal. It is a vector artefact and is removed with the tag.
#:
#: Source: the deposited construct of PDB entry 8IA6, entity 1, compared
#: with the human IGHV3-23 germline sequence.
VECTOR_N_TERMINUS = "A"

#: The two ends of the antibody proper, used to check the trim.
#:
#: ``EVQLLESGGG`` is the start of the human IGHV3-23 germline open reading
#: frame and ``FGQGTKVEIK`` is the end of the human IGKC constant region.
ANTIBODY_N_TERMINUS = "EVQLLESGGG"
ANTIBODY_C_TERMINUS = "FGQGTKVEIK"


def write_e113_payload() -> Path:
    """Extract the E113 single-chain variable fragment from the PDB construct.

    The deposited construct is a crystallisation construct: the antibody,
    with a vector residue at the N-terminus and a hexahistidine expression
    tag at the C-terminus. This function removes both and writes the
    antibody itself, which is what a therapeutic construct encodes.

    Args:
        None.

    Returns:
        The path of the written file.

    Raises:
        ValueError: If the deposited construct does not have the expected
            ends, which would mean the trim has to be revisited.
    """
    records = read_fasta(RAW / "rcsb_8IA6_scfv_e113.fasta")
    if not records:
        msg = "rcsb_8IA6_scfv_e113.fasta contains no records"
        raise ValueError(msg)
    _, deposited = records[0]

    if not deposited.startswith(VECTOR_N_TERMINUS):
        msg = "the deposited E113 construct does not start with the expected vector residue"
        raise ValueError(msg)
    if not deposited.endswith(EXPRESSION_TAG):
        msg = "the deposited E113 construct does not end with the expected expression tag"
        raise ValueError(msg)

    payload = deposited[len(VECTOR_N_TERMINUS) : -len(EXPRESSION_TAG)]
    if not payload.startswith(ANTIBODY_N_TERMINUS):
        msg = "the trimmed E113 payload does not start with the IGHV3-23 germline sequence"
        raise ValueError(msg)
    if not payload.endswith(ANTIBODY_C_TERMINUS):
        msg = "the trimmed E113 payload does not end with the IGKC constant region"
        raise ValueError(msg)

    header = (
        "E113 single-chain variable fragment, anti-Echis carinatus PLA2 | "
        "source: PDB entry 8IA6, entity 1, deposited by Kumar A, Madni ZK, "
        "Salunke DM | the deposited construct minus the N-terminal vector "
        "residue and the C-terminal expression tag | length 239 aa = "
        "VH (117 aa) + linker (15 aa) + VL (107 aa)"
    )
    target = OUT / "e113_scfv.fasta"
    write_fasta(target, [(header, payload)])
    return target


def write_ecarpholin_s() -> Path:
    """Normalise the downloaded UniProt ecarpholin S record.

    Ecarpholin S is the Lys49 myotoxic phospholipase A2 of *Echis
    carinatus*, the saw-scaled viper of India's Big Four, and it is the
    antigen that the E113 scFv binds.
    """
    records = read_fasta(RAW / "uniprot_P48650_ecarpholin_S.fasta")
    if not records:
        msg = "uniprot_P48650_ecarpholin_S.fasta contains no records"
        raise ValueError(msg)
    header, sequence = records[0]
    target = OUT / "ecarpholin_S_P48650.fasta"
    write_fasta(target, [(header, sequence)])
    return target


def write_albumin() -> Path:
    """Normalise the downloaded UniProt human serum albumin record.

    The signal peptide of this precursor is the leader used by the venom
    constructs.
    """
    records = read_fasta(RAW / "uniprot_P02768_albumin.fasta")
    if not records:
        msg = "uniprot_P02768_albumin.fasta contains no records"
        raise ValueError(msg)
    header, sequence = records[0]
    target = OUT / "albumin_P02768.fasta"
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


def split_rcsb_records(records: list[tuple[str, str]]) -> list[tuple[str, str, str]]:
    """Split an RCSB FASTA download into its polymer entities.

    The download gives one record per polymer entity, with a header of the
    form::

        >9GCN_1|Chains A, F[auth B]|Protein. Variable Domain of ...|Vicugna pacos

    The fields are the entity identifier, the chains, the description of the
    polymer and the source organism, in that order. The part before the
    first pipe is the entity identifier, so the records can be addressed by
    entity rather than by chain.

    Args:
        records: The records as returned by :func:`mrna.fasta.read_fasta`.

    Returns:
        A list of ``(entity_id, description, sequence)`` triples.

    Raises:
        ValueError: If a header does not carry an entity identifier.
    """
    entities: list[tuple[str, str, str]] = []
    for header, sequence in records:
        fields = header.split("|")
        entity_id = fields[0].strip()
        if "_" not in entity_id:
            msg = f"RCSB record header {header!r} does not carry an entity identifier"
            raise ValueError(msg)
        description = fields[2].strip() if len(fields) > 2 else ""
        entities.append((entity_id, description, sequence))
    return entities


def write_tpl1158_vhh() -> Path:
    """Extract the TPL1158_01_C09 nanobody from the PDB entry.

    This is the only single-domain candidate in the M4 table whose variable
    region can be located in a public database, and it is the payload that
    makes a second cocktail member constructible without engineering.

    Args:
        None.

    Returns:
        The path of the written file.

    Raises:
        ValueError: If the downloaded entry does not hold the expected
            entity.
    """
    records = read_fasta(RAW / "rcsb_9GCN_vhh_tpl1158_01_c09.fasta")
    entities = split_rcsb_records(records)
    vhh = next(
        (seq for entity, description, seq in entities if "Heavy-Chain only" in description),
        None,
    )
    if vhh is None:
        msg = "rcsb_9GCN_vhh_tpl1158_01_c09.fasta holds no VHH entity"
        raise ValueError(msg)
    header = (
        "TPL1158_01_C09 nanobody, anti-Naja kaouthia alpha-cobratoxin | "
        "source: PDB entry 9GCN, entity 1, a co-crystal with alpha-cobratoxin "
        "from Naja kaouthia | released 2025-08-06 | length "
        f"{len(vhh)} aa"
    )
    target = OUT / "tpl1158_01_c09_vhh.fasta"
    write_fasta(target, [(header, vhh)])
    return target


#: The PDB entries whose antibody chains are recorded for the M5 cocktail
#: plan, with the entities that hold the antibody rather than the toxin.
CANDIDATE_ANTIBODY_ENTITIES: tuple[tuple[str, str, tuple[int, ...]], ...] = (
    (
        "rcsb_8DA1_centi_lnx_d09.fasta",
        "Centi-LNX-D09 Fab, anti-long-neurotoxin, from PDB entry 8DA1",
        (1, 2),
    ),
    (
        "rcsb_8V13_centi_snx_b03.fasta",
        "Centi-SNX-B03 Fab, anti-short-neurotoxin, from PDB entry 8V13",
        (1, 2),
    ),
    (
        "rcsb_9CLP_h11_fab.fasta",
        "H11 Fab, anti-ecarin, from PDB entry 9CLP",
        (3, 4),
    ),
    (
        "rcsb_9HUO_a01.fasta",
        "A01 monoclonal antibody, anti-alpha-cobratoxin, from PDB entry 9HUO",
        (2, 3),
    ),
)


def write_candidate_fab_chains() -> Path:
    """Record the antibody chains of the multi-chain candidates.

    These candidates cannot be transcribed as one polypeptide, so their
    chains are kept here so that the cost of converting them is a number
    rather than an adjective.

    Args:
        None.

    Returns:
        The path of the written file.

    Raises:
        ValueError: If a downloaded entry does not hold the expected
            entities.
    """
    records: list[tuple[str, str]] = []
    for filename, label, entity_numbers in CANDIDATE_ANTIBODY_ENTITIES:
        entities = split_rcsb_records(read_fasta(RAW / filename))
        by_number = {
            int(entity.rsplit("_", 1)[1]): (entity, desc, seq)
            for entity, desc, seq in entities
        }
        for number in entity_numbers:
            if number not in by_number:
                msg = f"{filename} holds no entity {number}"
                raise ValueError(msg)
            entity, description, sequence = by_number[number]
            header = (
                f"{label} | entity {entity} | {description} | length {len(sequence)} aa"
            )
            records.append((header, sequence))
    target = OUT / "candidate_fab_chains.fasta"
    write_fasta(target, records)
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

    for path in (
        write_bnt162b2_cds(),
        write_reference_protein(),
        write_reference_gfp(),
        write_e113_payload(),
        write_ecarpholin_s(),
        write_albumin(),
        write_tpl1158_vhh(),
        write_candidate_fab_chains(),
    ):
        print(f"wrote {path.relative_to(ROOT)}")

    print("\nsha256 of raw downloads:")
    for raw in sorted(RAW.glob("*")):
        print(f"  {raw.name}: {sha256(raw)}")


if __name__ == "__main__":
    main()
