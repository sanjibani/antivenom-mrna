"""Primitive DNA sequence operations for mRNA construct design.

This module is the foundation layer of the antivenom-mrna design engine. It holds
no design policy. Every function is a measurement or a mechanical transformation
that a reviewer can verify by hand on paper. Design policy (GC windows, CAI
targets, motif blacklists, UTR choices) lives in :mod:`mrna.params` and
:mod:`mrna.design`, and every value there carries a citation.

Biological background
---------------------
An mRNA drug is one linear chain, read 5' to 3' by the ribosome::

    5'-cap | 5' UTR | CDS | 3' UTR | poly(A)- 3'

The ribosome has to find the single ``AUG`` start codon inside the 5' UTR before
it can translate anything, so the UTRs are engineered to be functionally quiet
while the coding region is engineered to be translated efficiently. That split
is why this module separates *measurement* (GC content, translation) from
*choice* (which UTR, which codon, which poly(A) length).

Sequence conventions
--------------------
Sequences are handled as DNA strings of ``A``, ``C``, ``G``, ``T`` (the ``U`` of
RNA is written as ``T`` throughout, matching GenBank and patent practice). Input
is case-insensitive and internal whitespace is ignored. Nucleoside
modifications that matter for manufacturing (N1-methylpseudouridine, Cap 1,
poly(A) tail length) are *not* visible in a DNA string and are tracked
separately in the construct specification.
"""

from __future__ import annotations

from collections.abc import Mapping

__all__ = [
    "AA_TO_CODONS",
    "CODON_TABLE",
    "STOP_CODONS",
    "back_translate",
    "gc_content",
    "translate",
]

#: The standard genetic code, NCBI translation table 1.
#:
#: Source: NCBI Genetic Codes, NCBI Bookshelf NBK26827 (The Genetic Codes,
#: 2008). Table 1 is used by every human cell and by every human vaccine
#: construct, so it is the only table this engine supports. Mitochondrial
#: table 11 and the bacterial/prokaryotic tables 4, 11 and 16 are deliberately
#: not implemented: a venom-coding construct destined for a mammalian host
#: must be read with table 1.
CODON_TABLE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

#: The three codons that terminate translation in the standard code.
STOP_CODONS: frozenset[str] = frozenset(
    codon for codon, aa in CODON_TABLE.items() if aa == "*"
)

#: Sense codons grouped by amino acid, each group sorted alphabetically.
#:
#: Sorting makes the group order deterministic, which matters because
#: :func:`back_translate` breaks ties between equally used codons by picking the
#: alphabetically first one. Reproducibility is a design requirement here, not
#: a nicety.
AA_TO_CODONS: dict[str, tuple[str, ...]] = {}
for _codon, _aa in sorted(CODON_TABLE.items()):
    if _aa != "*":
        AA_TO_CODONS.setdefault(_aa, ())
        AA_TO_CODONS[_aa] = AA_TO_CODONS[_aa] + (_codon,)
del _codon, _aa


def _normalise(seq: str, *, caller: str) -> str:
    """Uppercase a sequence and drop whitespace.

    Args:
        seq: Raw sequence, possibly with embedded spaces, tabs or newlines.
        caller: Name of the calling function, used in error messages.

    Returns:
        The normalised, uppercase sequence.

    Raises:
        ValueError: If the sequence is empty after normalisation.
    """
    cleaned = "".join(seq.split()).upper()
    if not cleaned:
        msg = f"{caller}() received an empty sequence"
        raise ValueError(msg)
    return cleaned


def translate(
    seq: str,
    *,
    allow_ambiguous: bool = False,
    stop_at_stop: bool = False,
) -> str:
    """Translate a coding sequence into its protein, using the standard code.

    Translation is the inverse of codon assignment. For an mRNA drug it is how
    the engine proves that a nucleotide sequence actually encodes the intended
    antigen: a design whose CDS does not translate back to the target protein
    is broken, and this is the check that catches it before synthesis.

    Args:
        seq: Coding sequence in DNA alphabet (``A``, ``C``, ``G``, ``T``).
            Case-insensitive; whitespace is ignored.
        allow_ambiguous: If True, any character outside ``ACGT`` becomes ``X``
            in the protein and the codon contributes a single ``X``. If False
            (the default) such characters raise an error, because a silently
            mistranslated construct is worse than a failed one.
        stop_at_stop: If True, the returned protein is truncated at the first
            stop codon, keeping the ``*`` marker. If False, translation
            continues through every stop codon and each one is rendered as
            ``*``.

    Returns:
        The protein sequence, one letter per codon.

    Raises:
        ValueError: If the sequence is empty, contains characters outside
            ``ACGT`` while ``allow_ambiguous`` is False, or its length is not
            a multiple of three (a truncated codon cannot be translated).

    Examples:
        >>> translate("ATGGGCCATTGTGATGGCTGG")
        'MGHCDGW'
    """
    cleaned = _normalise(seq, caller="translate")

    invalid = sorted(set(cleaned) - set("ACGT"))
    if invalid and not allow_ambiguous:
        msg = (
            f"translate() found non-ACGT characters {''.join(invalid)}. "
            "Pass allow_ambiguous=True to map them to 'X'."
        )
        raise ValueError(msg)

    remainder = len(cleaned) % 3
    if remainder:
        msg = (
            f"translate() needs a sequence whose length is a multiple of 3, "
            f"got {len(cleaned)} nt (trailing {remainder} nt cannot form a codon)"
        )
        raise ValueError(msg)

    protein_parts: list[str] = []
    for index in range(0, len(cleaned), 3):
        codon = cleaned[index : index + 3]
        if allow_ambiguous and not set(codon) <= set("ACGT"):
            protein_parts.append("X")
        else:
            protein_parts.append(CODON_TABLE[codon])

    protein = "".join(protein_parts)

    if stop_at_stop and "*" in protein:
        protein = protein[: protein.index("*") + 1]
    return protein


def back_translate(protein: str, usage_table: Mapping[str, float]) -> str:
    """Back-translate a protein into DNA using a host codon usage table.

    For each amino acid the codon with the highest synonymous usage in the
    target host is chosen. Because relative adaptiveness
    ``w_i = f_i / f_max`` (Sharp and Li, 1987) is a monotone function of the
    usage fraction ``f_i``, maximising ``f_i`` also maximises ``w_i``, which
    is the quantity the Codon Adaptation Index actually measures. Ties are
    broken alphabetically so that the same inputs always produce byte-identical
    output; a design that changes when you re-run it is not auditable.

    This function is deliberately naive. It ignores GC balance, long homopolymer
    runs and cryptic splice or polyadenylation motifs, all of which matter for
    a real construct. Those constraints belong to the codon optimiser added in
    milestone M1, which will call this function as one of its candidate
    generators.

    Args:
        protein: Protein sequence, one-letter codes. Must not contain ``*``;
            stop codons are the job of the CDS assembler.
        usage_table: Mapping of codon to synonymous codon usage for the target
            host, as relative fractions (the sum does not need to be 1).
            A human table can be taken from the Kazusa Codon Usage Database
            (Nakamura et al., 2000).

    Returns:
        The back-translated coding sequence.

    Raises:
        ValueError: If the protein is empty, contains ``*`` or contains
            characters that are not standard amino acid codes, or if the usage
            table is missing any codon needed by the protein.

    Examples:
        >>> table = {"ATG": 1.0, "TGG": 1.0, "CTG": 1.0, "CAG": 1.0}
        >>> back_translate("MWL", table)
        'ATGTGGCTGCAG'
    """
    cleaned = "".join(protein.split()).upper()
    if not cleaned:
        msg = "back_translate() received an empty protein sequence"
        raise ValueError(msg)

    if "*" in cleaned:
        msg = (
            "back_translate() does not emit stop codons; '*' was found in the "
            "protein. Handle termination in the CDS assembler instead."
        )
        raise ValueError(msg)

    unknown = sorted(set(cleaned) - set(AA_TO_CODONS))
    if unknown:
        msg = f"back_translate() found non-standard amino acids {''.join(unknown)}"
        raise ValueError(msg)

    needed = {codon for aa in set(cleaned) for codon in AA_TO_CODONS[aa]}
    missing = sorted(needed - set(usage_table))
    if missing:
        msg = (
            f"back_translate() needs {len(missing)} codons that the usage table "
            f"does not contain: {', '.join(missing)}"
        )
        raise ValueError(msg)

    codons_out: list[str] = []
    for aa in cleaned:
        ranked = sorted(AA_TO_CODONS[aa], key=lambda c: (-float(usage_table[c]), c))
        codons_out.append(ranked[0])
    return "".join(codons_out)


def gc_content(seq: str, window: int | None = None, *, step: int = 1) -> float | list[float]:
    """Measure GC content, overall or across a sliding window.

    GC content is the fraction of guanine and cytosine bases in a sequence. It
    matters twice over for mRNA drugs. First, a high-GC first coding window
    promotes translation initiation, because the GC-rich region around the
    start codon eases the scanning of the small ribosomal subunit
    (Kudla et al., 1999). Second, extreme GC in either direction destabilises
    the transcript: GC-rich stretches stack into secondary structure that
    slows the ribosome, while very AT-rich stretches make the message look
    more like a degradation intermediate.

    The overall number tells you where a construct sits. The sliding window
    tells you where it is *uneven*, which is the thing the overall number hides
    and the thing that actually causes a failed batch.

    Ambiguous characters (for example ``N``) are counted in the denominator
    but contribute no G or C, so the result is a lower bound on true GC content
    for sequences containing ambiguity codes.

    Args:
        seq: Sequence in DNA alphabet. Case-insensitive; whitespace is ignored.
        window: If None, return the overall GC fraction as a float. If an
            integer, return the GC fraction in every sliding window of that
            size, as a list.
        step: Stride between consecutive windows when ``window`` is given.
            The default of 1 advances one base at a time.

    Returns:
        A float when ``window`` is None, otherwise a list of floats with
        ``(len(seq) - window) // step + 1`` entries.

    Raises:
        ValueError: If the sequence is empty, if ``window`` is not positive,
            if ``window`` is longer than the sequence, or if ``step`` is not
            positive.

    Examples:
        >>> gc_content("ATGC")
        0.5
        >>> gc_content("ATGCGC", window=3)
        [0.3333333333333333, 0.6666666666666666, 1.0, 1.0]
    """
    cleaned = _normalise(seq, caller="gc_content")

    if window is None:
        gc = sum(cleaned.count(base) for base in "GC")
        return gc / len(cleaned)

    if isinstance(window, bool) or not isinstance(window, int):
        msg = f"gc_content() window must be an int or None, got {type(window).__name__}"
        raise TypeError(msg)
    if window < 1:
        msg = f"gc_content() window must be positive, got {window}"
        raise ValueError(msg)
    if window > len(cleaned):
        msg = (
            f"gc_content() window of {window} nt is longer than the sequence "
            f"({len(cleaned)} nt)"
        )
        raise ValueError(msg)
    if isinstance(step, bool) or not isinstance(step, int):
        msg = f"gc_content() step must be an int, got {type(step).__name__}"
        raise TypeError(msg)
    if step < 1:
        msg = f"gc_content() step must be positive, got {step}"
        raise ValueError(msg)

    prefix_sums = [0]
    for base in cleaned:
        prefix_sums.append(prefix_sums[-1] + (1 if base in "GC" else 0))

    return [
        (prefix_sums[start + window] - prefix_sums[start]) / window
        for start in range(0, len(cleaned) - window + 1, step)
    ]
