"""Structural analysis of mRNA constructs.

Everything in this module is a *measurement* of a nucleotide sequence: where
the open reading frame is, whether the start codon sits in a good initiation
context, how GC-rich the coding region is, how well the coding region matches
the host's codon preferences, and which amino acids differ from a reference
protein.

The measurements exist to answer one question at every step of the project:
"is this construct what its author says it is?" The reference case is the
licensed COVID-19 vaccines, whose published sequences are the only mRNA drugs
whose designs are public. Reverse-engineering them is how the engine earns
the right to be trusted on an antivenom construct nobody has published.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from mrna.params import GC_MAX, GC_MIN, GC_WINDOW_NT, KOZAK_MINUS3_ALLOWED, KOZAK_PLUS4
from mrna.sequence import CODON_TABLE, STOP_CODONS, gc_content, translate

__all__ = [
    "ConstructReport",
    "KozakContext",
    "LayoutCheck",
    "LayoutElement",
    "Orf",
    "Variant",
    "align_proteins",
    "analyze_cds",
    "codon_adaptation_index",
    "codon_counts",
    "find_orf",
    "kozak_context",
    "protein_variants",
    "relative_adaptiveness",
    "split_construct",
    "verify_layout",
]


@dataclass(frozen=True)
class Orf:
    """An open reading frame located inside a nucleotide sequence.

    Attributes:
        start: Zero-based index of the first base of the start codon.
        end: Zero-based index one past the last base of the stop codon.
        cds: The coding sequence, start codon through stop codon inclusive.
        protein: The translation of :attr:`cds`, with ``*`` marking the stop
            codon.
    """

    start: int
    end: int
    cds: str
    protein: str


@dataclass(frozen=True)
class KozakContext:
    """The nucleotide context immediately around an initiation codon.

    Attributes:
        minus3: The base three positions before the A of AUG.
        minus2: The base two positions before the A of AUG.
        minus1: The base immediately before the A of AUG.
        plus4: The base immediately after the G of AUG.
        context: The nine bases spanning positions -3 to +4.
    """

    minus3: str
    minus2: str
    minus1: str
    plus4: str
    context: str

    @property
    def matches_consensus(self) -> bool:
        """Whether both key positions follow the Kozak consensus.

        Kozak 1986 found that position -3 (a purine) and position +4 (a G) are
        the two positions that matter most for initiation efficiency.
        """
        return self.minus3 in KOZAK_MINUS3_ALLOWED and self.plus4 == KOZAK_PLUS4

    @property
    def deviations(self) -> tuple[str, ...]:
        """Human-readable names of the consensus positions that deviate."""
        problems: list[str] = []
        if self.minus3 not in KOZAK_MINUS3_ALLOWED:
            problems.append("-3 is not a purine (A/G)")
        if self.plus4 != KOZAK_PLUS4:
            problems.append("+4 is not G")
        return tuple(problems)


@dataclass(frozen=True)
class Variant:
    """A single amino-acid difference between two aligned proteins.

    Attributes:
        position: One-based position in the reference protein.
        reference_aa: The amino acid in the reference protein.
        variant_aa: The amino acid in the variant protein.
        variant_position: One-based position in the variant protein, or None
            when the variant protein carries an unaligned terminal extension.
    """

    position: int
    reference_aa: str
    variant_aa: str
    variant_position: int | None = None

    def describe(self) -> str:
        """Return a human-readable description such as ``S21P``."""
        return f"{self.reference_aa}{self.position}{self.variant_aa}"


@dataclass(frozen=True)
class ConstructReport:
    """Everything the engine measured about one coding sequence.

    Attributes:
        name: Name of the construct.
        source: Where the sequence came from, in one line.
        cds_length: Length of the coding sequence in nucleotides, including
            the stop codon.
        protein_length: Length of the encoded protein in amino acids,
            excluding the stop codon.
        gc_total: GC fraction over the whole coding sequence.
        gc_first_window: GC fraction over the first :data:`GC_WINDOW_NT`
            nucleotides, the region that controls translation initiation.
        gc_min_window: Lowest GC fraction found in any sliding window.
        gc_max_window: Highest GC fraction found in any sliding window.
        cai: Codon Adaptation Index against the supplied host usage table.
        kozak: Initiation context of the start codon, or None when the
            5' UTR was not available and the context cannot be read.
        variants: Amino-acid differences against the reference protein.
        protein: The encoded protein, without the trailing stop marker.
        alignment_offset: Number of extra N-terminal residues the variant
            carries relative to the reference protein, zero when the two are
            the same length.
    """

    name: str
    source: str
    cds_length: int
    protein_length: int
    gc_total: float
    gc_first_window: float
    gc_min_window: float
    gc_max_window: float
    cai: float
    kozak: KozakContext | None
    variants: tuple[Variant, ...]
    protein: str
    alignment_offset: int = 0

    @property
    def gc_is_in_band(self) -> bool:
        """Whether the overall GC content sits inside the accepted band."""
        return GC_MIN <= self.gc_total <= GC_MAX


def find_orf(seq: str, *, start_codon: str = "ATG") -> Orf:
    """Locate the first complete open reading frame in a nucleotide sequence.

    The search starts at the first occurrence of ``start_codon`` and ends at
    the first in-frame stop codon. A construct whose coding region does not
    close on a stop codon is a broken construct, so this function refuses to
    return a partial frame.

    Args:
        seq: Nucleotide sequence, whitespace tolerated.
        start_codon: The codon that opens the reading frame.

    Returns:
        The located :class:`Orf`.

    Raises:
        ValueError: If no start codon is found, if no in-frame stop codon
            follows it, or if the start codon is not a complete triplet in the
            standard code.
    """
    cleaned = "".join(seq.split()).upper()
    if start_codon not in CODON_TABLE:
        msg = f"find_orf() got {start_codon!r}, which is not a codon in the standard code"
        raise ValueError(msg)

    start = cleaned.find(start_codon)
    if start < 0:
        msg = f"find_orf() found no {start_codon} start codon in the sequence"
        raise ValueError(msg)

    protein_parts: list[str] = []
    for index in range(start, len(cleaned) - 2, 3):
        codon = cleaned[index : index + 3]
        amino_acid = CODON_TABLE[codon]
        protein_parts.append(amino_acid)
        if codon in STOP_CODONS:
            end = index + 3
            return Orf(start=start, end=end, cds=cleaned[start:end], protein="".join(protein_parts))

    msg = f"find_orf() found no in-frame stop codon after position {start}"
    raise ValueError(msg)


def kozak_context(seq: str, atg_index: int) -> KozakContext:
    """Read the initiation context around an ATG.

    Args:
        seq: Nucleotide sequence containing the start codon.
        atg_index: Zero-based index of the A of the ATG.

    Returns:
        The :class:`KozakContext` for that start codon.

    Raises:
        ValueError: If there are fewer than three bases before the start codon
            or fewer than one base after it.
    """
    cleaned = "".join(seq.split()).upper()
    if atg_index < 3:
        msg = f"kozak_context() needs three bases before position {atg_index}"
        raise ValueError(msg)
    if atg_index + 3 >= len(cleaned):
        msg = f"kozak_context() needs one base after the ATG at position {atg_index}"
        raise ValueError(msg)

    return KozakContext(
        minus3=cleaned[atg_index - 3],
        minus2=cleaned[atg_index - 2],
        minus1=cleaned[atg_index - 1],
        plus4=cleaned[atg_index + 3],
        context=cleaned[atg_index - 3 : atg_index + 6],
    )


def protein_variants(reference: str, variant: str) -> tuple[Variant, ...]:
    """Compare two proteins of the same length, position by position.

    Args:
        reference: The reference protein.
        variant: The protein to compare against the reference.

    Returns:
        The differences, ordered by position.

    Raises:
        ValueError: If the two proteins have different lengths.
    """
    ref = "".join(reference.split()).upper()
    var = "".join(variant.split()).upper()
    if len(ref) != len(var):
        msg = f"protein_variants() needs equal lengths, got {len(ref)} and {len(var)}"
        raise ValueError(msg)

    return tuple(
        Variant(position=index + 1, reference_aa=ref[index], variant_aa=var[index])
        for index in range(len(ref))
        if ref[index] != var[index]
    )


def relative_adaptiveness(usage_table: Mapping[str, float]) -> dict[str, float]:
    """Turn a codon usage table into relative adaptiveness weights.

    The relative adaptiveness of a codon is its usage frequency divided by the
    highest usage frequency among the codons for the same amino acid, so the
    best codon of each amino acid has weight 1.0 (Sharp and Li, 1987).

    Args:
        usage_table: Mapping of codon to synonymous usage frequency.

    Returns:
        Mapping of codon to relative adaptiveness, for every codon present in
        the table that encodes an amino acid.
    """
    best_per_amino_acid: dict[str, float] = {}
    for codon, frequency in usage_table.items():
        amino_acid = CODON_TABLE.get(codon.upper())
        if amino_acid is None or amino_acid == "*":
            continue
        current = best_per_amino_acid.get(amino_acid, 0.0)
        best_per_amino_acid[amino_acid] = max(current, float(frequency))

    weights: dict[str, float] = {}
    for codon, frequency in usage_table.items():
        amino_acid = CODON_TABLE.get(codon.upper())
        if amino_acid is None or amino_acid == "*":
            continue
        best = best_per_amino_acid[amino_acid]
        weights[codon.upper()] = float(frequency) / best if best else 0.0
    return weights


def codon_adaptation_index(cds: str, usage_table: Mapping[str, float]) -> float:
    """Compute the Codon Adaptation Index of a coding sequence.

    The CAI is the geometric mean of the relative adaptiveness of the codons
    in the sequence, ignoring the initiator ATG and any stop codon (Sharp and
    Li, 1987). A value near 1.0 means the sequence uses the codons the host
    prefers; a value near 0.2 means it actively fights the host's translation
    machinery.

    This is an independent implementation, written so that the engine does not
    have to trust a third-party package for a number it publishes. The
    validation notebook cross-checks it against the ``CAI`` package.

    Args:
        cds: The coding sequence.
        usage_table: Mapping of codon to synonymous usage frequency for the
            target host.

    Returns:
        The CAI, between 0 and 1.

    Raises:
        ValueError: If the sequence length is not a multiple of three, or if
            the usage table does not cover every codon used.
    """
    cleaned = "".join(cds.split()).upper()
    if len(cleaned) % 3:
        msg = f"codon_adaptation_index() needs a length divisible by 3, got {len(cleaned)}"
        raise ValueError(msg)

    weights = relative_adaptiveness(usage_table)
    codons = [cleaned[i : i + 3] for i in range(0, len(cleaned), 3)]
    sense_codons = [codon for codon in codons if codon not in STOP_CODONS]
    # The initiator ATG is excluded from the CAI, because its weight is
    # always 1.0 and including it would bias every sequence the same way.
    # A sequence that is nothing but methionine keeps its initiator, since
    # otherwise there would be nothing left to score.
    if len(sense_codons) > 1 and sense_codons[0] == "ATG":
        sense_codons = sense_codons[1:]
    if not sense_codons:
        msg = "codon_adaptation_index() found no sense codons to score"
        raise ValueError(msg)

    missing = sorted({codon for codon in sense_codons if codon not in weights})
    if missing:
        msg = f"codon_adaptation_index() is missing usage for {', '.join(missing)}"
        raise ValueError(msg)

    log_sum = sum(math.log(weights[codon]) for codon in sense_codons)
    return math.exp(log_sum / len(sense_codons))


def codon_counts(cds: str) -> dict[str, int]:
    """Count how often each codon is used in a coding sequence.

    Args:
        cds: The coding sequence.

    Returns:
        Mapping of codon to number of occurrences, sorted by codon.

    Raises:
        ValueError: If the sequence length is not a multiple of three.
    """
    cleaned = "".join(cds.split()).upper()
    if len(cleaned) % 3:
        msg = f"codon_counts() needs a length divisible by 3, got {len(cleaned)}"
        raise ValueError(msg)

    counts: dict[str, int] = {}
    for index in range(0, len(cleaned), 3):
        codon = cleaned[index : index + 3]
        counts[codon] = counts.get(codon, 0) + 1
    return dict(sorted(counts.items()))


def split_construct(seq: str) -> tuple[str, str, str, str]:
    """Split a full mRNA sequence into its four architectural parts.

    The parts are the 5' UTR, the coding sequence (including its stop codon),
    the 3' UTR, and the poly(A) tail. The poly(A) tail is detected as the
    trailing run of adenines; the run must be at least ten bases long, because
    short adenine runs occur by chance inside UTRs.

    Args:
        seq: The full mRNA sequence.

    Returns:
        A tuple ``(utr5, cds, utr3, poly_a)``.

    Raises:
        ValueError: If no complete open reading frame can be found.
    """
    cleaned = "".join(seq.split()).upper()
    orf = find_orf(cleaned)
    utr5 = cleaned[: orf.start]
    utr3 = cleaned[orf.end :]

    poly_a = ""
    stripped = utr3.rstrip("A")
    trailing = len(utr3) - len(stripped)
    if trailing >= 10:
        poly_a = utr3[len(stripped) :]
        utr3 = stripped
    return utr5, orf.cds, utr3, poly_a


def align_proteins(
    reference: str,
    variant: str,
    *,
    max_offset: int = 20,
) -> tuple[int, tuple[Variant, ...]]:
    """Align a variant protein against a reference protein.

    Vaccine constructs are not always the same length as the wild-type
    protein they are derived from. Moderna's mRNA-1273, for example, encodes a
    spike that starts a few residues upstream of the canonical sequence. A
    naive position-by-position diff of two proteins of different lengths
    reports every residue after the difference as changed, which is noise.

    This function searches offsets in ``[-max_offset, max_offset]``, where a
    positive offset means the variant carries extra residues at its N
    terminus, and returns the offset with the fewest mismatches. Ties are
    broken towards the smallest absolute offset, that is, towards the
    interpretation that changes the least.

    Args:
        reference: The reference protein.
        variant: The protein to align.
        max_offset: Largest terminal offset to consider.

    Returns:
        A tuple ``(offset, variants)``. The offset is the number of extra
        N-terminal residues in the variant.

    Raises:
        ValueError: If either protein is empty, or if no offset produces an
            overlap at all.
    """
    ref = "".join(reference.split()).upper()
    var = "".join(variant.split()).upper()
    if not ref or not var:
        msg = "align_proteins() needs two non-empty proteins"
        raise ValueError(msg)

    best_offset: int | None = None
    best_mismatches = -1
    best_variants: tuple[Variant, ...] = ()

    for offset in range(-max_offset, max_offset + 1):
        overlap = min(len(ref), len(var) - offset)
        if offset < 0:
            overlap = min(len(ref) + offset, len(var))
        if overlap <= 0:
            continue

        mismatches: list[Variant] = []
        for index in range(overlap):
            reference_aa = ref[index]
            variant_aa = var[index + offset]
            if reference_aa != variant_aa:
                mismatches.append(
                    Variant(
                        position=index + 1,
                        reference_aa=reference_aa,
                        variant_aa=variant_aa,
                        variant_position=index + offset + 1,
                    ),
                )

        unaligned = abs(len(ref) - len(var)) + abs(offset)
        score = len(mismatches) + unaligned
        if best_offset is None or score < best_mismatches:
            best_offset = offset
            best_mismatches = score
            best_variants = tuple(mismatches)

    if best_offset is None:
        msg = "align_proteins() found no offset with any overlap"
        raise ValueError(msg)
    return best_offset, best_variants


@dataclass(frozen=True)
class LayoutElement:
    """One published element of a construct's architecture.

    Attributes:
        name: Name of the element, for example ``5'-UTR``.
        start: First nucleotide of the element, 1-based and inclusive.
        end: Last nucleotide of the element, 1-based and inclusive.
        description: What the element is, in the words of the source.
    """

    name: str
    start: int
    end: int
    description: str = ""

    @property
    def length(self) -> int:
        """Length of the element in nucleotides."""
        return self.end - self.start + 1


@dataclass(frozen=True)
class LayoutCheck:
    """The result of checking one published element against a sequence.

    Attributes:
        element: The published element.
        observed_length: Length of the slice of the sequence that the element
            claims, as read from the sequence itself.
        observed_sequence: The nucleotides of that slice.
    """

    element: LayoutElement
    observed_length: int
    observed_sequence: str

    @property
    def matches(self) -> bool:
        """Whether the sequence is long enough to contain the element."""
        return self.observed_length == self.element.length


def verify_layout(seq: str, elements: Sequence[LayoutElement]) -> tuple[LayoutCheck, ...]:
    """Check a full sequence against a published table of elements.

    This is how the engine validates itself against a regulatory document. A
    published architecture gives a name and a coordinate range for every part
    of a construct. Slicing the sequence at those coordinates and measuring
    what comes back turns a drawing into a test.

    Args:
        seq: The full construct sequence.
        elements: The published elements, in order.

    Returns:
        One check per element.

    Raises:
        ValueError: If an element has an end before its start.
    """
    cleaned = "".join(seq.split()).upper()
    checks: list[LayoutCheck] = []
    for element in elements:
        if element.end < element.start:
            msg = f"layout element {element.name!r} ends before it starts"
            raise ValueError(msg)
        observed = cleaned[element.start - 1 : element.end]
        checks.append(
            LayoutCheck(
                element=element,
                observed_length=len(observed),
                observed_sequence=observed,
            ),
        )
    return tuple(checks)


def analyze_cds(
    name: str,
    source: str,
    cds: str,
    usage_table: Mapping[str, float],
    *,
    reference_protein: str | None = None,
    utr5: str | None = None,
) -> ConstructReport:
    """Measure a coding sequence and return everything worth reporting.

    Args:
        name: Name of the construct.
        source: One-line description of where the sequence came from.
        cds: The coding sequence, starting at the start codon.
        usage_table: Codon usage table for the target host.
        reference_protein: Optional reference protein to diff against. The
            diff is how the engine confirms that a construct encodes what it
            claims, for example that the only changes to a viral spike are the
            intended stabilising mutations.
        utr5: Optional 5' UTR. When present, the Kozak initiation context is
            read from the real junction between UTR and coding region. When
            absent, the context is reported as unknown rather than guessed.

    Returns:
        A :class:`ConstructReport` with every measurement filled in.

    Raises:
        ValueError: If the sequence is not a complete coding sequence or the
            reference protein has a different length.
    """
    cleaned = "".join(cds.split()).upper()
    if len(cleaned) % 3:
        msg = f"analyze_cds() needs a length divisible by 3, got {len(cleaned)}"
        raise ValueError(msg)

    protein = translate(cleaned)
    protein_without_stop = protein[:-1] if protein.endswith("*") else protein

    windows = gc_content(cleaned, window=GC_WINDOW_NT)
    total_gc = gc_content(cleaned)
    first_window = windows[0] if windows else total_gc

    variants: tuple[Variant, ...] = ()
    alignment_offset = 0
    if reference_protein is not None:
        alignment_offset, variants = align_proteins(reference_protein, protein_without_stop)

    kozak: KozakContext | None = None
    if utr5 is not None:
        junction = f"{utr5[-3:].upper()}{cleaned[:6]}"
        kozak = kozak_context(junction, 3)

    return ConstructReport(
        name=name,
        source=source,
        cds_length=len(cleaned),
        protein_length=len(protein_without_stop),
        gc_total=float(total_gc),
        gc_first_window=float(first_window),
        gc_min_window=min(windows) if windows else float(total_gc),
        gc_max_window=max(windows) if windows else float(total_gc),
        cai=codon_adaptation_index(cleaned, usage_table),
        kozak=kozak,
        variants=variants,
        protein=protein_without_stop,
        alignment_offset=alignment_offset,
    )
