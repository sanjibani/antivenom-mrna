"""The antibody payloads a construct can carry, and what each one costs.

A payload is the antibody itself. The sequences live in ``data/reference/``
because they are evidence rather than choices: every one of them is a
published sequence, copied from the database that publishes it, with the
deposited construct trimmed by script where the deposited construct is not
the molecule.

What this module adds to a sequence is the thing a construct has to know
about it: whether it is one polypeptide or several, and therefore whether
one transcript can carry it or whether it needs engineering before it can be
transcribed at all. That distinction is the whole of the M5 cocktail plan,
because it is the difference between a cocktail that can be built from
published material and one that cannot.
"""

from __future__ import annotations

from dataclasses import dataclass

from mrna.fasta import read_fasta
from mrna.reference import reference_dir

__all__ = [
    "PAYLOAD_SPECS",
    "Payload",
    "PayloadSpec",
    "load_payload",
    "load_payloads",
]


@dataclass(frozen=True)
class Payload:
    """One antibody candidate, with the chains it is made of.

    Attributes:
        name: The clone name as the source publishes it.
        chains: The protein chains of the payload, in order, each a
            ``(label, sequence)`` pair.
        molecule_format: The antibody format, one of
            :data:`mrna.candidates.FORMATS`.
        species: The species the payload targets.
        toxin_family: The toxin family the payload targets.
        toxin_target: The specific toxin or venom component.
        citation: Where the payload and its sequence come from.
    """

    name: str
    chains: tuple[tuple[str, str], ...]
    molecule_format: str
    species: tuple[str, ...]
    toxin_family: str
    toxin_target: str
    citation: str

    @property
    def is_single_polypeptide(self) -> bool:
        """Whether the payload is one chain, and so one open reading frame."""
        return len(self.chains) == 1

    @property
    def length(self) -> int:
        """Total number of residues across every chain of the payload."""
        return sum(len(sequence) for _, sequence in self.chains)

    @property
    def sequence(self) -> str:
        """The payload as one string, which is only meaningful when single-chain."""
        return "".join(sequence for _, sequence in self.chains)

    def chain_lengths(self) -> dict[str, int]:
        """Return the length of each chain, keyed by its label.

        Returns:
            Mapping of chain label to residue count.
        """
        return {label: len(sequence) for label, sequence in self.chains}


@dataclass(frozen=True)
class PayloadSpec:
    """Where a payload comes from and what it is, before it is loaded.

    Attributes:
        name: The clone name as the source publishes it.
        reference: The file under ``data/reference/`` holding the sequence.
        chain_labels: The label for each chain, in order.
        selectors: The substring that identifies each record in the
            reference file. Empty means the file holds exactly one record
            and it is the payload.
        molecule_format: The antibody format.
        species: The species the payload targets.
        toxin_family: The toxin family the payload targets.
        toxin_target: The specific toxin or venom component.
        citation: Where the payload and its sequence come from.
    """

    name: str
    reference: str
    chain_labels: tuple[str, ...]
    selectors: tuple[str, ...]
    molecule_format: str
    species: tuple[str, ...]
    toxin_family: str
    toxin_target: str
    citation: str


#: The payloads the repository knows about.
#:
#: The two single-chain entries are the ones that can be transcribed from a
#: published sequence as they stand. The four multi-chain entries cannot:
#: they are Fab fragments or an antibody with two chains, so each one needs
#: to be turned into a single polypeptide before one transcript can carry
#: it, and that conversion is an engineering step that changes the molecule
#: the published data belongs to.
PAYLOAD_SPECS: tuple[PayloadSpec, ...] = (
    PayloadSpec(
        name="E113",
        reference="e113_scfv.fasta",
        chain_labels=("VH-linker-VL",),
        selectors=(),
        molecule_format="scFv",
        species=("Echis carinatus",),
        toxin_family="PLA2",
        toxin_target="ecarpholin S, the Ser49 myotoxic phospholipase A2",
        citation=(
            "Kumar A, Madni ZK, Chaturvedi S, Salunke DM. 2024. Mol Immunol "
            "165:55-67. doi:10.1016/j.molimm.2023.12.006. PDB entry 8IA6, "
            "entity 1, with the crystallisation vector residue and the "
            "hexahistidine tag removed by scripts/prepare_reference_data.py."
        ),
    ),
    PayloadSpec(
        name="TPL1158_01_C09",
        reference="tpl1158_01_c09_vhh.fasta",
        chain_labels=("VHH",),
        selectors=(),
        molecule_format="VHH",
        species=("Naja kaouthia",),
        toxin_family="3FTx-L",
        toxin_target="alpha-cobratoxin",
        citation=(
            "PDB entry 9GCN, entity 1: a co-crystal of the camelid "
            "single-domain antibody TPL1158_01_C09 with alpha-cobratoxin from "
            "Naja kaouthia. Released 2025-08-06."
        ),
    ),
    PayloadSpec(
        name="Centi-LNX-D09",
        reference="candidate_fab_chains.fasta",
        chain_labels=("light", "heavy"),
        selectors=("8DA1_1", "8DA1_2"),
        molecule_format="Fab",
        species=("Naja naja", "Bungarus caeruleus"),
        toxin_family="3FTx-L",
        toxin_target="long neurotoxin",
        citation=(
            "PDB entry 8DA1, a complex of the Centi-LNX-D09 Fab with "
            "alpha-bungarotoxin from Bungarus multicinctus. Glanville J, Kwong "
            "K, et al. 2025. Cell. Patent US20240294618A1."
        ),
    ),
    PayloadSpec(
        name="Centi-SNX-B03",
        reference="candidate_fab_chains.fasta",
        chain_labels=("heavy", "light"),
        selectors=("8V13_1", "8V13_2"),
        molecule_format="Fab",
        species=("Naja naja",),
        toxin_family="3FTx-S",
        toxin_target="short neurotoxin",
        citation=(
            "PDB entry 8V13, a complex of the Centi-SNX-B03 Fab with a short "
            "neurotoxin from Dendroaspis polylepis polylepis. Glanville J, "
            "Kwong K, et al. 2025. Cell. Patent US20240294618A1."
        ),
    ),
    PayloadSpec(
        name="H11",
        reference="candidate_fab_chains.fasta",
        chain_labels=("heavy", "light"),
        selectors=("9CLP_3", "9CLP_4"),
        molecule_format="Fab",
        species=("Echis carinatus",),
        toxin_family="SVMP",
        toxin_target="ecarin, a zinc metalloproteinase-disintegrin",
        citation=(
            "PDB entry 9CLP: ecarin from the venom of the Kenyan saw-scaled "
            "viper in complex with the Fab of the neutralizing antibody H11. "
            "Released 2024-09-11."
        ),
    ),
    PayloadSpec(
        name="A01",
        reference="candidate_fab_chains.fasta",
        chain_labels=("light", "heavy"),
        selectors=("9HUO_2", "9HUO_3"),
        molecule_format="mAb-IgG",
        species=("Naja kaouthia",),
        toxin_family="3FTx-L",
        toxin_target="alpha-cobratoxin",
        citation=(
            "PDB entries 9HUO, 9FYT and 9HXO: the A01 monoclonal antibody "
            "bound to alpha-cobratoxin at pH 5.5, 4.5 and 6.0. Released "
            "2025-11-19."
        ),
    ),
)


def load_payload(spec: PayloadSpec) -> Payload:
    """Load one payload from the committed reference data.

    Args:
        spec: The specification of the payload.

    Returns:
        The loaded :class:`Payload`.

    Raises:
        ValueError: If the reference file is missing a record, or if the
            number of chains does not match the specification.
    """
    path = reference_dir() / spec.reference
    records = read_fasta(path)
    if spec.selectors:
        chains: list[tuple[str, str]] = []
        for selector in spec.selectors:
            match = next((record for record in records if selector in record[0]), None)
            if match is None:
                msg = f"{path.name} holds no record matching {selector!r}"
                raise ValueError(msg)
            chains.append((selector, match[1]))
    else:
        if not records:
            msg = f"{path.name} holds no records"
            raise ValueError(msg)
        chains = [(spec.name, records[0][1])]
    if len(chains) != len(spec.chain_labels):
        msg = (
            f"{spec.name} has {len(chains)} chains but the specification "
            f"names {len(spec.chain_labels)}"
        )
        raise ValueError(msg)
    return Payload(
        name=spec.name,
        chains=tuple(
            (label, sequence)
            for label, (_, sequence) in zip(spec.chain_labels, chains, strict=True)
        ),
        molecule_format=spec.molecule_format,
        species=spec.species,
        toxin_family=spec.toxin_family,
        toxin_target=spec.toxin_target,
        citation=spec.citation,
    )


def load_payloads() -> tuple[Payload, ...]:
    """Load every payload the repository knows about.

    Returns:
        The payloads, in specification order.
    """
    return tuple(load_payload(spec) for spec in PAYLOAD_SPECS)
