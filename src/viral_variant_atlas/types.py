from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    sequence: str


@dataclass(frozen=True)
class Placement:
    query_id: str
    orientation: str
    offset: int
    mismatch_count: int
    accepted: bool


@dataclass(frozen=True)
class Mutation:
    reference_id: str
    query_id: str
    orientation: str
    reference_position: int
    reference_base: str
    alternate_base: str
    placement_offset: int

    @property
    def substitution(self) -> str:
        return f"{self.reference_base}{self.reference_position}{self.alternate_base}"

    def row(self) -> dict[str, object]:
        row = asdict(self)
        row["substitution"] = self.substitution
        return row
