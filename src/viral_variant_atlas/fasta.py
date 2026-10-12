from __future__ import annotations

from pathlib import Path

from .types import FastaRecord

IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")


class FastaFormatError(ValueError):
    """Raised when FASTA content is malformed or unsupported."""


def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Parse a FASTA file, retaining record IDs and validating DNA characters."""
    source = Path(path)
    records: list[FastaRecord] = []
    identifier: str | None = None
    fragments: list[str] = []

    def flush() -> None:
        nonlocal fragments
        if identifier is None:
            return
        sequence = "".join(fragments).upper()
        if not sequence:
            raise FastaFormatError(
                f"{source}: record '{identifier}' has no sequence"
            )
        invalid = sorted(set(sequence) - IUPAC_DNA)
        if invalid:
            raise FastaFormatError(
                f"{source}: record '{identifier}' has unsupported bases: {''.join(invalid)}"
            )
        records.append(
            FastaRecord(identifier=identifier, sequence=sequence)
        )
        fragments = []

    for line_number, raw_line in enumerate(
        source.read_text(encoding="utf-8").splitlines(), 1
    ):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header_parts = line[1:].split()
            if not header_parts:
                raise FastaFormatError(
                    f"{source}:{line_number}: FASTA header has no identifier"
                )
            identifier = header_parts[0]
        elif identifier is None:
            raise FastaFormatError(
                f"{source}:{line_number}: sequence appears before a FASTA header"
            )
        else:
            fragments.append("".join(line.split()))
    flush()
    if not records:
        raise FastaFormatError(f"{source}: no FASTA records found")
    return records
