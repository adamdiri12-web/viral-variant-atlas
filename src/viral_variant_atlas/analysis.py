from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict

from .types import FastaRecord, Mutation, Placement

_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def allowed_mismatches(
    query_length: int, rate: float, absolute_cap: int | None
) -> int:
    if not 0 <= rate <= 1:
        raise ValueError("max mismatch rate must be between 0 and 1")
    if absolute_cap is not None and absolute_cap < 0:
        raise ValueError("absolute mismatch cap must be non-negative")
    threshold = int(query_length * rate)
    return (
        min(threshold, absolute_cap)
        if absolute_cap is not None
        else threshold
    )


def best_placement(
    reference: str,
    query: FastaRecord,
    rate: float,
    absolute_cap: int | None,
) -> Placement:
    """Find the minimum-Hamming-distance placement across both strands."""
    threshold = allowed_mismatches(len(query.sequence), rate, absolute_cap)
    if len(query.sequence) > len(reference):
        return Placement(
            query.identifier, "+", -1, len(query.sequence), False
        )
    candidates: list[tuple[int, int, str]] = []
    for orientation, sequence in (
        ("+", query.sequence),
        ("-", reverse_complement(query.sequence)),
    ):
        for offset in range(len(reference) - len(sequence) + 1):
            mismatches = sum(
                a != b
                for a, b in zip(
                    reference[offset : offset + len(sequence)], sequence
                )
            )
            candidates.append((mismatches, offset, orientation))
    mismatch_count, offset, orientation = min(
        candidates, key=lambda item: (item[0], item[1], item[2] != "+")
    )
    return Placement(
        query.identifier,
        orientation,
        offset,
        mismatch_count,
        mismatch_count <= threshold,
    )


def call_mutations(
    reference: FastaRecord, query: FastaRecord, placement: Placement
) -> list[Mutation]:
    if not placement.accepted:
        return []
    sequence = (
        query.sequence
        if placement.orientation == "+"
        else reverse_complement(query.sequence)
    )
    window = reference.sequence[
        placement.offset : placement.offset + len(sequence)
    ]
    return [
        Mutation(
            reference.identifier,
            query.identifier,
            placement.orientation,
            placement.offset + index + 1,
            ref_base,
            alt_base,
            placement.offset,
        )
        for index, (ref_base, alt_base) in enumerate(zip(window, sequence))
        if ref_base != alt_base
    ]


def hotspot_rows(mutations: list[Mutation]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, int, str], list[Mutation]] = defaultdict(list)
    for mutation in mutations:
        grouped[
            (
                mutation.reference_id,
                mutation.reference_position,
                mutation.reference_base,
            )
        ].append(mutation)
    rows: list[dict[str, object]] = []
    for (reference_id, position, base), events in grouped.items():
        substitutions = Counter(event.substitution for event in events)
        rows.append(
            {
                "reference_id": reference_id,
                "reference_position": position,
                "reference_base": base,
                "sample_recurrence": len(
                    {event.query_id for event in events}
                ),
                "event_count": len(events),
                "substitution_counts": "; ".join(
                    f"{call}:{count}"
                    for call, count in sorted(substitutions.items())
                ),
            }
        )
    return sorted(
        rows,
        key=lambda row: (
            -int(row["sample_recurrence"]),
            -int(row["event_count"]),
            int(row["reference_position"]),
        ),
    )


def placement_row(placement: Placement) -> dict[str, object]:
    return asdict(placement)
