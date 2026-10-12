from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from .analysis import (
    best_placement,
    call_mutations,
    hotspot_rows,
    placement_row,
)
from .fasta import read_fasta


def _write_csv(
    path: Path, rows: list[dict[str, object]], fields: list[str]
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare FASTA queries and rank recurrent mutation hotspots."
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA containing exactly one record.",
    )
    parser.add_argument(
        "--queries",
        required=True,
        nargs="+",
        help="One or more query FASTA files.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for CSV and JSON results.",
    )
    parser.add_argument(
        "--max-mismatch-rate",
        type=float,
        default=0.01,
        help="Maximum allowed mismatch rate (default: 0.01).",
    )
    parser.add_argument(
        "--absolute-mismatch-cap",
        type=int,
        default=None,
        help="Optional maximum mismatch count.",
    )
    parser.add_argument(
        "--top-hotspots",
        type=int,
        default=20,
        help="Hotspots retained in run summary (default: 20).",
    )
    args = parser.parse_args()
    if not 0 <= args.max_mismatch_rate <= 1:
        parser.error("--max-mismatch-rate must be between 0 and 1")
    if (
        args.absolute_mismatch_cap is not None
        and args.absolute_mismatch_cap < 0
    ):
        parser.error("--absolute-mismatch-cap must be non-negative")
    if args.top_hotspots < 0:
        parser.error("--top-hotspots must be non-negative")
    return args


def main() -> None:
    args = parse_args()
    references = read_fasta(args.reference)
    if len(references) != 1:
        raise ValueError("Reference FASTA must contain exactly one record")
    reference = references[0]
    queries = [
        record
        for file_name in args.queries
        for record in read_fasta(file_name)
    ]
    if len({record.identifier for record in queries}) != len(queries):
        raise ValueError(
            "Query IDs must be unique; remove duplicate samples first"
        )
    placements = [
        best_placement(
            reference.sequence,
            query,
            args.max_mismatch_rate,
            args.absolute_mismatch_cap,
        )
        for query in queries
    ]
    mutations = [
        mutation
        for query, placement in zip(queries, placements)
        for mutation in call_mutations(reference, query, placement)
    ]
    hotspots = hotspot_rows(mutations)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(
        output_dir / "mutations.csv",
        [item.row() for item in mutations],
        [
            "reference_id",
            "query_id",
            "orientation",
            "reference_position",
            "reference_base",
            "alternate_base",
            "placement_offset",
            "substitution",
        ],
    )
    _write_csv(
        output_dir / "hotspots.csv",
        hotspots,
        [
            "reference_id",
            "reference_position",
            "reference_base",
            "sample_recurrence",
            "event_count",
            "substitution_counts",
        ],
    )
    summary = {
        "reference_id": reference.identifier,
        "reference_length": len(reference.sequence),
        "query_count": len(queries),
        "accepted_query_count": sum(p.accepted for p in placements),
        "mutation_count": len(mutations),
        "parameters": {
            "max_mismatch_rate": args.max_mismatch_rate,
            "absolute_mismatch_cap": args.absolute_mismatch_cap,
        },
        "placements": [placement_row(item) for item in placements],
        "top_hotspots": hotspots[: args.top_hotspots],
    }
    (output_dir / "run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(
        f"Processed {len(queries)} query records; wrote results to {output_dir}"
    )


if __name__ == "__main__":
    main()
