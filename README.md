# Viral Variant Atlas

**A reproducible Python CLI for scanning viral gene sequences, identifying nucleotide substitutions, and ranking recurrent mutation hotspots.**

Viral Variant Atlas turns a collection of FASTA sequences into an auditable mutation report. It compares each query sequence against a reference, searches for the best ungapped placement, records substitutions using one-based genomic coordinates, and aggregates recurrent events into a hotspot table.

> This project is intended for exploratory sequence comparison. It uses an ungapped sliding-window matcher and is not a replacement for a validated multiple-sequence-alignment workflow or clinical interpretation.

## Why it is stronger than a one-off script

- **Reproducible CLI:** configuration is explicit and every run writes machine-readable results.
- **Defensive FASTA parsing:** validates sequence alphabet, whitespace, headers, and empty inputs.
- **Traceable calls:** each reported event includes reference/query IDs, placement offset, one-based position, and the observed substitution.
- **Hotspot analytics:** ranks recurrent positions across samples and calculates sample-level recurrence.
- **Portable by default:** built with the Python standard library; no notebook or external dependency required.
- **Tested core logic:** unit tests cover parsing, coordinates, reverse-complement search, thresholds, and aggregation.

## Quick start

Requires Python 3.10+.

From this directory, install the package first: `python -m pip install -e .`.

```bash
python -m viral_variant_atlas \
  --reference data/example/reference_orf7b.fasta \
  --queries data/example/variant_a.fasta data/example/variant_b.fasta \
  --output-dir results \
  --max-mismatch-rate 0.10 \
  --top-hotspots 10
```

The command produces:

```text
results/
├── mutations.csv       # one row per substitution
├── hotspots.csv        # recurrent positions, ranked by sample recurrence
└── run_summary.json    # inputs, parameters, placement metadata, and counts
```

## Method

For every query sequence, the tool evaluates every equally sized reference window across both the forward strand and reverse complement. The minimum-Hamming-distance placement is retained. A sample is accepted when its mismatch count is no greater than both the rate-derived threshold and optional absolute cap. Substitutions are represented as `REF<position>ALT` (for example, `C51T`), where position is **one-based within the reference sequence**.

Hotspots are aggregated by reference position. `sample_recurrence` counts unique query samples carrying a substitution at a position; it avoids inflating a hotspot when duplicate records occur for the same sample.

## Repository layout

```text
src/viral_variant_atlas/   CLI and analysis package
tests/                     unit tests for the analytical core
data/example/              small synthetic FASTA inputs for a smoke test
```

## Example output

```csv
reference_id,query_id,orientation,reference_position,substitution
ORF7b_reference,variant_a,+,5,T5C
ORF7b_reference,variant_b,+,5,T5C
```

```csv
reference_id,reference_position,reference_base,sample_recurrence,substitution_counts
ORF7b_reference,5,T,2,T5C:2
```

## Quality controls and limitations

- The tool rejects sequences containing characters outside the IUPAC DNA alphabet.
- Queries longer than the reference are not silently truncated; they are reported as unplaceable.
- Ties resolve deterministically (lowest offset, then forward orientation).
- Use a gapped aligner when indels, structural variation, or phylogenetic inference are in scope.
- Interpret recurrence alongside sequencing quality, coverage, metadata, and biological context.

## Run tests

```bash
python -m unittest discover -s tests -v
```

## Resume-ready project description

Built a reproducible Python bioinformatics CLI that parses FASTA inputs, performs bidirectional sliding-window sequence comparisons, produces coordinate-correct substitution calls, and ranks recurrent mutation hotspots through CSV/JSON reporting and automated tests.
