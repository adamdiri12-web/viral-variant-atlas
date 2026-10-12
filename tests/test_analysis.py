import unittest

from viral_variant_atlas.analysis import (
    best_placement,
    call_mutations,
    hotspot_rows,
    reverse_complement,
)
from viral_variant_atlas.types import FastaRecord


class AnalysisTests(unittest.TestCase):
    def test_calls_one_based_substitution_at_best_offset(self):
        reference = FastaRecord("ref", "GGATGCTACCTT")
        query = FastaRecord("sample", "ATGCTTC")
        placement = best_placement(
            reference.sequence, query, rate=0.20, absolute_cap=None
        )
        self.assertEqual(
            (
                placement.offset,
                placement.orientation,
                placement.mismatch_count,
            ),
            (2, "+", 1),
        )
        self.assertEqual(
            [
                event.substitution
                for event in call_mutations(reference, query, placement)
            ],
            ["A8T"],
        )

    def test_considers_reverse_complement(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        placement = best_placement(
            "CCGCATGG",
            FastaRecord("sample", "ATGC"),
            rate=0,
            absolute_cap=None,
        )
        self.assertEqual(
            (
                placement.orientation,
                placement.offset,
                placement.mismatch_count,
            ),
            ("-", 2, 0),
        )

    def test_hotspots_count_unique_samples(self):
        reference = FastaRecord("ref", "ATGC")
        queries = [FastaRecord("a", "ATAC"), FastaRecord("b", "ATAC")]
        mutations = [
            event
            for query in queries
            for event in call_mutations(
                reference,
                query,
                best_placement(reference.sequence, query, 0.5, None),
            )
        ]
        hotspot = hotspot_rows(mutations)[0]
        self.assertEqual(
            (
                hotspot["reference_position"],
                hotspot["sample_recurrence"],
                hotspot["substitution_counts"],
            ),
            (3, 2, "G3A:2"),
        )

    def test_rejects_over_threshold_placement(self):
        placement = best_placement(
            "AAAA",
            FastaRecord("sample", "CCCC"),
            rate=0.25,
            absolute_cap=None,
        )
        self.assertFalse(placement.accepted)

    def test_reverse_complement_exact_match_is_accepted(self):
        result = best_placement(
            "AAAA", FastaRecord("sample", "TTTT"), 0, None
        )
        self.assertTrue(result.accepted)
        self.assertEqual(result.orientation, "-")

    def test_tie_prefers_lowest_offset_and_forward(self):
        result = best_placement(
            "ATATAT", FastaRecord("sample", "ATAT"), 0, None
        )
        self.assertEqual((result.offset, result.orientation), (0, "+"))

    def test_long_query_is_unplaceable(self):
        result = best_placement(
            "ATGC", FastaRecord("sample", "ATGCA"), 0, None
        )
        self.assertEqual((result.offset, result.accepted), (-1, False))

    def test_absolute_cap_rejects_otherwise_accepted_query(self):
        result = best_placement("AAAA", FastaRecord("sample", "AACC"), 1, 1)
        self.assertFalse(result.accepted)


if __name__ == "__main__":
    unittest.main()
