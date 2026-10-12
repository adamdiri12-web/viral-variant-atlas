import tempfile
import unittest
from pathlib import Path

from viral_variant_atlas.fasta import FastaFormatError, read_fasta


class FastaTests(unittest.TestCase):
    def parse(self, text):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.fasta"
            path.write_text(text, encoding="utf-8")
            return read_fasta(path)

    def test_multiple_records_and_whitespace(self):
        records = self.parse(">a description\nat g\tc\n>b\nNNRY\n")
        self.assertEqual(
            [(r.identifier, r.sequence) for r in records],
            [("a", "ATGC"), ("b", "NNRY")],
        )

    def test_malformed_inputs_are_rejected(self):
        for text in ("", ">\nACGT", ">a\n", "ACGT", ">a\nACGT!"):
            with self.subTest(text=text), self.assertRaises(
                FastaFormatError
            ):
                self.parse(text)
