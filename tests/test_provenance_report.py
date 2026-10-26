"""Tests for g2f_provenance_report (synthetic provenance records only)."""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

import util

SCRIPT = os.path.join(util.PIPELINE_DIR, "scripts", "g2f_provenance_report.py")


class ProvenanceReportTest(unittest.TestCase):
    def test_table(self):
        d = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        os.makedirs(os.path.join(d, "SAMPLE_A", "provenance", "hifiasm"))
        rec = {"schema": "g2f-provenance/1", "sample": "SAMPLE_A", "step": "hifiasm", "method_id": "m1",
               "status": "completed", "exit_code": 0, "start": "2000-01-01T00:00:00", "end": None,
               "pipeline_git": {"commit": "abc", "dirty": False},
               "tools": {"hifiasm": {"reported": "0.25.0-r726"}, "seqkit": {"reported": "seqkit v2"}}}
        json.dump(rec, open(os.path.join(d, "SAMPLE_A", "provenance", "hifiasm", "m1.json"), "w"))
        json.dump({"schema": "other"}, open(os.path.join(d, "ignored.json"), "w"))
        open(os.path.join(d, "broken.json"), "w").write("{")
        out = os.path.join(d, "t.tsv")
        r = subprocess.run([sys.executable, SCRIPT, d, "--out", out], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        rows = list(csv.DictReader(open(out), delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["tools_reported"], "hifiasm=0.25.0-r726;seqkit=seqkit v2")
        self.assertEqual(rows[0]["end"], "")
        self.assertEqual(rows[0]["pipeline_dirty"], "False")


if __name__ == "__main__":
    unittest.main()
