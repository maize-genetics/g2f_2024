"""Tests for the supplementary fingerprinting diagnostics (synthetic hash sets only)."""
import random
import unittest

import util  # noqa: F401  (puts scripts/ on sys.path)
from g2f import fingerprint as fp
import g2f_fp_diagnostics as dg


def _counts(kmers, depth, rng):
    """Poisson-like abundances of `kmers` sequenced at `depth` (absent if 0)."""
    out, n = {}, 40
    for h in kmers:
        c = sum(1 for _ in range(n) if rng.random() < depth / n)
        if c:
            out[h] = c
    return out


class DiagnosticsTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = random.Random(5)
        base = list(range(1, 30001))
        cls.geno, nxt = {}, 100000
        for g in ("A", "B", "C", "D"):
            ks = set(rng.sample(base, 24000))       # shared ancestry
            ks |= set(range(nxt, nxt + 1500))       # private sequence
            nxt += 10000
            cls.geno[g] = ks
        P = {"min_count_floor": 2, "max_copy": 3.0}
        refs, counts_of = {}, {}
        for g, ks in cls.geno.items():
            counts_of[g] = _counts(ks, 12, rng)
            refs[g] = fp.Reference(g, counts_of[g], P, ["r1"])
        cls.panel = fp.Panel(refs)
        cls.solid_n = {}
        for g in cls.panel.genotypes:
            for h in refs[g].solid:
                cls.solid_n[h] = cls.solid_n.get(h, 0) + 1
        cls.classes = dg.marker_classes(cls.panel, counts_of, cls.solid_n, "A")
        cls.strata = dg.line_strata(cls.panel, cls.solid_n, cls.classes, "A")
        cls.rng = rng

    def _excess(self, q, versus):
        rows = dg.line_specific("q", "run", "A", q, self.panel, self.strata)
        return next(r for r in rows if r["versus"] == versus)

    def test_zero_truncated_poisson_inverse(self):
        for lam in (0.3, 1.0, 2.5, 6.0):
            m = lam / (1 - 2.718281828459045 ** -lam)
            self.assertAlmostEqual(dg.zt_lambda(m), lam, places=4)
        self.assertEqual(dg.zt_lambda(1.0), 0.0)

    def test_pure_query_has_no_line_specific_excess(self):
        q = _counts(self.geno["A"], 4, self.rng)
        for other in "BCD":       # no errors simulated: no foreign k-mers at all
            self.assertEqual(self._excess(q, other)["excess_hits"], 0)

    def test_low_level_admixture_has_low_dosage(self):
        q = _counts(self.geno["A"], 4, self.rng)
        fp.add_counts(q, _counts(self.geno["B"] - self.geno["A"], 0.4, self.rng))  # ~10 %
        r = self._excess(q, "B")
        self.assertGreater(r["excess_hits"], 300)
        self.assertLess(r["excess_dose_vs_core"], 0.3)
        # lines sharing alleles with the true source show a smaller "shadow" excess
        self.assertLess(self._excess(q, "C")["excess_hits"], r["excess_hits"] / 2)

    def test_regional_haplotype_has_genomic_dosage(self):
        # the query carries B's alleles at full depth, but only at a subset of loci
        b_only = sorted(self.geno["B"] - self.geno["A"])
        part = set(self.rng.sample(b_only, len(b_only) // 5))
        q = _counts(self.geno["A"] | part, 4, self.rng)
        r = self._excess(q, "B")
        self.assertGreater(r["excess_hits"], 300)
        self.assertGreater(r["excess_dose_vs_core"], 0.6)


if __name__ == "__main__":
    unittest.main()
