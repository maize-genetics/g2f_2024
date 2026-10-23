"""Tests for g2f_read_forensics (synthetic genomes and reads only).

Integration tests need pysam and bcftools (BCFTOOLS or the PATH); they are
skipped otherwise.
"""
import gzip
import os
import random
import shutil
import tempfile
import unittest

import util  # noqa: F401  (puts scripts/ on sys.path)
import g2f_read_forensics as rf

try:
    import pysam
except ImportError:          # pragma: no cover
    pysam = None

BCFTOOLS = os.environ.get("BCFTOOLS") or shutil.which("bcftools") \
    or ("/programs/bcftools-1.20/bin/bcftools" if os.path.exists("/programs/bcftools-1.20/bin/bcftools") else None)


class UnitTest(unittest.TestCase):
    def test_read_alleles_cigar(self):
        # ref 100..: 5M 2I 3M 2D 4M ; read seq indices
        seq = "AAAAA" + "GG" + "CCC" + "TTTT"
        cig = [(0, 5), (1, 2), (0, 3), (2, 2), (0, 4)]
        sites = [100, 104, 105, 107, 108, 109, 110, 113]
        out = rf.read_alleles(100, cig, seq, None, sites)
        self.assertEqual(out[100][0], "A")
        self.assertEqual(out[104][0], "A")
        self.assertEqual(out[105][0], "C")       # insertion skipped
        self.assertEqual(out[107][0], "C")
        self.assertEqual(out[108][0], "-")       # deletion
        self.assertEqual(out[109][0], "-")
        self.assertEqual(out[110][0], "T")
        self.assertEqual(out[113][0], "T")

    def test_cluster_and_count(self):
        cl = rf.cluster_sites([1, 5, 9, 12, 100, 5000, 5001, 5002, 5003, 9000], 4, 50)
        self.assertEqual(cl, [[1, 5, 9, 12], [5000, 5001, 5002, 5003]])
        n, sup = rf.count_haplotypes(["01", "01", "01", "10", "10", "10", "11"], 3, 0.1)
        self.assertEqual((n, sup), (2, [3, 3]))
        n, _ = rf.count_haplotypes(["01"] * 5 + ["10"] * 5 + ["00"] * 4, 3, 0.1)
        self.assertEqual(n, 3)

    def test_mixture_bic(self):
        rng = random.Random(1)

        def draw(p, n):
            return (sum(rng.random() < p for _ in range(n)), n)
        one = [draw(0.48, 34) for _ in range(3000)]
        one = [x for x in one if 0.1 <= x[0] / x[1] <= 0.9]
        two = [draw(0.25, 34) for _ in range(1500)] + [draw(0.5, 34) for _ in range(1500)]
        two = [x for x in two if 0.1 <= x[0] / x[1] <= 0.9]

        def best(d):
            fits = []
            for k in (1, 2, 3):
                ll, comp = rf.fit_binom_mixture(d, k)
                fits.append((-2 * ll + (2 * k - 1) * __import__("math").log(len(d)), k, comp))
            return min(fits)
        b1 = best(one)
        self.assertEqual(b1[1], 1)
        self.assertAlmostEqual(b1[2][0][1], 0.48, delta=0.02)
        b2 = best(two)
        self.assertGreaterEqual(b2[1], 2)
        ps = sorted(p for w, p in b2[2] if w > 0.2)
        self.assertAlmostEqual(ps[0], 0.25, delta=0.04)

    def test_parse_ad_and_classify(self):
        self.assertEqual(rf.parse_ad("C,T,<*>", "10,5,1,0"), {"REF": 10, "C": 5, "T": 1})
        self.assertEqual(rf.classify_hifi_site(15, 15, 20, 40, 0.2, 0.8), "het")
        self.assertEqual(rf.classify_hifi_site(1, 29, 20, 40, 0.2, 0.8), "homalt")
        self.assertEqual(rf.classify_hifi_site(5, 2, 20, 40, 0.2, 0.8), "offdepth")

    def test_windows(self):
        self.assertEqual(rf.make_windows([("a", 2500), ("b", 900)], 500, 1000),
                         [("a", 0, 500), ("a", 1000, 1500), ("a", 2000, 2500), ("b", 0, 500)])
        self.assertEqual(rf.merge_intervals([("a", 10, 20), ("a", 15, 30), ("a", 40, 50)], 5),
                         [("a", 10, 30), ("a", 40, 50)])


def _mut(seq, sites, rng):
    s = list(seq)
    for p in sites:
        s[p] = rng.choice([b for b in "ACGT" if b != s[p]])
    return "".join(s)


def _md(ref, read):
    out, run = [], 0
    for a, b in zip(ref, read):
        if a == b:
            run += 1
        else:
            out.append(str(run) + a)
            run = 0
    out.append(str(run))
    return "".join(out)


def _write_bam(path, ref_name, ref, sources, n_reads, rlen, rng, err=0.0, qual=40):
    """sources: list of (haplotype sequence, weight). M-only alignments + MD."""
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": ref_name, "LN": len(ref)}]}
    recs = []
    tot = sum(w for _, w in sources)
    for i in range(n_reads):
        x = rng.random() * tot
        for hap, w in sources:
            x -= w
            if x <= 0:
                break
        st = rng.randrange(0, len(ref) - rlen)
        s = list(hap[st:st + rlen])
        for j in range(rlen):
            if err and rng.random() < err:
                s[j] = rng.choice([b for b in "ACGT" if b != s[j]])
        recs.append((st, "".join(s), f"r{i}"))
    recs.sort()
    tmp = path + ".unsorted.bam"
    with pysam.AlignmentFile(tmp, "wb", header=hdr) as out:
        for st, s, name in recs:
            a = pysam.AlignedSegment()
            a.query_name = name
            a.query_sequence = s
            a.flag = 0
            a.reference_id = 0
            a.reference_start = st
            a.mapping_quality = 60
            a.cigartuples = [(0, len(s))]
            a.query_qualities = pysam.qualitystring_to_array(chr(33 + qual) * len(s))
            a.set_tag("MD", _md(ref[st:st + len(s)], s))
            out.write(a)
    pysam.sort("-o", path, tmp)
    pysam.index(path)
    os.remove(tmp)


@unittest.skipIf(pysam is None or BCFTOOLS is None, "needs pysam and bcftools")
class IntegrationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="rf_test.", dir=os.environ.get("TMPDIR"))
        rng = random.Random(7)
        L = 60000
        ref = "".join(rng.choice("ACGT") for _ in range(L))
        cls.ref_fa = os.path.join(cls.tmp, "ref.fa")
        with open(cls.ref_fa, "w") as fh:
            fh.write(">chr1\n")
            for i in range(0, L, 60):
                fh.write(ref[i:i + 60] + "\n")
        pysam.faidx(cls.ref_fa)
        h1 = ref                                              # reference haplotype
        v2 = sorted(rng.sample(range(100, L - 100), 300))
        h2 = _mut(ref, v2, rng)
        # h3 shares half of h2's variants (same allele): sibling-like plants
        shared = set(v2[::2])
        h3 = "".join(h2[i] if i in shared else c for i, c in enumerate(ref))
        h3 = _mut(h3, sorted(rng.sample([i for i in range(100, L - 100) if i not in set(v2)], 150)), rng)
        cls.dip = os.path.join(cls.tmp, "diploid.bam")        # one plant: h1/h2
        _write_bam(cls.dip, "chr1", ref, [(h1, 1), (h2, 1)], 400, 3000, rng)
        cls.pool = os.path.join(cls.tmp, "pool.bam")          # two plants: h1/h2 + h1/h3
        _write_bam(cls.pool, "chr1", ref, [(h1, 2), (h2, 1), (h3, 1)], 600, 3000, rng)
        cls.ont = os.path.join(cls.tmp, "ont_h3.bam")         # other-plant long reads: h3
        _write_bam(cls.ont, "chr1", ref, [(h3, 1)], 6, 20000, rng, err=0.001)
        cls.ont_same = os.path.join(cls.tmp, "ont_h2.bam")
        _write_bam(cls.ont_same, "chr1", ref, [(h2, 1)], 6, 20000, rng, err=0.001)
        cls.bed = os.path.join(cls.tmp, "all.bed")
        with open(cls.bed, "w") as fh:
            fh.write(f"chr1\t0\t{L}\n")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def _snv(self, bam, name):
        out = os.path.join(self.tmp, name + ".snv.tsv.gz")
        rf.main(["snv", "--aln", bam, "--ref", self.ref_fa, "--bed", self.bed, "--bcftools", BCFTOOLS,
                 "--threads", "2", "--tmpdir", self.tmp, "--out", out])
        return out

    def test_diploid_vs_pool(self):
        import json
        res = {}
        for name, bam, depth in (("dip", self.dip, 20), ("pool", self.pool, 30)):
            snv = self._snv(bam, name)
            od = os.path.join(self.tmp, name + ".af")
            rf.main(["af", "--sites", snv, "--expected-depth", str(depth), "--outdir", od,
                     "--window", "20000"])
            with open(os.path.join(od, "af_summary.json")) as fh:
                s = json.load(fh)
            hp = os.path.join(self.tmp, name + ".haps.tsv")
            rf.main(["haps", "--aln", bam, "--ref", self.ref_fa, "--sites", snv, "--expected-depth",
                     str(depth), "--max-span", "1500", "--out", hp])
            with open(hp) as fh:
                rows = [ln.split("\t") for ln in fh][1:]
            res[name] = (s, [int(r[6]) for r in rows])
        s, nh = res["dip"]
        self.assertEqual(s["best_k_by_bic"], 1)
        self.assertAlmostEqual(s["mixture_fits"][0]["components"][0]["p"], 0.5, delta=0.04)
        self.assertTrue(nh and max(nh) <= 2)
        s, nh = res["pool"]
        # pool: private variants at 0.25, shared ones at 0.5
        self.assertGreaterEqual(s["best_k_by_bic"], 2)
        fit = s["mixture_fits"][s["best_k_by_bic"] - 1]["components"]
        self.assertTrue(any(abs(c["p"] - 0.25) < 0.05 and c["weight"] > 0.3 for c in fit))
        self.assertGreater(sum(1 for x in nh if x >= 3), 0)

    def test_ont_other_plant(self):
        snv = self._snv(self.dip, "dip2")
        rates = {}
        for name, bam in (("h3", self.ont), ("h2", self.ont_same)):
            od = os.path.join(self.tmp, name + ".obs")
            rf.main(["ontobs", "--aln", bam, "--ref", self.ref_fa, "--sites", snv, "--outdir", od,
                     "--min-aligned", "1000"])
            pos = os.path.join(od, "subs.pos")
            with gzip.open(os.path.join(od, "obs.tsv.gz"), "rt") as fh, open(pos, "w") as fo:
                next(fh)
                for line in fh:
                    f = line.split("\t")
                    if f[6].strip() == "sub":
                        fo.write(f"{f[1]}\t{f[2]}\n")
            lk = os.path.join(od, "lookup.tsv.gz")
            rf.main(["lookup", "--aln", self.dip, "--ref", self.ref_fa, "--bed", self.bed, "--positions", pos,
                     "--bcftools", BCFTOOLS, "--threads", "2", "--tmpdir", self.tmp, "--out", lk])
            sd = os.path.join(self.tmp, name + ".sum")
            rf.main(["ontsum", "--sites", snv, "--lookup", lk, "--obsdir", od, "--expected-depth", "20",
                     "--outdir", sd])
            with gzip.open(os.path.join(sd, "reads_summary.tsv.gz"), "rt") as fh:
                hdr = next(fh).rstrip("\n").split("\t")
                rows = [dict(zip(hdr, ln.rstrip("\n").split("\t"))) for ln in fh]
            rates[name] = (sum(int(r["sub_only"]) for r in rows), sum(int(r["het_alt"]) for r in rows),
                           sum(int(r["het_ref"]) for r in rows))
        # h3 reads carry ~300/60 kb alleles absent from the h1/h2 library; h2 reads only errors
        self.assertGreater(rates["h3"][0], 3 * rates["h2"][0] + 20)
        # at h1/h2 heterozygous sites h2 reads carry the alt allele; h3 shares half of them
        self.assertGreater(rates["h2"][1], 10 * max(1, rates["h2"][2]))
        f3 = rates["h3"][1] / (rates["h3"][1] + rates["h3"][2])
        self.assertTrue(0.3 < f3 < 0.7, f3)


if __name__ == "__main__":
    unittest.main()
