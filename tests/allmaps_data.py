"""Synthetic genome, contig assembly and anchor markers for tests/test_allmaps.py (invented data).

Truth: three chromosomes (chrA..chrC as LGs 1..3), cut into contigs with known order and orientation.
  - ctg_chim joins the end of LG1 to the start of LG3 (a chimeric contig: must be flagged, never split)
  - a 6-kb segment of LG1 is copied into LG2, so its markers map twice (must be rejected as non-unique)
  - every 7th marker is listed under two more IDs (redundant records, identical sequence and position)
  - 5 markers carry 7 clustered substitutions (aligned with NM 7; rejected as divergent with max_nm 5)
  - ctg_chim has 40 kb of LG1 and 60 kb of LG3: placed by its LG3 majority, flagged, not split
  - ctg_short (< min_contig_bp) and ctg_nomark (no markers) stay unplaced
Map positions are chromosome coordinates / 1e5, as in the supplied maize anchor set, EXCEPT LG2 80-150 kb (the span
of ctg05), where they are reversed: an "older reference version" inverted relative to the true chromosome.
ref.fa holds the true chromosomes (R1..R3) as a reference sample for the reference-position variant.
"""
import gzip
import os
import random

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def write_fasta(path, recs, width=60):
    with open(path, "w") as fh:
        for n, s in recs:
            fh.write(f">{n}\n" + "\n".join(s[i:i + width] for i in range(0, len(s), width)) + "\n")
    # .fai (uniform line lengths)
    with open(path + ".fai", "w") as fai:
        off = 0
        for n, s in recs:
            off += len(f">{n}\n")
            fai.write(f"{n}\t{len(s)}\t{off}\t{width}\t{width + 1}\n")
            off += len(s) + (len(s) + width - 1) // width


def make(root, seed=7):
    rnd = random.Random(seed)
    L = {"1": 240000, "2": 210000, "3": 180000}
    chrom = {k: "".join(rnd.choice("ACGT") for _ in range(v)) for k, v in L.items()}
    seg = chrom["1"][30000:36000]                     # repeat: LG1 30-36 kb copied into LG2 at 100 kb
    chrom["2"] = chrom["2"][:100000] + seg + chrom["2"][106000:]
    # contigs: (name, lg, start, end, orientation) in true chromosome order
    plan = [("ctg01", "1", 0, 60000, "+"), ("ctg02", "1", 60000, 130000, "-"), ("ctg03", "1", 130000, 200000, "+"),
            ("ctg04", "2", 0, 80000, "-"), ("ctg05", "2", 80000, 150000, "+"), ("ctg06", "2", 150000, 200000, "-"),
            ("ctg07", "3", 60000, 120000, "+"), ("ctg08", "3", 120000, 180000, "-")]
    recs, truth = [], {}
    for n, lg, a, b, o in plan:
        s = chrom[lg][a:b]
        recs.append((n, s if o == "+" else rc(s)))
        truth[n] = (lg, a, b, o)
    chim = chrom["1"][200000:240000] + chrom["3"][0:60000]  # LG1 tail + LG3 head
    recs.append(("ctg_chim", chim))
    truth["ctg_chim"] = ("1+3", 0, 100000, "+")
    recs.append(("ctg_short", chrom["2"][200000:208000]))  # 8 kb of LG2 beyond ctg06, with its own markers
    recs.append(("ctg_nomark", "".join(rnd.choice("ACGT") for _ in range(30000))))
    asm = os.path.join(root, "asm.fa")
    write_fasta(asm, recs)
    # markers: 100-bp windows every 1 kb on every chromosome
    mfa, meta, k = [], [], 0
    divergent = 0
    for lg in ("1", "2", "3"):
        for p in range(750 if lg == "2" else 500, L[lg] - 600, 1000):  # LG2 grid offset: repeat markers differ from LG1
            s = chrom[lg][p:p + 100]
            k += 1
            if divergent < 5 and lg == "3" and p > 130000:
                s = list(s)
                for j in range(84, 98, 2):  # 7 substitutions in the last 16 bp: seeds still match
                    s[j] = {"A": "C", "C": "G", "G": "T", "T": "A"}[s[j]]
                s = "".join(s)
                divergent += 1
            ids = [f"pg_{k}"] + ([f"pg_{k}r1", f"pg_{k}r2"] if k % 7 == 0 else [])
            mp = (80000 + 150000 - p - 100) if (lg == "2" and 80000 <= p < 150000 - 100) else p
            for m in ids:
                mfa.append((m, s))
                meta.append((lg, m, f"{mp / 1e5:.5f}"))
    with gzip.open(os.path.join(root, "markers.fa.gz"), "wt") as fh:
        for m, s in mfa:
            fh.write(f">{m}\n{s}\n")
    with gzip.open(os.path.join(root, "markers.txt.gz"), "wt") as fh:
        for lg, m, p in meta:
            fh.write(f"{lg}\t{m}\t{p}\n")
    ref = os.path.join(root, "ref.fa")
    write_fasta(ref, [(f"R{k}", chrom[k]) for k in ("1", "2", "3")])
    return {"ref_fa": ref, "asm_fa": asm, "markers_fa": os.path.join(root, "markers.fa.gz"),
            "markers_txt": os.path.join(root, "markers.txt.gz"), "truth": truth, "n_records": len(mfa),
            "n_unique": k, "chrom": chrom}
