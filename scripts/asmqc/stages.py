"""Stage definitions: what each stage runs, needs, costs and validates."""
import gzip
import json
import os
import re
import subprocess

from . import graph as G
from . import pairalign, plots, readsupport, refalign, regions, evidence
from .config import tool
from .runner import Stage, StageError, copy_if_exists, readtext

TYPE_ORDER = ["prepare", "graph", "ref_index", "ref_align", "ref_analysis", "pair_align", "map_hifi",
              "read_support", "regions", "evidence", "pair_analysis"]
HEAVY = {"map_hifi", "ref_index", "ref_align", "pair_align"}


# ---------------------------------------------------------------- helpers
def asm_dir(cfg, s, a, sub):
    return os.path.join(cfg["root"], "samples", s, a, sub)


def ref_dir(cfg, r):
    return os.path.join(cfg["root"], "references", r)


def pair_dir(cfg, s, a, b, sub):
    return os.path.join(cfg["root"], "samples", s, "pairs", f"{a}_vs_{b}", sub)


def fasta_link_name(path):
    return "asm.fa.gz" if path.endswith(".gz") else "asm.fa"


def link_fasta(src, workdir, name):
    """Symlink a FASTA (and existing .fai/.gzi) into workdir; index if needed.
    Never writes next to the original."""
    dst = os.path.join(workdir, name)
    os.symlink(os.path.realpath(src), dst)
    for ext in (".fai", ".gzi"):
        if os.path.exists(src + ext):
            os.symlink(os.path.realpath(src + ext), dst + ext)
    return dst


def read_fai(path):
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.split("\t")
            out.append((p[0], int(p[1])))
    return out


def nstats(lengths):
    ls = sorted(lengths, reverse=True)
    tot = sum(ls)
    res = {"bp": tot, "contigs": len(ls), "largest": ls[0] if ls else 0}
    for q in (50, 90):
        acc = 0
        for i, l in enumerate(ls, 1):
            acc += l
            if acc >= tot * q / 100:
                res[f"N{q}"], res[f"L{q}"] = l, i
                break
    return res


ASM_LOG_PATTERNS = {
    "version": r"\[M::main\] Version: (\S+)",
    "command": r"\[M::main\] CMD: (.+)",
    "peak_hom": r"peak_hom: (\d+)",
    "peak_het": r"peak_het: (\d+)",
    "hom_cov_threshold": r"homozygous read coverage threshold: (\d+)",
    "purge_cov_threshold": r"purge duplication coverage threshold: (\d+)",
    "hifi_bases": r"\[M::pec\] # bases: (\d+)",
}


def parse_assembler_log(path):
    out = {}
    if not path or not os.path.exists(path):
        return out
    text = readtext(path, errors="replace")
    for k, rx in ASM_LOG_PATTERNS.items():
        m = re.findall(rx, text)
        if m:
            out[k] = m[-1] if k in ("version", "command") else sorted(set(m), key=m.index)
    return out


def _gb(path):
    try:
        return os.path.getsize(path) / 1e9
    except OSError:
        return 0.0


def _asm_gb(a):
    # bgzip FASTA is ~3.5x smaller than the sequence
    return _gb(a["fasta"]) * (3.5 if a["fasta"].endswith(".gz") else 1.0)


# ---------------------------------------------------------------- builders
def build(cfg, only_samples=None):
    sci = cfg["science"]
    st = {}

    def add(stage):
        st[stage.sid] = stage

    for rname, r in cfg["references"].items():
        d = ref_dir(cfg, rname)
        add(Stage(f"ref_index:{rname}", d, run=_run_ref_index(cfg, r), inputs=lambda r=r: [r["fasta"]],
                  params=lambda: {"preset": sci["ref_align"]["preset"]}, heavy=True,
                  estimate=lambda ctx, r=r: {"threads": ctx.res["threads"], "ram_gb": 6 * _gb(r["fasta"]) + 10,
                                             "disk_final_gb": 4 * _gb(r["fasta"]), "hours": 0.3},
                  describe="minimap2 index + reference N-gap bins"))

    for sname, s in cfg["samples"].items():
        if only_samples and sname not in only_samples:
            continue
        for aname, a in s["assemblies"].items():
            P = asm_dir(cfg, sname, aname, "prepare")
            add(Stage(f"prepare:{sname}:{aname}", P, run=_run_prepare(cfg, s, a),
                      inputs=lambda a=a: [a["fasta"]] + ([a["assembler_log"]] if a["assembler_log"] else []),
                      describe="link FASTA, index, contig N50/L50, assembler log"))
            ref_fa = os.path.join(P, fasta_link_name(a["fasta"]))
            M = asm_dir(cfg, sname, aname, "map_hifi")
            aln = os.path.join(M, f"{sname}.{aname}.hifi.{cfg['resources']['alignment_format']}")
            add(Stage(f"map_hifi:{sname}:{aname}", M, deps=[f"prepare:{sname}:{aname}"], heavy=True,
                      run=_run_map(cfg, s, a, ref_fa), inputs=lambda s=s, a=a: s["hifi_reads"] + [a["fasta"]],
                      params=lambda: sci["mapping"], estimate=_est_map(s, a),
                      describe="minimap2 map-hifi -> sorted, indexed alignments"))
            R = asm_dir(cfg, sname, aname, "read_support")
            add(Stage(f"read_support:{sname}:{aname}", R, deps=[f"map_hifi:{sname}:{aname}"],
                      run=_run_read_support(cfg, aln, ref_fa, P), inputs=lambda aln=aln: [aln],
                      params=lambda: sci["read_support"],
                      estimate=lambda ctx: {"threads": ctx.res["threads"], "ram_gb": 2 * ctx.res["threads"],
                                            "hours": 1.0},
                      describe="depth / continuity / clipping / mixed sites per window"))
            Rg = asm_dir(cfg, sname, aname, "regions")
            add(Stage(f"regions:{sname}:{aname}", Rg, deps=[f"read_support:{sname}:{aname}"],
                      run=_run_regions(cfg, s, R, P), inputs=lambda R=R: [os.path.join(R, "VALIDATED.json")],
                      params=lambda s=s: {"regions": sci["regions"], "expected_depth": s["expected_depth"],
                                          "kmer_peak": s["kmer_peak"]},
                      describe="classify windows and breakpoints"))
            if a["ctg_gfa"] and a["utg_gfa"]:
                Gd = asm_dir(cfg, sname, aname, "graph")
                add(Stage(f"graph:{sname}:{aname}", Gd, deps=[f"prepare:{sname}:{aname}"],
                          run=_run_graph(s, a, P), inputs=lambda a=a: [a["ctg_gfa"], a["utg_gfa"]],
                          params=lambda s=s: {"kmer_peak": s["kmer_peak"]},
                          estimate=lambda ctx: {"ram_gb": 8, "hours": 0.2},
                          describe="unitig topology, contig->unitig paths"))
            ref_dirs = {}
            for rname, r in cfg["references"].items():
                Ad = asm_dir(cfg, sname, aname, f"ref_{rname}_align")
                add(Stage(f"ref_align:{sname}:{aname}:{rname}", Ad, heavy=True,
                          deps=[f"prepare:{sname}:{aname}", f"ref_index:{rname}"],
                          run=_run_ref_align(cfg, r, ref_fa, sname, aname),
                          inputs=lambda a=a, r=r: [a["fasta"], r["fasta"]],
                          params=lambda: {"preset": sci["ref_align"]["preset"]},
                          estimate=lambda ctx, a=a, r=r: {"threads": ctx.res["threads"],
                                                          "ram_gb": 8 * _gb(r["fasta"]) + 3 * _asm_gb(a) + 10,
                                                          "disk_final_gb": 0.5, "hours": 1.0},
                          describe=f"minimap2 {sci['ref_align']['preset']} vs {rname}"))
                An = asm_dir(cfg, sname, aname, f"ref_{rname}_analysis")
                ref_dirs[rname] = An
                add(Stage(f"ref_analysis:{sname}:{aname}:{rname}", An,
                          deps=[f"ref_align:{sname}:{aname}:{rname}"],
                          run=_run_ref_analysis(cfg, r, Ad, P, sname, a),
                          inputs=lambda Ad=Ad: [os.path.join(Ad, "VALIDATED.json")],
                          params=lambda: sci["ref_align"], describe="assignment, junctions, dot plot"))
            Ev = asm_dir(cfg, sname, aname, "evidence")
            deps = [f"regions:{sname}:{aname}"] + [f"ref_analysis:{sname}:{aname}:{r}" for r in cfg["references"]]
            gdir = asm_dir(cfg, sname, aname, "graph") if (a["ctg_gfa"] and a["utg_gfa"]) else None
            if gdir:
                deps.append(f"graph:{sname}:{aname}")
            add(Stage(f"evidence:{sname}:{aname}", Ev, deps=deps,
                      run=_run_evidence(cfg, a, Rg, R, ref_dirs, gdir, f"{sname}.{aname}"),
                      inputs=lambda deps=deps: [os.path.join(st[d].outdir, "VALIDATED.json") for d in deps],
                      params=lambda: {"evidence": sci["evidence"], "graph": sci["graph"]},
                      describe="ranked candidates + graph neighbourhoods"))
        for pr in s["pairs"]:
            A, B = s["assemblies"][pr["a"]], s["assemblies"][pr["b"]]
            pa = pair_dir(cfg, sname, pr["a"], pr["b"], "align")
            extra = B["extra_targets"] if pr["use_b_extra_targets"] else []
            add(Stage(f"pair_align:{sname}:{pr['a']}:{pr['b']}", pa, heavy=True,
                      deps=[f"prepare:{sname}:{pr['a']}", f"prepare:{sname}:{pr['b']}"],
                      run=_run_pair_align(cfg, sname, pr, A, B, extra),
                      inputs=lambda A=A, B=B, extra=extra: [A["fasta"], B["fasta"]] + extra,
                      params=lambda: sci["pair_align"],
                      estimate=lambda ctx, A=A, B=B: {"threads": ctx.res["threads"],
                                                      "ram_gb": 12 * _asm_gb(B) + 4 * _asm_gb(A) + 10,
                                                      "disk_final_gb": 1, "hours": 2.0},
                      describe="minimap2 asm5 (base-level) A -> B and A -> B haplotypes"))
            pn = pair_dir(cfg, sname, pr["a"], pr["b"], "analysis")
            deps = [f"pair_align:{sname}:{pr['a']}:{pr['b']}", f"regions:{sname}:{pr['a']}",
                    f"regions:{sname}:{pr['b']}"]
            have_graph = all(x["ctg_gfa"] and x["utg_gfa"] for x in (A, B))
            if have_graph:
                deps += [f"graph:{sname}:{pr['a']}", f"graph:{sname}:{pr['b']}"]
            add(Stage(f"pair_analysis:{sname}:{pr['a']}:{pr['b']}", pn, deps=deps,
                      run=_run_pair_analysis(cfg, sname, pr, A, B, pa, have_graph),
                      inputs=lambda deps=deps: [os.path.join(st[d].outdir, "VALIDATED.json") for d in deps],
                      params=lambda: sci["pair_align"],
                      describe="sequence classes, depth cross-tab, read retention"))
    return st


def order(stages):
    def key(sid):
        return (TYPE_ORDER.index(sid.split(":")[0]), sid)
    return sorted(stages, key=key)


# ---------------------------------------------------------------- runners
def _run_prepare(cfg, s, a):
    def run(ctx, wd):
        if not os.path.exists(a["fasta"]):
            raise StageError(f"assembly FASTA missing: {a['fasta']}")
        name = fasta_link_name(a["fasta"])
        dst = link_fasta(a["fasta"], wd, name)
        if a["fasta"].endswith(".gz") and not os.path.exists(a["fasta"] + ".gzi"):
            raise StageError("compressed FASTA needs bgzip + .fai + .gzi next to it (samtools faidx)")
        if not os.path.exists(dst + ".fai"):
            ctx.sh(f"{tool(cfg, 'samtools')} faidx {name}", wd, pin=False)
        contigs = read_fai(dst + ".fai")
        if not contigs:
            raise StageError("empty assembly")
        stats = nstats([l for _, l in contigs])
        stats["assembler_log"] = parse_assembler_log(a["assembler_log"])
        stats["label"], stats["method_id"] = a["label"], a["method_id"]
        with open(os.path.join(wd, "stats.json"), "w") as fh:
            json.dump(stats, fh, indent=2)
        with open(os.path.join(wd, "contigs.tsv"), "w") as fh:
            fh.write("contig\tlength\n")
            for c, l in contigs:
                fh.write(f"{c}\t{l}\n")
        copy_if_exists(a["assembler_log"], os.path.join(wd, "assembler.stderr.log"))
        return {"outputs": ["stats.json", "contigs.tsv"], "validation": {"contigs": len(contigs), "bp": stats["bp"]}}
    return run


def _est_map(s, a):
    def est(ctx):
        reads_gb = sum(_gb(p) for p in s["hifi_reads"])
        fmt = ctx.res["alignment_format"]
        return {"threads": ctx.res["threads"] + 4 + ctx.res["sort_threads"],
                "ram_gb": 15 + 6 * _asm_gb(a) + ctx.res["sort_threads"] * _mem_gb(ctx.res["sort_mem_per_thread"]) * 1.3,
                "disk_final_gb": reads_gb * (0.3 if fmt == "cram" else 0.75),
                "disk_tmp_gb": reads_gb * 0.9,
                "hours": round(reads_gb / 38 * 3.0 * 32 / max(8, ctx.res["threads"]), 1)}
    return est


def _mem_gb(x):
    m = re.match(r"(\d+)([KMG]?)", str(x))
    v, u = int(m.group(1)), m.group(2)
    return v / {"K": 1e6, "M": 1e3, "G": 1}.get(u, 1e9)


def read_stream_cmd(cfg, path, drop_quals=True):
    """Command writing one read file to stdout (FASTA if drop_quals, else as is)."""
    if path.endswith((".bam", ".cram")):
        return f"{tool(cfg, 'samtools')} {'fasta' if drop_quals else 'fastq'} -@ 4 {path} 2>/dev/null"
    if re.search(r"\.(fastq|fq)(\.gz)?$", path) and drop_quals:
        return f"{tool(cfg, 'seqkit')} fq2fa -j 4 {path}"
    return f"zcat -f {path}"


def _run_map(cfg, s, a, ref_fa):
    def run(ctx, wd):
        sci, res = ctx.sci["mapping"], ctx.res
        fmt = res["alignment_format"]
        out = f"{s['name']}.{a['name']}.hifi.{fmt}"
        src = "{ " + " ".join(read_stream_cmd(cfg, r, sci["drop_base_qualities"]) + ";" for r in s["hifi_reads"]) + " }"
        rg = f"'@RG\\tID:{s['name']}\\tSM:{s['name']}\\tPL:PACBIO'"
        tmp = os.path.join(cfg["tmpdir"], f"sort.{s['name']}.{a['name']}.{os.getpid()}")
        os.makedirs(tmp, exist_ok=True)
        ofmt = "-O CRAM --output-fmt-option version=3.0" if fmt == "cram" else "-O BAM"
        cmd = (f"{src} | {tool(cfg, 'minimap2')} -t {res['threads']} -a -x {sci['preset']} "
               f"{' '.join(sci['extra_args'])} -R {rg} {ref_fa} - 2> minimap2.log "
               f"| {tool(cfg, 'samtools')} sort -@ {res['sort_threads']} -m {res['sort_mem_per_thread']} "
               f"-T {tmp}/part --reference {ref_fa} {ofmt} -o {out} -")
        ctx.sh(cmd, wd)
        ctx.sh(f"{tool(cfg, 'samtools')} index -@ 4 {out}", wd)
        ctx.sh(f"{tool(cfg, 'samtools')} quickcheck {out}", wd)
        ctx.sh(f"{tool(cfg, 'samtools')} flagstat -@ 8 -O json --input-fmt-option reference={ref_fa} {out} > flagstat.json", wd)
        if os.path.isdir(tmp) and not os.listdir(tmp):
            os.rmdir(tmp)
        with open(os.path.join(wd, "flagstat.json")) as fh:
            fs = json.load(fh)["QC-passed reads"]
        val = {"primary": fs["primary"], "primary_mapped": fs["primary mapped"],
               "primary_mapped_pct": fs["primary mapped %"], "supplementary": fs["supplementary"],
               "secondary": fs["secondary"]}
        exp = s.get("hifi_read_count")
        if exp is not None:
            if int(fs["primary"]) != int(exp):
                raise StageError(f"primary records {fs['primary']} != expected HiFi read count {exp}")
            val["read_count_check"] = "primary records == configured HiFi read count"
        else:
            val["read_count_check"] = "not checked (hifi_read_count not configured)"
        idx = out + (".crai" if fmt == "cram" else ".bai")
        return {"outputs": [out, os.path.basename(idx), "flagstat.json", "minimap2.log"], "validation": val}
    return run


def _contigs(prep_dir):
    return [(r.split("\t")[0], int(r.split("\t")[1])) for r in readtext(os.path.join(prep_dir, "contigs.tsv")).splitlines()[1:]]


def _run_read_support(cfg, aln, ref_fa, prep_dir):
    def run(ctx, wd):
        contigs = _contigs(prep_dir)
        cpus = ctx.res.get("cpus")
        if cpus:  # pin the Python workers like external tools
            os.sched_setaffinity(0, _cpu_list(cpus))
        summ = readsupport.run(aln, ref_fa, contigs, ctx.sci["read_support"], ctx.res["threads"], wd + "/")
        if summ["assembly_bp"] != sum(l for _, l in contigs):
            raise StageError(f"windows cover {summ['assembly_bp']} bp != assembly {sum(l for _, l in contigs)} bp")
        return {"outputs": ["windows.tsv.gz", "clip_clusters.tsv.gz", "depth_hist.tsv", "summary.json"],
                "validation": {"windows": summ["windows"], "depth_mode_all": summ["depth_mode_all"],
                               "covers_all_bases": True}}
    return run


def _cpu_list(spec):
    out = set()
    for part in str(spec).split(","):
        if "-" in part:
            a, b = part.split("-")
            out.update(range(int(a), int(b) + 1))
        elif part.strip():
            out.add(int(part))
    return out


def _run_regions(cfg, s, read_dir, prep_dir):
    def run(ctx, wd):
        lengths = dict(_contigs(prep_dir))
        summ = regions.run(os.path.join(read_dir, "windows.tsv.gz"), os.path.join(read_dir, "clip_clusters.tsv.gz"),
                           os.path.join(read_dir, "summary.json"), lengths, ctx.sci, wd,
                           kmer_peak=s["kmer_peak"], expected_setting=s["expected_depth"])
        hist = [l.split("\t") for l in readtext(os.path.join(read_dir, "depth_hist.tsv")).splitlines()[1:]]
        plots.depth_histograms([(s["name"], [int(h[0]) for h in hist], [int(h[1]) for h in hist],
                                 summ["expected_depth"])], os.path.join(wd, "depth_histogram"))
        plots.class_bars([(s["name"], summ["class_bp"])], os.path.join(wd, "depth_classes"))
        return {"outputs": ["regions.tsv", "breakpoints.tsv", "window_classes.tsv.gz", "summary.json",
                            "depth_histogram.svg"],
                "validation": {"expected_depth": summ["expected_depth"], "warning": summ["warning"]}}
    return run


def _peak_hom(prep_dir):
    with open(os.path.join(prep_dir, "stats.json")) as fh:
        ph = json.load(fh)["assembler_log"].get("peak_hom")
    return int(ph[-1]) if ph else None  # last estimate (ha_pt_gen, corrected reads)


def _run_graph(s, a, prep_dir):
    def run(ctx, wd):
        summ = G.run(a["ctg_gfa"], a["utg_gfa"], wd, peak_hom=_peak_hom(prep_dir), kmer_peak=s["kmer_peak"])
        if summ["contig_reads_found_in_utg"] < 0.9:
            raise StageError(f"only {summ['contig_reads_found_in_utg']:.1%} of contig reads found in the unitig "
                             "graph: the two GFA files do not belong together")
        return {"outputs": ["summary.json", "contig_utg_path.tsv", "path_transitions.tsv", "utg_segments.tsv"],
                "validation": {"contig_reads_found_in_utg": summ["contig_reads_found_in_utg"],
                               "rd_reference": summ["rd_reference"], "rd_reference_source": summ["rd_reference_source"]}}
    return run


def _run_ref_index(cfg, r):
    def run(ctx, wd):
        name = "ref.fa.gz" if r["fasta"].endswith(".gz") else "ref.fa"
        dst = link_fasta(r["fasta"], wd, name)
        if not os.path.exists(dst + ".fai"):
            ctx.sh(f"{tool(cfg, 'samtools')} faidx {name}", wd, pin=False)
        # N-gap fraction per coverage bin
        binsz = ctx.sci["ref_align"]["coverage_bin"]
        rx = re.compile(r["chromosomes_regex"])
        rows, opener = [], gzip.open if name.endswith(".gz") else open
        cur, pos, nb = None, 0, {}
        with opener(dst, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    cur, pos = line[1:].split()[0], 0
                    continue
                line = line.rstrip()
                if cur and rx.match(cur):
                    n = line.count("N") + line.count("n")
                    if n:
                        for i, ch in enumerate(line):
                            if ch in "Nn":
                                k = (cur, (pos + i) // binsz * binsz)
                                nb[k] = nb.get(k, 0) + 1
                pos += len(line)
        lengths = dict(read_fai(dst + ".fai"))
        with open(os.path.join(wd, "ngap_bins.tsv"), "w") as fh:
            fh.write("chrom\tbin_start\tn_frac\n")
            for (c, b), n in sorted(nb.items()):
                fh.write(f"{c}\t{b}\t{n / min(binsz, lengths[c] - b):.4f}\n")
        preset = ctx.sci["ref_align"]["preset"]
        ctx.sh(f"{tool(cfg, 'minimap2')} -x {preset} -t {ctx.res['threads']} -d ref.{preset}.mmi {name} 2> minimap2_index.log", wd)
        return {"outputs": [name, "ngap_bins.tsv", f"ref.{preset}.mmi"], "validation": {"sequences": len(lengths)}}
    return run


def _run_ref_align(cfg, r, asm_fa, sname, aname):
    def run(ctx, wd):
        preset = ctx.sci["ref_align"]["preset"]
        mmi = os.path.join(ref_dir(cfg, r["name"]), f"ref.{preset}.mmi")
        ctx.sh(f"{tool(cfg, 'minimap2')} -x {preset} -t {ctx.res['threads']} --secondary=no {mmi} {asm_fa} "
               f"2> minimap2.log | gzip -c > aln.paf.gz", wd)
        n = sum(1 for _ in gzip.open(os.path.join(wd, "aln.paf.gz"), "rt"))
        if n == 0:
            raise StageError("no alignments")
        return {"outputs": ["aln.paf.gz", "minimap2.log"], "validation": {"paf_records": n}}
    return run


def _ref_lengths_and_gaps(cfg, r):
    d = ref_dir(cfg, r["name"])
    name = "ref.fa.gz" if r["fasta"].endswith(".gz") else "ref.fa"
    lengths = dict(read_fai(os.path.join(d, name + ".fai")))
    gaps = {}
    for line in readtext(os.path.join(d, "ngap_bins.tsv")).splitlines()[1:]:
        c, b, f = line.split("\t")
        gaps[(c, int(b))] = float(f)
    return lengths, gaps


def _run_ref_analysis(cfg, r, align_dir, prep_dir, sname, a):
    def run(ctx, wd):
        lengths, gaps = _ref_lengths_and_gaps(cfg, r)
        contigs = dict(_contigs(prep_dir))
        summ, by_q, assign = refalign.run(os.path.join(align_dir, "aln.paf.gz"), contigs, lengths, gaps,
                                          ctx.sci, r["chromosomes_regex"], wd)
        plots.dotplot_genome(by_q, assign, lengths, r["chromosomes_regex"], ctx.sci["ref_align"]["plot_min_len"],
                             f"{sname} {a['label']} vs {r['label']}", r["label"], f"{sname} {a['label']}",
                             os.path.join(wd, "dotplot_genome"))
        return {"outputs": ["contig_assignment.tsv", "junctions.tsv", "synteny_segments.tsv",
                            "reference_coverage.tsv", "summary.json", "dotplot_genome.svg", "dotplot_genome.pdf"],
                "validation": {"aligned_bp": summ["aligned_bp"], "junctions": summ["junctions"]}}
    return run


def _run_evidence(cfg, a, regions_dir, read_dir, ref_dirs, graph_dir, label):
    def run(ctx, wd):
        summ = evidence.run(regions_dir, read_dir, ref_dirs, graph_dir, a["utg_gfa"], ctx.sci, wd, label)
        return {"outputs": ["candidates.tsv", "summary.json"], "validation": {"candidates": summ["candidates"]}}
    return run


def _run_pair_align(cfg, sname, pr, A, B, extra):
    def run(ctx, wd):
        p = ctx.sci["pair_align"]
        a_fa = os.path.join(asm_dir(cfg, sname, pr["a"], "prepare"), fasta_link_name(A["fasta"]))
        b_fa = os.path.join(asm_dir(cfg, sname, pr["b"], "prepare"), fasta_link_name(B["fasta"]))
        mm = tool(cfg, "minimap2")
        t = ctx.res["threads"]
        ctx.sh(f"{mm} -c -x {p['preset']} -t {t} --secondary=no {b_fa} {a_fa} 2> minimap2.log | gzip -c > a_vs_b.paf.gz", wd)
        outs = ["a_vs_b.paf.gz", "minimap2.log"]
        if extra:
            for x in extra:
                if not os.path.exists(x):
                    raise StageError(f"extra target missing: {x}")
            ctx.sh(f"{mm} -c -x {p['preset']} -t {t} --secondary=no <(zcat -f {' '.join(extra)}) {a_fa} "
                   f"2> minimap2_extra.log | gzip -c > a_vs_b_extra.paf.gz", wd)
            outs += ["a_vs_b_extra.paf.gz", "minimap2_extra.log"]
        return {"outputs": outs, "validation": {"extra_targets": extra}}
    return run


def _run_pair_analysis(cfg, sname, pr, A, B, align_dir, have_graph):
    def run(ctx, wd):
        a_len = dict(_contigs(asm_dir(cfg, sname, pr["a"], "prepare")))
        b_len = dict(_contigs(asm_dir(cfg, sname, pr["b"], "prepare")))
        ac = pairalign.WindowClasses(os.path.join(asm_dir(cfg, sname, pr["a"], "regions"), "window_classes.tsv.gz"))
        bc = pairalign.WindowClasses(os.path.join(asm_dir(cfg, sname, pr["b"], "regions"), "window_classes.tsv.gz"))
        extra = os.path.join(align_dir, "a_vs_b_extra.paf.gz")
        summ = pairalign.run(os.path.join(align_dir, "a_vs_b.paf.gz"), a_len, b_len, ctx.sci, wd,
                             extra_paf=extra if os.path.exists(extra) else None, a_classes=ac, b_classes=bc)
        plots.pair_class_bar(summ, A["label"], B["label"], os.path.join(wd, "sequence_classes"))
        ih = [tuple(map(float, l.split("\t"))) for l in readtext(os.path.join(wd, "redundant_identity_hist.tsv")).splitlines()[1:]]
        if ih:
            plots.identity_hist(ih, os.path.join(wd, "redundant_identity"), f"{sname}: extra {A['label']} copies")
        outs = ["summary.json", "redundant_copies.tsv", "a_unaligned.tsv", "a_contig_classes.tsv",
                "sequence_classes.svg"]
        if have_graph:
            with open(os.path.join(asm_dir(cfg, sname, pr["a"], "graph"), "summary.json")) as fh:
                ref_rd = json.load(fh)["rd_reference"]
            r = G.retention(A["ctg_gfa"], B["ctg_gfa"], A["utg_gfa"], wd, ref_rd=ref_rd)
            outs += ["read_retention_summary.json", "read_retention_by_contig.tsv", "dropped_intervals.tsv"]
            summ["read_retention"] = r
        return {"outputs": outs, "validation": {"class_bp": summ["class_bp"]}}
    return run
