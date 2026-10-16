"""hifiasm assembly-graph analysis from GFA files without sequences
(*.noseq.gfa are sufficient: S lines carry LN and rd, A lines the reads).

- Unitig-graph topology: segments, links, tips, branching ends, components,
  simple bubbles (two parallel segments between the same pair of ends), and
  the length-weighted distribution of the rd (read depth) tag.
- Contig -> unitig paths: hifiasm A lines list the reads that make up each
  contig and each unitig. Joining them on read name maps contig coordinates
  to unitig segments, which is how candidate regions are located in the graph.
  Transitions between consecutive segments are labelled linked / unlinked
  (no GFA link between them) and branching (the exit end has >= 2 links).
- Local neighbourhoods of candidate regions: subgraph GFA + SVG/PDF drawing.

The rd tag is hifiasm's read depth of a segment. Depth classes are relative
to a homozygous (single-copy, both haplotypes) reference depth: hifiasm's own
peak_hom estimate from its log, else the sample's k-mer peak, else the
length-weighted rd mode. The mode is NOT safe on its own: in a sample that is
heterozygous over much of the genome the mode is the haploid (one-allele)
depth. In inbred lines rd modes match peak_hom closely.
"""
import json
import os

import networkx as nx


def parse_gfa(path, want_a=True):
    seg, links, alines = {}, [], {}
    with open(path) as fh:
        for line in fh:
            t = line[:2]
            if t == "S\t":
                p = line.rstrip("\n").split("\t")
                tags = {x[:2]: x[5:] for x in p[3:]}
                ln = int(tags["LN"]) if "LN" in tags else (len(p[2]) if p[2] != "*" else 0)
                seg[p[1]] = {"len": ln, "rd": int(tags["rd"]) if "rd" in tags else None}
            elif t == "L\t":
                p = line.split("\t")
                links.append((p[1], p[2], p[3], p[4], p[5]))
            elif t == "A\t" and want_a:
                p = line.split("\t")
                rlen = int(p[6]) - int(p[5])
                alines.setdefault(p[1], []).append((int(p[2]), rlen, p[4]))
    for v in alines.values():
        v.sort()
    return {"segments": seg, "links": links, "alines": alines}


def _flip(o):
    return "-" if o == "+" else "+"


def oriented_adjacency(g):
    """out[(u,o)] = set of (v,o2) reachable by one link, using both link directions."""
    out = {}
    for u, uo, v, vo, _ in g["links"]:
        out.setdefault((u, uo), set()).add((v, vo))
        out.setdefault((v, _flip(vo)), set()).add((u, _flip(uo)))
    return out


def topology(g):
    seg, out = g["segments"], oriented_adjacency(g)
    tips = branching = 0
    for s in seg:
        for o in "+-":
            n = len(out.get((s, o), ()))
            tips += n == 0
            branching += n >= 2
    und = nx.Graph()
    und.add_nodes_from(seg)
    und.add_edges_from((u, v) for u, _, v, _, _ in g["links"])
    comps = list(nx.connected_components(und))
    singletons = sum(1 for c in comps if len(c) == 1)
    # simple bubbles: (s,o) -> {x, y} -> (t,o') with x, y having that single in/out
    bubbles, bubble_bp, seen = 0, 0, set()
    for start, nbrs in out.items():
        if len(nbrs) < 2:
            continue
        exits = {}
        for x in nbrs:
            succ = out.get(x, set())
            pred = out.get((x[0], _flip(x[1])), set())
            if len(succ) == 1 and len(pred) == 1:
                exits.setdefault(next(iter(succ)), []).append(x)
        for end, xs in exits.items():
            if len(xs) >= 2:
                key = frozenset([start, (end[0], _flip(end[1]))])
                if key in seen:
                    continue
                seen.add(key)
                bubbles += 1
                bubble_bp += sum(seg[x[0]]["len"] for x in xs)
    lengths = sorted((s["len"] for s in seg.values()), reverse=True)
    total = sum(lengths)
    acc, n50 = 0, 0
    for l in lengths:
        acc += l
        if acc >= total / 2:
            n50 = l
            break
    undirected_links = len({frozenset([(u, uo), (v, vo)]) for u, uo, v, vo, _ in g["links"]})
    return {"segments": len(seg), "links": undirected_links, "total_bp": total, "segment_n50": n50,
            "tips": tips, "branching_ends": branching, "components": len(comps),
            "singleton_components": singletons, "simple_bubbles": bubbles, "bubble_branch_bp": bubble_bp,
            "rd_mode": rd_mode(seg)}


def rd_mode(seg, min_len=100000):
    """Length-weighted mode of rd over long segments (all segments if none is long)."""
    for floor in (min_len, 0):  # small graphs: fall back to all segments
        w = {}
        for s in seg.values():
            if s["rd"] is not None and s["len"] >= floor:
                w[s["rd"]] = w.get(s["rd"], 0) + s["len"]
        if w:
            return max(w, key=w.get)
    return None


def rd_class(rd, mode):
    if rd is None or not mode:
        return "na"
    r = rd / mode
    if r < 0.3:
        return "low"
    if r < 0.7:
        return "half"
    if r < 1.4:
        return "normal"
    if r < 2.5:
        return "dup"
    return "multi"


def rd_class_bp(seg, mode):
    out = {}
    for s in seg.values():
        c = rd_class(s["rd"], mode)
        out[c] = out.get(c, 0) + s["len"]
    return out


def read_index(g):
    """read name -> (segment, offset) for a unitig graph's A lines."""
    idx = {}
    for s, lst in g["alines"].items():
        for off, _, r in lst:
            idx[r] = (s, off)
    return idx


def contig_paths(ctg, utg, uidx=None):
    """Runs of consecutive contig reads that belong to the same unitig.
    Returns rows (contig, start, end, utg, n_reads) and the fraction of contig
    reads found in the unitig graph."""
    uidx = uidx or read_index(utg)
    rows, found, total = [], 0, 0
    for c, lst in ctg["alines"].items():
        clen = ctg["segments"][c]["len"]
        cur = None
        for i, (off, rlen, r) in enumerate(lst):
            total += 1
            u = uidx.get(r, (None,))[0]
            found += u is not None
            u = u or "NA"
            if cur and cur[3] == u:
                cur[2] = max(cur[2], off + rlen)
                cur[4] += 1
            else:
                if cur:
                    cur[2] = max(cur[1] + 1, min(cur[2], off))
                    rows.append(cur)
                cur = [c, off, off + rlen, u, 1]
        if cur:
            cur[2] = min(max(cur[2], cur[1] + 1), clen)
            rows.append(cur)
    return rows, (found / total if total else 0.0)


def annotate_transitions(rows, utg):
    out = oriented_adjacency(utg)
    linked = {}
    for (u, _), nb in out.items():
        for v, _ in nb:
            linked.setdefault(u, set()).add(v)
    ann = []
    for a, b in zip(rows, rows[1:]):
        if a[0] != b[0] or a[3] == "NA" or b[3] == "NA":
            continue
        is_linked = b[3] in linked.get(a[3], ())
        deg = max(len(out.get((a[3], "+"), ())), len(out.get((a[3], "-"), ())))
        ann.append({"contig": a[0], "pos": b[1], "from": a[3], "to": b[3], "linked": int(is_linked),
                    "exit_degree": deg})
    return ann


def segments_in(rows_by_contig, contig, s, e):
    return sorted({r[3] for r in rows_by_contig.get(contig, []) if r[1] < e and r[2] > s and r[3] != "NA"})


def neighborhood(utg, seeds, radius, max_nodes):
    und = nx.Graph()
    for u, _, v, _, _ in utg["links"]:
        und.add_edge(u, v)
    nodes = set(s for s in seeds if s in utg["segments"])
    frontier = set(nodes)
    for _ in range(radius):
        nxt = set()
        for n in frontier:
            if n in und:
                nxt.update(und.neighbors(n))
        nxt -= nodes
        if len(nodes) + len(nxt) > max_nodes:
            nxt = set(sorted(nxt)[: max(0, max_nodes - len(nodes))])
        nodes |= nxt
        frontier = nxt
        if not frontier:
            break
    return nodes


def write_subgraph(utg, nodes, path):
    with open(path, "w") as fh:
        fh.write("H\tVN:Z:1.0\n")
        for n in sorted(nodes):
            s = utg["segments"][n]
            rd = f"\trd:i:{s['rd']}" if s["rd"] is not None else ""
            fh.write(f"S\t{n}\t*\tLN:i:{s['len']}{rd}\n")
        for u, uo, v, vo, ov in utg["links"]:
            if u in nodes and v in nodes:
                fh.write(f"L\t{u}\t{uo}\t{v}\t{vo}\t{ov.strip()}\n")


CLASS_COLORS = {"low": "#c7c7c7", "half": "#4c78a8", "normal": "#8c8c8c", "dup": "#f58518",
                "multi": "#e45756", "na": "#ffffff"}


def render(utg, nodes, seeds, path_nodes, mode, title, out_base):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    g = nx.Graph()
    g.add_nodes_from(nodes)
    for u, _, v, _, _ in utg["links"]:
        if u in nodes and v in nodes and u != v:
            g.add_edge(u, v)
    pos = nx.kamada_kawai_layout(g) if 1 < len(g) <= 120 else nx.spring_layout(g, seed=1)
    fig, ax = plt.subplots(figsize=(8, 7))
    nx.draw_networkx_edges(g, pos, ax=ax, edge_color="#999999", width=1.0)
    for n in g.nodes:
        s = utg["segments"][n]
        size = 60 + 90 * max(0.0, min(6.0, __import__("math").log10(max(s["len"], 10)) - 3))
        ax.scatter(*pos[n], s=size, c=CLASS_COLORS[rd_class(s["rd"], mode)],
                   edgecolors="black" if n in seeds else ("#333333" if n in path_nodes else "none"),
                   linewidths=2.5 if n in seeds else (1.0 if n in path_nodes else 0), zorder=3)
        ax.annotate(f"{n.replace('utg', '').lstrip('0').rstrip('lc')}\n{s['len'] / 1000:.0f} kb",
                    pos[n], fontsize=6, ha="center", va="bottom", xytext=(0, 6), textcoords="offset points")
    handles = [Line2D([], [], marker="o", ls="", markerfacecolor=c, markeredgecolor="#555555", label=k)
               for k, c in CLASS_COLORS.items() if k != "na"]
    handles.append(Line2D([], [], marker="o", ls="", markerfacecolor="white", markeredgecolor="black",
                          markeredgewidth=2.5, label="candidate segment"))
    ax.legend(handles=handles, title=f"rd / homozygous depth ({mode})", fontsize=7, title_fontsize=7, loc="best")
    ax.set_title(title, fontsize=9)
    ax.axis("off")
    fig.tight_layout()
    for ext in ("svg", "pdf"):
        fig.savefig(f"{out_base}.{ext}")
    plt.close(fig)


def reference_rd(seg, peak_hom=None, kmer_peak=None):
    """(value, source) of the homozygous reference depth for rd classes."""
    if peak_hom:
        return int(peak_hom), "hifiasm peak_hom (assembler log)"
    if kmer_peak:
        return int(kmer_peak), "configured k-mer peak"
    return rd_mode(seg), "rd mode (may be the haploid depth in a heterozygous sample)"


def rd_hist(seg, min_len=100000, step=2):
    """Length-weighted histogram of rd over long segments, in bins of `step`."""
    h = {}
    for s in seg.values():
        if s["rd"] is not None and s["len"] >= min_len:
            b = s["rd"] // step * step
            h[b] = h.get(b, 0) + s["len"]
    return dict(sorted(h.items()))


def run(ctg_gfa, utg_gfa, out_dir, peak_hom=None, kmer_peak=None):
    """Graph metrics and contig->unitig paths for one assembly."""
    ctg = parse_gfa(ctg_gfa)
    utg = parse_gfa(utg_gfa)
    topo_u = topology(utg)
    topo_c = topology(ctg)
    rows, found = contig_paths(ctg, utg)
    trans = annotate_transitions(rows, utg)
    mode, mode_source = reference_rd(utg["segments"], peak_hom, kmer_peak)
    with open(os.path.join(out_dir, "contig_utg_path.tsv"), "w") as fh:
        fh.write("contig\tstart\tend\tutg\tn_reads\tutg_len\tutg_rd\trd_class\n")
        for r in rows:
            s = utg["segments"].get(r[3], {"len": 0, "rd": None})
            fh.write("\t".join(map(str, r)) + f"\t{s['len']}\t{s['rd']}\t{rd_class(s['rd'], mode)}\n")
    with open(os.path.join(out_dir, "path_transitions.tsv"), "w") as fh:
        fh.write("contig\tpos\tfrom\tto\tlinked\texit_degree\n")
        for t in trans:
            fh.write("\t".join(str(t[k]) for k in ("contig", "pos", "from", "to", "linked", "exit_degree")) + "\n")
    with open(os.path.join(out_dir, "utg_segments.tsv"), "w") as fh:
        fh.write("segment\tlength\trd\trd_class\n")
        for n, s in utg["segments"].items():
            fh.write(f"{n}\t{s['len']}\t{s['rd']}\t{rd_class(s['rd'], mode)}\n")
    # bp of contig sequence by unitig rd class
    ctg_rd = {}
    for r in rows:
        s = utg["segments"].get(r[3])
        c = rd_class(s["rd"], mode) if s else "na"
        ctg_rd[c] = ctg_rd.get(c, 0) + (r[2] - r[1])
    summary = {"utg": topo_u, "ctg": topo_c, "contig_reads_found_in_utg": round(found, 4),
               "rd_reference": mode, "rd_reference_source": mode_source,
               "utg_rd_hist_bp_long_segments": rd_hist(utg["segments"]),
               "utg_rd_class_bp": rd_class_bp(utg["segments"], mode), "contig_bp_by_utg_rd_class": ctg_rd,
               "path_transitions": len(trans), "unlinked_transitions": sum(1 for t in trans if not t["linked"]),
               "branching_transitions": sum(1 for t in trans if t["exit_degree"] >= 2),
               "ctg_links": len(ctg["links"])}
    with open(os.path.join(out_dir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary


def load_paths(path):
    by = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            by.setdefault(p[0], []).append([p[0], int(p[1]), int(p[2]), p[3], int(p[4])])
    return by


def retention(ctg_a_gfa, ctg_b_gfa, utg_a_gfa, out_dir, ref_rd=None):
    """Which reads of A's contigs are also in B's contigs (same reads, same
    assembler run except purging): per-A-contig retention and the unitig rd
    class of the dropped sequence."""
    a = parse_gfa(ctg_a_gfa)
    b = parse_gfa(ctg_b_gfa)
    utg = parse_gfa(utg_a_gfa)
    mode = ref_rd or rd_mode(utg["segments"])
    uidx = read_index(utg)
    b_reads = {r for lst in b["alines"].values() for _, _, r in lst}
    a_reads = {r for lst in a["alines"].values() for _, _, r in lst}
    dropped_by_class, kept_by_class, rows, dropped_iv = {}, {}, [], []
    for c, lst in a["alines"].items():
        kept = 0
        for i, (off, rlen, r) in enumerate(lst):
            nxt = lst[i + 1][0] if i + 1 < len(lst) else a["segments"][c]["len"]
            span = max(0, min(nxt, a["segments"][c]["len"]) - off)
            u = uidx.get(r, (None,))[0]
            cls = rd_class(utg["segments"][u]["rd"], mode) if u else "na"
            if r in b_reads:
                kept += span
                kept_by_class[cls] = kept_by_class.get(cls, 0) + span
            else:
                dropped_by_class[cls] = dropped_by_class.get(cls, 0) + span
                last = dropped_iv[-1] if dropped_iv else None
                if last and last[0] == c and last[2] == off and last[3] == cls and last[4] == (u or "NA"):
                    last[2] = off + span
                else:
                    dropped_iv.append([c, off, off + span, cls, u or "NA"])
        rows.append((c, a["segments"][c]["len"], kept))
    with open(os.path.join(out_dir, "dropped_intervals.tsv"), "w") as fh:
        fh.write("a_contig\tstart\tend\tutg_rd_class\tutg\n")
        for iv in dropped_iv:
            if iv[2] > iv[1]:
                fh.write("\t".join(map(str, iv)) + "\n")
    with open(os.path.join(out_dir, "read_retention_by_contig.tsv"), "w") as fh:
        fh.write("a_contig\tlength\tbp_with_reads_in_b\tretained_frac\n")
        for c, l, k in sorted(rows, key=lambda x: -x[1]):
            fh.write(f"{c}\t{l}\t{k}\t{k / l if l else 0:.4f}\n")
    summary = {"a_contig_reads": len(a_reads), "b_contig_reads": len(b_reads),
               "b_reads_also_in_a": len(b_reads & a_reads) / len(b_reads) if b_reads else 0,
               "a_reads_also_in_b": len(b_reads & a_reads) / len(a_reads) if a_reads else 0,
               "rd_reference": mode, "kept_bp_by_utg_rd_class": kept_by_class,
               "dropped_bp_by_utg_rd_class": dropped_by_class}
    with open(os.path.join(out_dir, "read_retention_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    return summary
