"""Build per-sample summaries from validated JSON records only."""
import csv
import json

COLUMNS = [
    "sample", "method_id", "science_version", "status",
    "hifi_reads", "hifi_bases", "hifi_approx_cov",
    "ont_min_read_len", "ont_raw_reads", "ont_raw_bases", "ont_retained_reads",
    "ont_retained_bases", "ont_retained_read_n50", "ont_retained_approx_cov",
    "kmer_main_peak", "kmer_screen_flagged",
    "contigs", "total_bp", "contig_n50", "contig_l50", "max_contig", "gc_pct", "gap_bp",
    "alt_contigs", "alt_total_bp",
    "busco_lineage", "busco_dataset", "busco_C", "busco_S", "busco_D", "busco_F", "busco_M",
    "busco_n", "busco_internal_stop_pct",
    "merqury_qv", "merqury_kmer_completeness_pct",
    "hifiasm_version", "hifiasm_peak_rss_gb",
]


def _load(p):
    with open(p) as fh:
        return json.load(fh)


def build(sample, method_id, science_version, hifi, ont, kmer, asm_p, asm_a, busco, merqury,
          hifiasm_log):
    for name, rec in (("hifi", hifi), ("assembly", asm_p), ("alt assembly", asm_a),
                      ("kmer", kmer), *(("busco", b) for b in busco),
                      *((("merqury", merqury),) if merqury else ()),
                      *((("ont", ont),) if ont else ())):
        if not rec.get("ok"):
            raise ValueError(f"{sample}/{method_id}: {name} record is not validated")
    rows = []
    for b in busco or [None]:
        row = {
            "sample": sample, "method_id": method_id, "science_version": science_version,
            "status": "validated",
            "hifi_reads": hifi["total_reads"], "hifi_bases": hifi["total_bases"],
            "hifi_approx_cov": hifi["approx_coverage"],
            "kmer_main_peak": kmer.get("main_peak_coverage"),
            "kmer_screen_flagged": kmer.get("flagged"),
            "contigs": asm_p["n_contigs"], "total_bp": asm_p["total_bp"],
            "contig_n50": asm_p["n50"], "contig_l50": asm_p["l50"], "max_contig": asm_p["max_len"],
            "gc_pct": asm_p["gc_pct"], "gap_bp": asm_p["gap_bp"],
            "alt_contigs": asm_a["n_contigs"], "alt_total_bp": asm_a["total_bp"],
            "hifiasm_version": hifiasm_log.get("version"),
            "hifiasm_peak_rss_gb": hifiasm_log.get("peak_rss_gb"),
        }
        if ont:
            row.update(ont_min_read_len=ont["min_read_len"], ont_raw_reads=ont["raw_reads"],
                       ont_raw_bases=ont["raw_bases"], ont_retained_reads=ont["retained_reads"],
                       ont_retained_bases=ont["retained_bases"],
                       ont_retained_read_n50=ont["retained_read_n50"],
                       ont_retained_approx_cov=ont["retained_approx_coverage"])
        if b:
            ds = b.get("lineage_dataset_cfg") or {}
            row.update(busco_lineage=b["lineage"],
                       busco_dataset=f"OrthoDB {ds.get('OrthoDB_version', '?')}, created "
                                     f"{ds.get('creation_date', '?')}",
                       busco_C=b["complete_pct"], busco_S=b["single_pct"],
                       busco_D=b["duplicated_pct"], busco_F=b["fragmented_pct"],
                       busco_M=b["missing_pct"], busco_n=b["n_markers"],
                       busco_internal_stop_pct=b.get("internal_stop_codon_percent"))
        if merqury:
            row.update(merqury_qv=merqury.get("qv"),
                       merqury_kmer_completeness_pct=merqury.get("kmer_completeness_pct"))
        rows.append({c: row.get(c, "") for c in COLUMNS})
    return rows


def write_tsv(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))
