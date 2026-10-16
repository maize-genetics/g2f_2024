"""Run configuration for assembly validation: loading, defaults, validation.

A run configuration (YAML) names the output directory, tools, resources, the
scientific parameter file and the samples/assemblies to validate. Relative
paths are resolved against the pipeline/ directory, like the rest of the
pipeline. See config/asmqc.example.yaml.
"""
import copy
import hashlib
import json
import os
import re

import yaml

PIPELINE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")

DEFAULT_TOOLS = {
    "minimap2": {"cmd": "minimap2", "version_cmd": "minimap2 --version"},
    "samtools": {"cmd": "samtools", "version_cmd": "samtools --version | head -1"},
    "seqkit": {"cmd": "seqkit", "version_cmd": "seqkit version"},
}

DEFAULT_RESOURCES = {
    "threads": 16,          # threads for minimap2 and parallel Python workers
    "cpus": "",             # optional taskset CPU list, e.g. "64-95"
    "nice": 10,
    "sort_threads": 8,
    "sort_mem_per_thread": "2G",
    "reserve_disk_gb": 100,  # free space that must remain after a stage's estimate
    "min_free_ram_gb": 64,
    "alignment_format": "cram",  # cram (needs the assembly FASTA to read) or bam
}


class ConfigError(Exception):
    pass


def resolve(path, base=PIPELINE_DIR):
    if path is None or path == "":
        return None
    path = os.path.expanduser(str(path))
    return path if os.path.isabs(path) else os.path.normpath(os.path.join(base, path))


def _hash(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True).encode()).hexdigest()[:12]


def _need(d, key, where):
    if key not in d or d[key] in (None, ""):
        raise ConfigError(f"{where}: missing '{key}'")
    return d[key]


def _check_name(name, where):
    if not NAME_RE.match(str(name)):
        raise ConfigError(f"{where}: invalid name {name!r} (letters, digits, . _ - only)")


def load(path):
    path = resolve(path, os.getcwd())
    if not os.path.exists(path):
        raise ConfigError(f"configuration file not found: {path}")
    with open(path) as fh:
        raw = yaml.safe_load(fh) or {}
    cfg = {"config_path": path}

    cfg["profile"] = str(_need(raw, "profile", "config"))
    _check_name(cfg["profile"], "profile")
    cfg["outdir"] = resolve(_need(raw, "outdir", "config"))
    cfg["root"] = os.path.join(cfg["outdir"], cfg["profile"])

    science_path = resolve(raw.get("science", "config/asmqc_science.yaml"))
    with open(science_path) as fh:
        science = yaml.safe_load(fh)
    for k, v in (raw.get("science_overrides") or {}).items():
        if k not in science or not isinstance(v, dict):
            raise ConfigError(f"science_overrides: unknown section {k!r}")
        science[k].update(v)
    cfg["science"] = science
    cfg["science_path"] = science_path
    cfg["science_hash"] = _hash(science)

    tools = copy.deepcopy(DEFAULT_TOOLS)
    if raw.get("tools_file"):
        with open(resolve(raw["tools_file"])) as fh:
            tools.update((yaml.safe_load(fh) or {}).get("tools", {}))
    tools.update(raw.get("tools") or {})
    cfg["tools"] = tools

    res = dict(DEFAULT_RESOURCES)
    res.update(raw.get("resources") or {})
    if res["alignment_format"] not in ("cram", "bam"):
        raise ConfigError("resources.alignment_format must be cram or bam")
    cfg["resources"] = res
    cfg["tmpdir"] = resolve(raw.get("tmpdir") or os.path.join(cfg["outdir"], "tmp"))

    refs = {}
    for name, r in (raw.get("references") or {}).items():
        _check_name(name, "references")
        fasta = resolve(_need(r, "fasta", f"references.{name}"))
        refs[name] = {
            "name": name,
            "fasta": fasta,
            "label": r.get("label", name),
            "chromosomes_regex": r.get("chromosomes_regex", r"^chr\d+$"),
        }
    cfg["references"] = refs

    samples = {}
    for sname, s in (_need(raw, "samples", "config") or {}).items():
        _check_name(sname, "samples")
        where = f"samples.{sname}"
        reads = [resolve(p) for p in _need(s, "hifi_reads", where)]
        asms = {}
        for aname, a in _need(s, "assemblies", where).items():
            _check_name(aname, f"{where}.assemblies")
            aw = f"{where}.assemblies.{aname}"
            asms[aname] = {
                "name": aname,
                "label": a.get("label", aname),
                "method_id": a.get("method_id", ""),
                "fasta": resolve(_need(a, "fasta", aw)),
                "ctg_gfa": resolve(a.get("ctg_gfa")),
                "utg_gfa": resolve(a.get("utg_gfa")),
                "assembler_log": resolve(a.get("assembler_log")),
                "extra_targets": [resolve(p) for p in a.get("extra_targets", [])],
                "existing_metrics": a.get("existing_metrics") or {},
            }
        samples[sname] = {
            "name": sname,
            "hifi_reads": reads,
            "hifi_read_count": s.get("hifi_read_count"),
            "kmer_peak": s.get("kmer_peak"),
            "expected_depth": s.get("expected_depth"),
            "assemblies": asms,
            "pairs": [],
        }
        for p in s.get("pairs", []) or []:
            a, b = _need(p, "a", f"{where}.pairs"), _need(p, "b", f"{where}.pairs")
            for x in (a, b):
                if x not in asms:
                    raise ConfigError(f"{where}.pairs: unknown assembly {x!r}")
            samples[sname]["pairs"].append({"a": a, "b": b, "use_b_extra_targets": bool(p.get("use_b_extra_targets", True))})
    cfg["samples"] = samples
    cfg["existing_metrics_tsv"] = resolve(raw.get("existing_metrics_tsv"))
    return cfg


def tool(cfg, name):
    try:
        return cfg["tools"][name]["cmd"]
    except KeyError:
        raise ConfigError(f"tool {name!r} is not configured")
