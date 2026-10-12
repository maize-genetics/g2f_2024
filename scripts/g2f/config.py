"""Load and validate the layered workflow configuration."""
import hashlib
import json
import os
import re

import yaml

SAMPLE_ID_RE = r"^[A-Za-z0-9][A-Za-z0-9_.-]*$"


class ConfigError(Exception):
    pass


def _load_yaml(path, required=True):
    if not os.path.exists(path):
        if required:
            raise ConfigError(f"Configuration file not found: {path}")
        return {}
    with open(path) as fh:
        data = yaml.safe_load(fh)
    return data or {}


def _resolve(base, path):
    if path is None:
        return None
    path = os.path.expanduser(str(path))
    return os.path.normpath(path if os.path.isabs(path) else os.path.join(base, path))


def _inside(child, parent):
    child, parent = os.path.realpath(child), os.path.realpath(parent)
    return child == parent or child.startswith(parent.rstrip(os.sep) + os.sep)


def _as_list(value):
    if not value:
        return []
    if isinstance(value, str):
        return [v.strip() for v in value.split(",") if v.strip()]
    return [str(v) for v in value]


def load_methods(methods_dir):
    methods = {}
    for name in sorted(os.listdir(methods_dir)):
        if not name.endswith((".yaml", ".yml")):
            continue
        m = _load_yaml(os.path.join(methods_dir, name))
        mid = m.get("method_id")
        if not mid:
            raise ConfigError(f"{name}: missing method_id")
        if os.path.splitext(name)[0] != mid:
            raise ConfigError(f"{name}: file name must equal method_id ({mid})")
        for key in ("assembler", "expected_version", "read_types", "args"):
            if key not in m:
                raise ConfigError(f"{name}: missing required key '{key}'")
        if "ont" in m["read_types"]:
            ont = m.get("ont") or {}
            if not isinstance(ont.get("min_read_len"), int) or ont["min_read_len"] < 0:
                raise ConfigError(f"{name}: ont.min_read_len must be a non-negative integer")
        m["_file"] = os.path.join(methods_dir, name)
        methods[mid] = m
    return methods


def load_config(config, pipeline_dir):
    """Expand the top-level Snakemake config dict into a validated settings dict."""
    base = pipeline_dir
    cfg = {
        "pipeline_dir": base,
        "data_dir": _resolve(base, config.get("data_dir", "../data")),
        "results_dir": _resolve(base, config.get("results_dir", "../results")),
        "legacy_dirs": [_resolve(base, p) for p in config.get("legacy_dirs", [])],
        "only_samples": _as_list(config.get("only_samples")),
    }
    cfg["samples"] = (_load_yaml(_resolve(base, config.get("samples", "config/samples.yaml")))
                      .get("samples") or {})
    cfg["science"] = _load_yaml(_resolve(base, config.get("science", "config/science.yaml")))
    cfg["methods"] = load_methods(_resolve(base, config.get("methods_dir", "config/methods")))
    cfg["resources"] = _load_yaml(_resolve(base, config.get("resources", "config/resources/default.yaml")))
    cfg["local"] = _load_yaml(_resolve(base, config.get("local", "config/local/local.yaml")),
                              required=False)
    tools_cfg = _load_yaml(_resolve(base, config.get("tools", "config/tools/biohpc.yaml")))
    containers = _resolve(base, cfg["local"].get("containers_dir")
                          or tools_cfg.get("containers_dir") or ".")
    cfg["tools"] = {
        name: {k: str(v).replace("{containers}", str(containers)) for k, v in spec.items()}
        for name, spec in (tools_cfg.get("tools") or {}).items()
    }
    cfg["fingerprint"] = (_load_yaml(_resolve(base, config.get("fingerprint", "config/fingerprint.yaml")),
                                     required=False).get("fingerprint") or {})
    cfg["tmpdir"] = _resolve(base, config.get("tmpdir") or cfg["resources"].get("tmpdir", "../tmp"))
    cfg["external_assemblies"] = load_external_assemblies(
        _resolve(base, config.get("external_assemblies", "config/external_assemblies.yaml")), base)
    validate(cfg)
    return cfg


EXTERNAL_KEYS = ("p_ctg_gfa", "p_ctg_fa", "hifiasm_log")
EXTERNAL_QC = ("stats", "busco", "merqury")


def load_external_assemblies(path, base):
    """Read-only historical assemblies registered for QC benchmarking."""
    data = _load_yaml(path, required=False)
    if not data:
        return {}
    root = _resolve(base, data.get("root", "."))
    out = {}
    for sample, entries in sorted((data.get("assemblies") or {}).items()):
        for aid, e in sorted((entries or {}).items()):
            if not re.match(r"^[A-Za-z0-9_]+$", aid):
                raise ConfigError(f"external assembly id {aid!r}: use letters, digits, '_'")
            for k in EXTERNAL_KEYS:
                if not e.get(k):
                    raise ConfigError(f"external assembly {sample}/{aid}: missing {k}")
            qc = list(e.get("qc") or [])
            bad = set(qc) - set(EXTERNAL_QC)
            if bad:
                raise ConfigError(f"external assembly {sample}/{aid}: unknown qc {sorted(bad)}")
            rec = {k: _resolve(root, e[k]) for k in EXTERNAL_KEYS}
            rec.update(sample=sample, aid=aid, akey=f"imported-{aid}", qc=qc,
                       description=e.get("description", ""), register=e.get("register"),
                       read_types=list(e.get("read_types") or []))
            out.setdefault(sample, {})[aid] = rec
    return out


def _canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


def method_fingerprint(m):
    """Hash of the scientific content of a method (excludes description/file path)."""
    core = {k: v for k, v in m.items() if not k.startswith("_") and k != "description"}
    return hashlib.sha256(_canonical(core).encode()).hexdigest()[:12]


def file_hash(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()[:12]


def validation_identity(cfg, code_files):
    """Identity of everything that determines validation results and summaries.

    Changes to science.yaml content, any method definition, or validation code
    give a new identity; resources and tool paths do not."""
    parts = {
        "science": hashlib.sha256(_canonical(cfg["science"]).encode()).hexdigest()[:12],
        "methods": {mid: method_fingerprint(m) for mid, m in sorted(cfg["methods"].items())},
        "code": {os.path.relpath(p, cfg["pipeline_dir"]): file_hash(p) for p in sorted(code_files)},
    }
    vid = hashlib.sha256(_canonical(parts).encode()).hexdigest()[:10]
    return vid, parts


def validate(cfg):
    res = cfg["results_dir"]
    for d in cfg["legacy_dirs"] + [cfg["data_dir"]]:
        if _inside(res, d) or _inside(d, res):
            raise ConfigError(f"results_dir {res} must not overlap {d}")
    sci = cfg["science"]
    if not isinstance(sci.get("science_version"), int):
        raise ConfigError("science.yaml: science_version must be an integer")
    policy = sci.get("method_policy") or {}
    for key in ("ont_present", "ont_absent"):
        mid = policy.get(key)
        if mid not in cfg["methods"]:
            raise ConfigError(f"science.yaml: method_policy.{key} = {mid!r} is not a defined method")
    if "ont" not in cfg["methods"][policy["ont_present"]]["read_types"]:
        raise ConfigError("method_policy.ont_present must name a method that uses ONT reads")
    if "ont" in cfg["methods"][policy["ont_absent"]]["read_types"]:
        raise ConfigError("method_policy.ont_absent must name a HiFi-only method")
    gs = sci.get("genome_size_bp")
    if not isinstance(gs, int) or gs <= 0:
        raise ConfigError("science.yaml: genome_size_bp must be a positive integer")
    for sid, entry in cfg["samples"].items():
        if not re.match(SAMPLE_ID_RE, sid):
            raise ConfigError(f"samples.yaml: invalid sample ID {sid!r}")
        entry = entry or {}
        ov = entry.get("method_override")
        if ov:
            if ov not in cfg["methods"]:
                raise ConfigError(f"samples.yaml: {sid}: unknown method_override {ov}")
            if not entry.get("override_reason"):
                raise ConfigError(f"samples.yaml: {sid}: method_override requires override_reason")
        for name, why in (entry.get("exclude_inputs") or {}).items():
            if not why:
                raise ConfigError(f"samples.yaml: {sid}: exclude_inputs[{name}] needs a reason")
    for tool in ("hifiasm", "seqkit", "samtools"):
        if tool not in cfg["tools"]:
            raise ConfigError(f"tools config: missing required tool '{tool}'")


def rule_resources(cfg, rule):
    rules = cfg["resources"].get("rules") or {}
    out = dict(rules.get("default") or {"threads": 1, "mem_mb": 4000})
    out.update(rules.get(rule) or {})
    return out
