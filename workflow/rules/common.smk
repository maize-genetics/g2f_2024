# Shared setup: configuration, input discovery, path helpers.
import json
import re
import shlex

from snakemake.exceptions import WorkflowError

from g2f import discovery
from g2f.config import (ConfigError, load_config, method_fingerprint, rule_resources,
                        validation_identity)

try:
    CFG = load_config(config, PIPELINE_DIR)
except ConfigError as exc:
    raise WorkflowError(f"Configuration error: {exc}")

SCI = CFG["science"]
RES = CFG["results_dir"]
GENOME_SIZE = SCI["genome_size_bp"]
K = SCI["kmer"]["k"]
PY = sys.executable
TOOL = f"{shlex.quote(PY)} {shlex.quote(os.path.join(PIPELINE_DIR, 'scripts', 'g2f_tool.py'))}"
SAMPLES_YAML = os.path.join(PIPELINE_DIR, config.get("samples", "config/samples.yaml"))
if not os.path.exists(SAMPLES_YAML):     # the samples file is optional: read as empty
    SAMPLES_YAML = os.devnull

# ---------------------------------------------------------------- discovery
INPUTS, ISSUES = discovery.discover(CFG)
for _i in ISSUES:
    if _i.severity == "warning":
        sys.stderr.write(_i.render() + "\n")
_errors = [i for i in ISSUES if i.severity == "error"]
if _errors:
    raise WorkflowError(
        "Input discovery found problems that need a human decision. Nothing was run.\n\n"
        + "\n\n".join(i.render() for i in _errors)
        + "\n\nReview inputs with: python3 scripts/g2f_tool.py discover")

SAMPLES = sorted(INPUTS)
INV = {s: discovery.inventory_records(s, INPUTS[s], CFG) for s in SAMPLES}
METHOD = {s: INV[s]["selected_method"] for s in SAMPLES}


def inventory_path(s):
    return f"{RES}/{s}/inputs/inventory.json"


# Inputs must not change silently underneath existing results.
for _s in SAMPLES:
    _p = inventory_path(_s)
    if os.path.exists(_p):
        with open(_p) as _fh:
            _old = json.load(_fh)
        if discovery.stable_view(_old) != discovery.stable_view(INV[_s]):
            raise WorkflowError(
                f"Inputs for sample {_s} differ from the inventory recorded at\n  {_p}\n"
                "which existing results were built from (files added, removed, changed, or the "
                "selected method changed). Review with `python3 scripts/g2f_tool.py discover "
                f"--only {_s}`. If the change is intended, move results/{_s}/ aside so the sample "
                "is processed from scratch; the workflow will not mix input sets.")

# ---------------------------------------------------------------- identity
# Validation identity: science.yaml content + all method definitions + the
# code that validates and summarizes. A change gives a new identity, which
# (a) re-runs the cheap validation/summary jobs (they depend on the stamp
# file below) and (b) writes validated markers and summary tables under new
# names, so an outdated validation can never be mistaken for a current one.
# Heavy jobs (conversion, k-mer counting, hifiasm, BUSCO, Merqury) do not
# depend on it, and resource/tool-path settings are not part of it.
_SCRIPTS = os.path.join(PIPELINE_DIR, "scripts")
VALIDATION_CODE = [os.path.join(_SCRIPTS, "g2f", f) for f in
                   ("validate.py", "gfa.py", "kmer.py", "summary.py")] + \
                  [os.path.join(_SCRIPTS, "g2f_tool.py")]
VID, VID_PARTS = validation_identity(CFG, VALIDATION_CODE)
STAMP = f"{RES}/.identity/validation-{VID}.json"
METHOD_FP = {m: method_fingerprint(d) for m, d in CFG["methods"].items()}


def method_record_path(s, m):
    return f"{RES}/{s}/assembly/{m}/{s}.{m}.method.json"


# Existing assemblies must have been produced by the current method definition.
for _s in SAMPLES:
    for _m in CFG["methods"]:
        _p = method_record_path(_s, _m)
        if os.path.exists(_p):
            with open(_p) as _fh:
                _fp = json.load(_fh)["fingerprint"]
            if _fp != METHOD_FP[_m]:
                raise WorkflowError(
                    f"The method definition {_m} was changed after it was used for sample {_s} "
                    f"(recorded fingerprint {_fp}, current {METHOD_FP[_m]}; see {_p}). Method "
                    "files are immutable: restore the original file, or create a new method_id "
                    "for the changed parameters.")

# ------------------------------------------------- imported (legacy) assemblies
_only = set(CFG["only_samples"])
IMPORTED = {}
for _s, _entries in CFG["external_assemblies"].items():
    if _only and _s not in _only:
        continue
    for _aid, _rec in _entries.items():
        IMPORTED[(_s, _rec["akey"])] = _rec


def import_record_path(s, aid):
    return f"{RES}/{s}/imported/{aid}/import.json"


for (_s, _akey), _rec in IMPORTED.items():
    _p = import_record_path(_s, _rec["aid"])
    if os.path.exists(_p):
        with open(_p) as _fh:
            _old = json.load(_fh)
        for _k, _src in (("gfa", _rec["p_ctg_gfa"]), ("fa", _rec["p_ctg_fa"]),
                         ("log", _rec["hifiasm_log"])):
            _f = _old["files"][_k]
            _ok = (os.path.isfile(_src) and _f["path"] == _src
                   and os.path.getsize(_src) == _f["size"]
                   and discovery.fingerprint(_src, _f["size"]) == _f["fingerprint"])
            if not _ok:
                raise WorkflowError(
                    f"Registered historical file for {_s}/{_rec['aid']} changed or moved since "
                    f"import: {_src}. Historical assemblies must stay unchanged; investigate "
                    f"before continuing (import record: {_p}).")

AKEYS = list(CFG["methods"]) + sorted({k for _, k in IMPORTED})

wildcard_constraints:
    sample="|".join(map(re.escape, sorted(set(SAMPLES) | {s for s, _ in IMPORTED}))) or "(?!)",
    method="|".join(map(re.escape, CFG["methods"])),
    akey="|".join(map(re.escape, AKEYS)),
    aid="|".join(map(re.escape, sorted({r["aid"] for r in IMPORTED.values()}))) or "(?!)",
    iid=r"[A-Za-z0-9_.-]+",
    ctg="p_ctg|a_ctg",
    minlen=r"\d+",
    lineage=r"[A-Za-z0-9_.-]+",
    reviewid=r"[0-9a-f]+|none",


# ---------------------------------------------------------------- helpers
def R(rule):
    return rule_resources(CFG, rule)


def used(s, tech):
    return [i for i in INPUTS[s] if i.tech == tech and i.status == "use"]


def get_input(s, tech, iid):
    for i in used(s, tech):
        if i.input_id == iid:
            return i
    raise WorkflowError(f"unknown {tech} input {iid} for sample {s}")


def raw_input(s, tech, iid):
    """Path of an original input file for rule inputs, or nothing once its local
    raw copy was released (deleted after verified preprocessing; see
    g2f.discovery.released_inputs). Jobs whose outputs exist are then simply
    not rerun; a job that would need the deleted file fails via released_guard."""
    i = get_input(s, tech, iid)
    return [] if discovery.is_released(i) else i.path


def released_guard(s, tech, iid, cmd):
    i = get_input(s, tech, iid)
    if not discovery.is_released(i):
        return cmd
    return (f"echo 'Cannot run: the local raw copy {i.path} was released (deleted after verified "
            "preprocessing) and the prepared file derived from it is missing. Restore the raw "
            "file from the source storage node, then rerun.' >&2; exit 1")


def method(m):
    return CFG["methods"][m]


def uses_ont(m):
    return "ont" in method(m)["read_types"]


def tool(name):
    try:
        return CFG["tools"][name]["cmd"]
    except KeyError:
        raise WorkflowError(f"tool '{name}' is not defined in the tools configuration")


def prov_path(s, step, name):
    return f"{RES}/{s}/provenance/{step}/{name}.json"


def wrap(step, cmd, sample, prov, tools=(), method_id=None, expect=None, inputs=(), outputs=()):
    """Shell command that runs `cmd` through the provenance recorder."""
    parts = [TOOL, "run", "--prov", shlex.quote(prov), "--step", step, "--sample", sample]
    if method_id:
        parts += ["--method", method_id]
    for t in tools:
        parts += ["--version", shlex.quote(f"{t}={CFG['tools'][t]['version_cmd']}")]
    for t, v in (expect or {}).items():
        parts += ["--expect", shlex.quote(f"{t}={v}")]
    if inputs:
        parts += ["--inputs"] + [shlex.quote(str(p)) for p in inputs]
    if outputs:
        parts += ["--outputs"] + [shlex.quote(str(p)) for p in outputs]
    return " ".join(parts + ["--", shlex.quote(cmd)])


def q(path):
    return shlex.quote(str(path))


# Optional cleanup of large, regenerable intermediates (computational setting,
# config/resources/*.yaml `cleanup:`; never changes results). Off by default.
CLEANUP = CFG["resources"].get("cleanup") or {}


def _flag(v):
    return str(v).strip().lower() in ("1", "true", "yes", "keep")


# Retention of the converted HiFi FASTQ (a prepared hifiasm input): kept unless
# set to false, in which case Snakemake deletes it once read QC, the prepared-
# input check, meryl and hifiasm have used it. --config retain_prepared_hifi=
# overrides the resources file (the launcher sets it per genotype from the
# disk projection).
RETAIN_HIFI = _flag(config.get("retain_prepared_hifi", CLEANUP.get("retain_prepared_hifi", True)))


shell.prefix(f"export TMPDIR={shlex.quote(CFG['tmpdir'])}; mkdir -p \"$TMPDIR\"; ")
