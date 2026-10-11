# Shared setup: configuration, input discovery, path helpers.
import json
import re
import shlex

from snakemake.exceptions import WorkflowError

from g2f import discovery
from g2f.config import ConfigError, load_config, rule_resources

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

wildcard_constraints:
    sample="|".join(map(re.escape, SAMPLES)) or "(?!)",
    method="|".join(map(re.escape, CFG["methods"])),
    iid=r"[A-Za-z0-9_.-]+",
    ctg="p_ctg|a_ctg",
    minlen=r"\d+",
    lineage=r"[A-Za-z0-9_.-]+",


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


shell.prefix(f"export TMPDIR={shlex.quote(CFG['tmpdir'])}; mkdir -p \"$TMPDIR\"; ")
