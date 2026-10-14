"""CPU blocks for the full-scale launcher: a hard thread cap per genotype.

Each concurrently running genotype is pinned (taskset) to its own block of
`threads` logical CPUs taken from the CPUs this process may use, one per
physical core before any hyperthread sibling. Every tool and
helper inherits the pin, so the run never uses more than jobs x threads CPUs,
whatever extra threads a tool starts.
"""
import os


def compress(cpus):
    """[0, 1, 2, 5] -> '0-2,5' (taskset -c list format)."""
    out, cpus = [], sorted(cpus)
    i = 0
    while i < len(cpus):
        j = i
        while j + 1 < len(cpus) and cpus[j + 1] == cpus[j] + 1:
            j += 1
        out.append(str(cpus[i]) if i == j else f"{cpus[i]}-{cpus[j]}")
        i = j + 1
    return ",".join(out)


def core_first(cpus, sysfs="/sys/devices/system/cpu"):
    """Order CPUs so one logical CPU of every physical core comes before any
    second hyperthread of a core (whatever the host's CPU numbering), then by
    number. Without topology information: numeric order."""
    rank = {}
    for c in cpus:
        try:
            with open(f"{sysfs}/cpu{c}/topology/thread_siblings_list") as fh:
                sib = sorted(int(x) for part in fh.read().strip().split(",")
                             for x in _expand(part))
            rank[c] = sib.index(c)
        except (OSError, ValueError):
            rank[c] = 0
    return sorted(cpus, key=lambda c: (rank[c], c))


def _expand(part):
    a, _, b = part.partition("-")
    return range(int(a), int(b or a) + 1)


def slots(jobs, threads, allowed=None, sysfs="/sys/devices/system/cpu"):
    """One disjoint CPU list per job slot (slot n = list n-1), physical cores first."""
    if jobs < 1 or threads < 1:
        raise ValueError("jobs and threads must be >= 1")
    allowed = core_first(sorted(os.sched_getaffinity(0) if allowed is None else allowed), sysfs)
    need = jobs * threads
    if need > len(allowed):
        raise ValueError(f"{jobs} genotypes x {threads} threads = {need} CPUs requested, but only "
                         f"{len(allowed)} are available to this process; lower --jobs or --threads")
    return [compress(allowed[k * threads:(k + 1) * threads]) for k in range(jobs)]
