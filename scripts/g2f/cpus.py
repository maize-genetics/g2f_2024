"""CPU blocks for the full-scale launcher: a hard thread cap per genotype.

Each concurrently running genotype is pinned (taskset) to its own block of
`threads` logical CPUs taken from the CPUs this process may use. Every tool and
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


def slots(jobs, threads, allowed=None):
    """One disjoint CPU list per job slot (slot n = list n-1)."""
    if jobs < 1 or threads < 1:
        raise ValueError("jobs and threads must be >= 1")
    allowed = sorted(os.sched_getaffinity(0) if allowed is None else allowed)
    need = jobs * threads
    if need > len(allowed):
        raise ValueError(f"{jobs} genotypes x {threads} threads = {need} CPUs requested, but only "
                         f"{len(allowed)} are available to this process; lower --jobs or --threads")
    return [compress(allowed[k * threads:(k + 1) * threads]) for k in range(jobs)]
