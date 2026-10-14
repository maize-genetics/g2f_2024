# Genotype fingerprinting (optional k-mer sample-identity QC)

**Question.** Is each ONT read set genetically consistent with the HiFi data
of the genotype it is labelled as, rather than another genotype in the
collection? The HiFi labels are treated as the reference; the ONT labels are
the hypotheses being tested. Every ONT read set is compared with **every**
HiFi genotype, not only its presumed partner.

**What it does not do.** It never moves, renames or relabels files and
does not feed into assembly. Its output is evidence for a human decision
(`samples.yaml`). Any future automatic reorganization would be a separate,
explicitly approved step.

Code: `scripts/g2f/fpsketch.py` (sketching), `scripts/g2f/fingerprint.py`
(statistics), `scripts/g2f/fpanalysis.py` (driver/report),
`scripts/g2f_fingerprint.py` (CLI), `workflow/rules/fingerprint.smk`
(Snakemake target `fingerprint`). Settings: `config/fingerprint.yaml`.
Tests: `tests/test_fingerprint.py`.

## Running it

Inside the workflow (uses the discovered inputs; independent of assembly):

```bash
snakemake fingerprint --profile profiles/local --cores 32 --resources mem_mb=200000
# a subset, or other settings (target first, then --config):
snakemake fingerprint --profile profiles/local --cores 32 \
    --config only_samples=SAMPLE_A,SAMPLE_B fingerprint=config/fingerprint_pilot.yaml
```

Outputs (never overwritten; the directory names are hashes of the settings
and code that produced them):

```
results/fingerprint/sketches-<SID>/<LINE>/<hifi|ont>/<input_id>.{sig.gz,json,prov.json}
results/fingerprint/report-<RID>/
    assignments.tsv                 per ONT run / replicate half / pooled genotype:
                                    best match, pairwise support, classification
    similarity_matrix.tsv           ONT x HiFi private-k-mer score (primary k)
    qc_summary.md                   human-readable report
    parameters_and_provenance.json  settings, tool versions, inputs, sampling, resources
    analysis.log                    progress + /usr/bin/time resource record
    k21/ k31/                       full results per k: scores_long.tsv,
                                    pairwise_long.tsv, containment_matrix.tsv,
                                    controls.tsv, panel.tsv, hifi_relatedness.tsv,
                                    calibration.json
```

Without the workflow (any lab, any layout; ONT label may be empty = unknown):

```bash
FPY=python3      # an interpreter with sourmash installed (the fingerprint_python tool)
$FPY scripts/g2f_fingerprint.py sketch --input reads.bam --tech hifi --label LINE_A \
    --run-id run1 --ksizes 21,31 --scaled 2000 --max-bases 15000000000 \
    --out-sig sk/A.run1.sig.gz --out-json sk/A.run1.json
$FPY scripts/g2f_fingerprint.py sketch --input ont_run.tar --tech ont --label "" \
    --run-id ontX --ksizes 21,31 --scaled 2000 --max-bases 10000000000 --min-read-len 1000 \
    --out-sig sk/ontX.sig.gz --out-json sk/ontX.json
python3 scripts/g2f_fingerprint.py analyze --records sk/*.json --outdir report --ksizes 21,31
```

One sketch job per file is independent, so many files can be sketched in
parallel (e.g. GNU parallel). At least two HiFi genotypes are needed;
leave-one-run-out controls need a genotype with two or more HiFi runs.

## Method

### Sketching (streaming, bounded, reproducible)

* Each **original** read file is streamed once: PacBio BAM via
  `samtools fasta -F 0x900`, FASTQ/FASTA (gz or not), and ONT tar archives
  (pass members streamed in place; fail members excluded). No reads are
  converted or copied to disk.
* **Bounded sampling:** the first `*_max_bases_per_input` bases. PacBio BAMs
  are in ZMW order and ONT reads in acquisition order, neither related to
  genome position, so a prefix is an unbiased genomic sample. A seeded
  read-level subsample (`read_fraction`) is available for inputs whose order
  might be genomic. Reads shorter than `min_read_len` (1 kb) are skipped.
* **FracMinHash** (sourmash 4.9.3; MurmurHash3 of canonical k-mers, seed 42):
  keeps every k-mer whose hash is < 2^64/`scaled`. With `scaled` = 2000 this
  is a deterministic 1/2000 sample of k-mer space, the same for every input,
  so sketches from different files are directly comparable and reusable
  (~1.1 M sampled 21-mers per maize genome). Abundances are tracked.
* **Replicate halves:** each read goes to half A or B by a seeded hash of its
  name; halves are stored separately (the full sketch is their sum). This
  gives within-run technical replicates without extra cost.
* Each sketch has a JSON record with the input's size and fingerprint, the
  exact sampling (reads/bases used, whether capped, a SHA-256 of the used read
  names), parameters, sourmash version and run time.

### k-mer size

21 is primary, with 31 computed in the same pass as a robustness check.
Cross-platform matching needs k-mers that are error-free in the noisier
reads. With ONT R10.4.1 HAC v4.2.0 reads (mean quality ~Q14 measured on this
project's data), about 0.96^21 = 42 % of true 21-mers are error-free, against
28 % for 31-mers. A random 21-mer is still essentially unique in a 2.3 Gb
genome (4^21 = 4.4 x 10^12). Errors only remove k-mers or create novel ones;
they cannot make one genotype's k-mers look like another's except by
recreating the other allele exactly, and that rate is measured empirically
(below).

### Reference panel (HiFi)

For each HiFi genotype all its runs are summed. Then:

* **present(G):** every sampled hash seen at least once (including errors).
* **solid(G):** hashes with abundance from the end of the error peak (first
  minimum of the abundance spectrum, at least 2) to 3 x the main single-copy
  peak. This removes sequencing errors and also high-copy repeats, organelle
  and high-coverage contaminant k-mers, which are not informative about
  nuclear genotype. The maize genome is mostly repeats, so this matters.
* **distinct(G vs H):** solid(G) hashes never seen in H's reads, i.e. G's
  alleles and G-specific sequence where the two genotypes differ.
* **private(G):** solid(G) hashes seen in no other genotype.
* **core:** hashes solid in every genotype.

The strict absence rule (not seen even once) requires enough HiFi depth:
at 6.5x, a k-mer is missed with probability ~e^-6.5 = 0.15 %. Between maize
inbreds, distinct sets are a large fraction of solid k-mers, so this
contaminates them by well under 1 %.

### Scores for a query read set Q

* **sensitivity** s(Q) = |core ∩ Q| / |core|: the fraction of k-mers shared
  by all genotypes that Q recovers. This captures coverage and per-k-mer
  accuracy together, so scores below are comparable between deep HiFi and
  shallow, noisy ONT, and between runs of very different size.
* **pairwise** has(Q; G vs H) = |distinct(G vs H) ∩ Q| / |distinct(G vs H)| / s(Q):
  about 1 if Q is G, about 0 (the error null) if Q is H.
* **private score** |private(G) ∩ Q| / |private(G)| / s(Q) (the similarity
  matrix: one interpretable number per ONT x HiFi cell) and **containment**
  |solid(G) ∩ Q| / |solid(G)| / s(Q) (genome-wide shared k-mer fraction;
  descriptive).

Q is **consistent with G** when, against every other genotype H, Q carries
G's distinct alleles (has(G vs H) above the threshold) and lacks H's
(has(H vs G) below it), with the whole 99.9 % binomial (Wilson) interval of
each score on the correct side. Pairwise tests resolve close relatives that
share most of their k-mers (few private k-mers), and they identify mixtures:
a mixture carries the distinct alleles of two genotypes.

### Calibration from the data (no fixed similarity cut-off)

* **Positive controls (truth known):** each HiFi run of a genotype with two
  or more runs, scored against a panel in which that genotype is represented
  only by its *other* runs (leave-one-run-out: independent SMRT cells).
  Plus, for every genotype, replicate half A scored against a reference built
  from half B (same libraries: **technical reproducibility only**, not
  independent genotype validation; reported separately as `split_half`).
  These give the "has" scores of true matches.
* **Negative controls:** the same control reads scored for the alleles of
  every *other* genotype.
* **Threshold** = midpoint between the largest negative-control score and
  the smallest positive-control score. If the two distributions overlap,
  the data cannot discriminate and **every** query is reported as ambiguous.
* **ONT null (reported, not used for the threshold):** each ONT run's scores
  for the genotypes other than its best match. This is the platform-matched
  error null. It is excluded from the threshold because a genuinely mixed ONT
  sample would inflate it.

Control references have less depth than the full panel (fewer runs, or half
the reads). That inflates their negative scores slightly (missed k-mers look
"distinct"), so the threshold is conservative.

### Classification

| class | meaning |
|---|---|
| `supported` | consistent with the expected genotype and no other |
| `potentially_mismatched` | consistent only with a different genotype, or the reads contradict the expected genotype (they lack its alleles and carry another's) |
| `ambiguous` | anything else, with the reason: insufficient data (an interval spans the threshold), several genotypes' alleles present (mixture/contamination, `possible_admixture`), no panel genotype consistent (genotype absent from panel?), genotypes unresolvable (too few informative k-mers), or the data do not discriminate at all |
| `assigned` | (no expected label) consistent with exactly one genotype |

Evidence is reported per ONT **run**, per replicate **half** of each run,
and **pooled** over all runs with the same expected label. `consistent_best`
says whether all runs and halves of a label agree. A pooled set mixing a
correct and a mislabelled run is (correctly) not `supported`. k = 31 results
are given alongside (`class_k31`); disagreement between k values is listed
in `qc_summary.md`.

### Factors accounted for

| factor | handling |
|---|---|
| ONT errors / read length | k = 21; sensitivity normalization; error null measured, not assumed |
| unequal coverage between runs | normalization by core recovery; intervals widen at low coverage |
| repeats, organelles | solid band excludes k-mers above 3x the single-copy peak |
| HiFi errors | solid band starts after the error peak |
| closely related genotypes | pairwise distinct sets; unresolvable pairs reported, not forced |
| contamination / mixtures | a second genotype's alleles above the control null → ambiguous with `possible_admixture`. Contamination by non-maize DNA only dilutes coverage; it is not detected as such |
| genotype not in panel | no consistent genotype → ambiguous (or `potentially_mismatched` if the expected genotype is clearly contradicted) |

### Limitations

* Only genotypes with HiFi data can be identified. An ONT set from a line
  outside the panel can only be reported as "not consistent with any".
* Two genotypes that are genetically identical (e.g. the same seed lot under
  two names) cannot be told apart. Their pairs are reported as unresolved.
* Prefix sampling assumes read order unrelated to genome position (true for
  PacBio BAM and MinKNOW output, not for coordinate-sorted BAM; use
  `read_fraction` there).
* The control threshold depends on having controls. With a single HiFi run
  per genotype, only split-half (technical) controls exist.

## Resources

The analysis is standard-library Python. Its memory is dominated by the
panel's hash sets: a few GB for about ten maize genotypes at scaled = 2000. It is
recorded in `parameters_and_provenance.json` and `analysis.log`. Sketching
takes one hashing thread plus decompression at ~12-20 Mbp/s per k-size pair,
so a 15 Gb HiFi prefix takes ~15-20 min. Sketches are 10-60 MB each.
