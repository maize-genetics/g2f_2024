# Genotype fingerprinting: supplementary diagnostics (2026-10-12)

This report re-examines the fingerprinting results in
`docs/fingerprinting_results_2026-10-12.md`, which is preserved unchanged. That
report covered `results/fingerprint/report-742a05a563/` (sketches
`sketches-90f6312d72`, code ff68d36). The question here is what the evidence
supports for GEMS-0061 and GEMS-0113 in particular. The production algorithm,
its thresholds and all classifications are unchanged. No input, assembly or
pipeline output was modified.

* Diagnostics: `results/fingerprint/diagnostics-2026-10-12/` (3.6 MB). The
  summary and per-k tables (everything except the per-read tables) are copied
  into the repository under `docs/fingerprinting_diagnostics_2026-10-12/`,
  with SHA-256 checksums of every output file in `MANIFEST.sha256`.
* Code (read-only diagnostics): `scripts/g2f_fp_diagnostics.py`,
  `scripts/g2f_fp_diagnostics_summary.py`. Tests: `tests/test_fp_diagnostics.py`.
* The panel is rebuilt with the unchanged `g2f.fingerprint` code and the
  report's own parameters.

**Decision (project lead, 2026-10-12).** The existing assemblies proceed
unchanged for all nine genotypes, GEMS-0061 and GEMS-0113 included. The
possible HiFi–ONT source difference of those two lines is documented here
and checked after assembly (`docs/post_assembly_qc_plan.md`). HiFi-only
assemblies and further mixing experiments are not run unless that QC
justifies them.

## 1. Conclusions

1. **The genotype assignments are supported for all nine genotypes.** Every usable ONT
   run, and both replicate halves of every run, matches its labelled HiFi
   genotype. The margin is large: the best alternative genotype's net
   support is at most −0.54. The identity threshold sits far outside the ONT
   noise: ONT "has" scores are 0.66–0.99, while the ONT-specific null is at
   most 0.056 per run (0.073 pooled). GEMS-0219 run 9 is insufficient data,
   and it is excluded from assembly (verified, §5).
2. **No evidence was found that DNA from another panel genotype is mixed
   into any ONT library.** For any panel line, about 1–3 % whole-genome
   admixture would have been detectable (2–5 % for the shallow PHBW8 set)
   (§4.4; indicative limits). This does not exclude contamination below
   those limits, by lines outside the panel, or by non-maize DNA.
3. **GEMS-0061 and GEMS-0113: the ONT and HiFi data show genetic
   discordance between sequencing sources.** The foreign alleles
   are real sequence. They replicate across flow cells, co-occur on single
   reads, and are unchanged between k = 21 and k = 31. They are present at
   near-genomic dosage (about 0.7–1.0 of core depth), yet only at a minority
   of loci. In the other direction, the ONT lacks 11–16 % of the HiFi
   sample's informative alleles (27 % of its half-dose alleles). Admixture
   alone adds alleles but does not remove them, and low-level contamination
   would give dosage close to the admixture fraction (spike-in calibrated,
   §4.3). Low-level admixture is therefore not a sufficient explanation.
   The data are most consistent with genetically different individuals, or
   pools of individuals, of a line that is not fixed. The HiFi samples of
   both lines are themselves measurably heterozygous or heterogeneous. The
   biological provenance has not been established.
4. **GEMS-0113's FBLL signal has two parts.** About half of FBLL's elevated
   score is a panel-composition effect, and the other half is genuinely
   FBLL-specific. FBLL's set of alleles distinct from GEMS-0113 is the
   smallest in the panel and the most enriched for widely shared alleles.
   The FBLL-specific excess is present in both flow cells and at both k. It
   occurs at 0.45–0.8 of genomic dosage over a limited part of the genome.
   This is more consistent with FBLL-like haplotype segments in the ONT
   source plant(s) than with low-level admixture of FBLL DNA: a 2 % FBLL
   spike-in gives similar scores but a dosage of about 0.05.
5. **Status of the contamination question.** Contamination of ONT libraries
   by another panel genotype is not supported by these data, but it is not
   categorically excluded. Still unresolved: (a) which DNA
   source represents the registered GEMS-0061 and GEMS-0113 lines; (b) the
   biological cause (residual heterozygosity, a heterogeneous seed lot, or
   historical outcrossing); (c) contamination below about 1–3 %, a
   minor contaminant hidden under the larger genetic signal (§4.4), and
   contamination by non-maize DNA, which this method does not test. Answering
   these needs seed-lot or plant provenance, or external genotype data, not
   more k-mer analysis.

Normalized scores are not contamination percentages. As §4.3 shows, the same
"lacks" value of about 0.17 arises from 2 % admixture or from full-dosage
genetic differences.

## 2. GEMS-0061

**Which genotypes contribute.** Foreign-allele scores against every other line
(pooled ONT, k = 21; hits / distinct set size):

| line | 3AZA1 | FBLL | GEMS-0113 | GEMS-0219 | GEMS-0227 | PHB47 | PHBW8 | PHP02 |
|---|---|---|---|---|---|---|---|---|
| score | 0.167 | 0.195 | 0.196 | 0.191 | 0.172 | 0.171 | 0.166 | 0.150 |
| hits | 10,110 | 7,495 | 10,825 | 10,448 | 9,623 | 10,551 | 9,283 | 9,845 |
| distinct set | 64,555 | 40,984 | 58,840 | 58,318 | 59,598 | 65,658 | 59,502 | 69,874 |

These eight scores are not eight independent signals. The pairwise distinct
sets overlap heavily, because a common allele is distinct for many lines at
once. Instead, recovery of foreign alleles rises steeply with the number of
panel lines that carry them (`k21/class_recovery.tsv`): 0.092 (one carrier)
→ 0.117 → 0.136 → 0.158 → 0.203 → 0.257 → 0.302 → 0.336 (all eight
carriers). In pure lines this curve is flat at the error null (0.04–0.06).
That gradient alone predicts every line's score to within ±0.02, FBLL
included (0.195 observed vs 0.199 predicted). No line-specific signal of
the kind that contamination by a panel line would produce was detected.

Testing within rarely shared alleles (`line_specific_excess.tsv`) gives a
modest excess for GEMS-0113 (+45 %) and GEMS-0219 (+19 %) haplotypes
(k = 31: +53 %, +26 %). It is at genomic dosage (0.64–0.76, pooled) and
present in all three runs. FBLL shows none. This is consistent with shared
ancestry among GEMS lines in some regions; it lacks the low dosage expected
from DNA mixing.

**Real sequence, not errors** (`replication.tsv`, `read_cooccurrence.tsv`):

* Foreign k-mers seen in one run reappear in another run at 0.68–0.79 of the
  core-k-mer replication rate. This includes run pairs across flow cells
  PAQ87161 and PAU07791 (0.70–0.77). The corresponding rate between runs is
  0.32–0.57 for pure-line null hits and 0.05–0.14 for novel (error)
  k-mers.
* Per read (run 8, all 1.37 Gb, 233,662 reads): reads carrying two or more
  foreign markers occur 2.36× more often than expected under independence
  (stratified by read quality). Reads whose informative markers are all
  foreign occur 6.1× more often than expected (758 reads). In pure-line
  controls the ratios are 0.65× and 0.90×.
* The signal is spread across 6,409 reads. The top 1 % of those reads hold
  only 4.1 % of the foreign markers, and the maximum is 12 markers on one
  read. Foreign-carrying reads have almost normal quality (median Q15.4 vs
  15.9). In pure lines the null hits sit on poorer reads (Q13.8 vs 15.2;
  Q16.4 vs 19.6).
* The excess over the ONT null is the same at k = 21 and k = 31 (run level
  0.074–0.124 vs 0.068–0.122). An error-driven signal would shrink at k = 31
  (expected relative factor about 0.96^10 = 0.66).

**Dosage: consistent with genetic difference rather than admixture** (`generic_foreign_excess.tsv`).
After subtracting the ONT null, the excess foreign k-mers (15,343 pooled,
k = 21) are at 0.76–0.89 of core depth (k = 31: 0.73–0.82). This is close
to the 0.81 dosage of GEMS-0061's own single-copy alleles. Whole-genome
admixture at 10 % would give about 0.1 (spike-in, §4.3).

**Missing alleles.** The ONT recovers only 84 % of the HiFi sample's
informative alleles. For the half-dosage ("heterozygous-like") alleles in
the HiFi spectrum, recovery is 73 %. Pure lines recover 96–100 %. Admixture
never removes alleles: in 50 % spike-ins this "has" score stays at ≥ 0.91.

**The HiFi data show the same phenomenon.** GEMS-0061's HiFi spectrum has a
clear half-dose shoulder: mass at abundance 4–6 relative to 8–13 is 0.444,
vs 0.20 in PHBW8 and PHP02 at the same 25 Gb. It carries 41,587 low-dose
informative k-mers, against 15–17k in clean inbreds. These low-dose alleles
are rarer in the panel than the line's single-copy alleles (mean carriers
3.5 vs 4.1). The split-half control is the highest of all split-halves
(lacks 0.0995). That is expected when heterozygous k-mers at about 2× depth
per half are missed. No independent second HiFi run exists, so HiFi
within-sample heterogeneity cannot be separated from heterozygosity.

**Run consistency.** All three runs (two flow cells) and all six halves give
the same picture: has 0.66–0.80, lacks 0.14–0.20, the same carrier gradient,
and the same GEMS-0113/0219 excess. Scores rise only slightly with
sensitivity (0.19 → 0.94). This points to one biological DNA source shared
by all libraries, not a run-specific or computational artifact.

**Assessment.** No option is a clean fit, so each is stated separately:

* **A (cross-genotype contamination, mixed-genotype libraries): not
  supported.** The dosage is genomic, alleles are missing, and no
  line-specific signature was found. A minor contaminant below the detection
  limit could still be present underneath the larger signal.
* **A′ (genuine reads from genetically different material of the same
  nominal line): best supported.** This could be a different plant or a pool
  of non-identical plants.
* **B (similarity / cross-reactivity): explains the pattern across lines,
  not the signal itself.** It accounts for the uniform spread across lines
  and the FBLL/GEMS-0113/GEMS-0219 ranking.
* **C (ONT errors or coverage): not supported** as the main cause
  (replication, read-level co-occurrence, k-independence).
* **D (insufficient evidence): applies only to the finer question.** The data
  cannot tell one different plant from a pool, or residual heterozygosity
  from seed-lot heterogeneity or outcrossing.

## 3. GEMS-0113 and FBLL

| | FBLL | other lines |
|---|---|---|
| pooled lacks, k = 21 / k = 31 | 0.178 / 0.166 | 0.078–0.108 / 0.062–0.094 |
| run 7 / run 7b (k = 21) | 0.163 / 0.148 | 0.071–0.100 / 0.061–0.089 |
| absolute hits (pooled) | 5,675 | 3,798–5,893 |
| distinct-from-GEMS-0113 set | **39,805** | 57,134–74,973 |
| mean panel carriers of that set | **4.37** | 3.11–3.97 |
| singletons (only that line carries) | **3.3 %** | 5–31 % |

**Why FBLL recurs.** FBLL is the panel's most "central" genotype. It has only
1,010 private k-mers, and 85 % of its solid k-mers occur in GEMS-0113. Any
non-HiFi signal made of common alleles is therefore concentrated in FBLL's
small distinct set. In absolute numbers, GEMS-0113 ONT carries about as many
GEMS-0061-distinct (5,893) or PHB47-distinct (5,036) k-mers as FBLL-distinct
ones (5,675). The same FBLL inflation appears in HiFi-only controls where
ONT plays no part:

| HiFi control | FBLL score | other lines |
|---|---|---|
| GEMS-0219 leave-one-out run bc2041 | 0.114 | 0.038–0.082 |
| GEMS-0227 leave-one-out run bc2004 | 0.096 | 0.040–0.066 |
| GEMS-0061 split-half | 0.100 | 0.054–0.084 |
| GEMS-0113 split-half | 0.084 | 0.044–0.058 |

The allele-sharing gradient predicts FBLL at 0.115 (observed 0.178) and the
other lines within ±0.03.

**The FBLL-specific part is real.** Among rarely shared foreign alleles (≤ 3
carriers), FBLL-carried ones are recovered at 0.171 vs 0.046 for the rest.
That is an excess of 1,347 k-mers (+243 %; k = 31: +360 %, 2,295). It is
present in both independent flow cells (run 7: 1,033; run 7b: 589). No
other line has a pooled excess above +18 % (at most +31 % in a single run at
k = 31). FBLL-private k-mers are recovered at 0.094,
vs a pooled-ONT private null of median 0.030 and maximum 0.050 (76 of 1,010).
Reads with two or more foreign markers that are all FBLL-carried: 436
observed vs 119 expected (run 7b).

**Dosage argues against FBLL DNA admixture.** The FBLL-specific excess is at 0.45–0.50 of
genomic dosage (k = 21) and 0.65–0.79 (k = 31), consistent in both runs.
In-silico FBLL spike-ins into clean ONT sets reproduce GEMS-0113's scores
at about 2 % FBLL (GEMS-0227 host: lacks 0.166, class "ambiguous: FBLL").
At that level the dosage is about 0.05, and a dosage of 0.5 needs about 50 %
FBLL, which pushes the FBLL score above 0.9. No uniform FBLL admixture
yields both the observed excess and the observed dosage. A limited set of
genomic segments that are FBLL-like at about heterozygous-to-homozygous
dosage does yield both.

GEMS-0113's HiFi low-dose k-mers are *not* FBLL-enriched (FBLL carries 58 %
of them vs 68 % of single-copy ones). The FBLL-like segments therefore
appear specific to the ONT source.

GEMS-0113 also shows a smaller generic foreign signal like GEMS-0061's
(4,234 excess k-mers; dosage 0.99–1.44 at k = 21 and 0.79–1.04 at k = 31,
imprecise at this depth). The ONT lacks 11 % of the HiFi
sample's informative alleles and 27 % of its half-dose alleles. The HiFi
half/full ratio is 0.341.

**Assessment.**

* **A (FBLL DNA admixed in the reads): not supported** on dosage; FBLL
  admixture below about 1–2 % cannot be excluded.
* **B: explains about half of FBLL's raw score** (composition). It does not
  explain the FBLL-specific excess.
* **C: not supported** as the main cause.
* **A′ (best supported):** the ONT source material carries FBLL-like
  haplotypes in some regions where the HiFi source does not. That fits
  FBLL-related germplasm segregating in GEMS-0113.
* **D: applies to the cause.** The data cannot tell pedigree from outcrossing
  in a seed increase. The pedigree was not checked and should not be
  assumed.

## 4. Methodology and thresholds

### 4.1 HiFi thresholds applied to ONT

The identity threshold (0.518) is appropriate for identity. ONT positives
are ≥ 0.66, and the ONT null is ≤ 0.073.

The HiFi negative maximum (0.114) is **not** a valid contamination criterion
for ONT, for two reasons:

* It is not an error null. Its maximum comes from GEMS-0219 run bc2041,
  which genuinely differs from the other GEMS-0219 runs (LOO has 0.923,
  FBLL 0.114). Its median (0.019) is below the ONT null median (0.034,
  k = 21).
* More fundamentally, a lacks score above any null shows that alleles are
  present. It does not say how they got there: about 2 % admixture and
  full-dosage genetic differences give the same values (§4.3).

The current "possible_admixture" wording ("mixture or contamination?")
therefore overstates what the score can show.

### 4.2 Dependence on depth, read length and k

Clean lines (`run_sensitivity.tsv`):

* has rises from 0.88–0.95 at sensitivity ≤ 0.3 to 0.98–0.99 at ≥ 0.9.
* The null rises from about 0.03 to 0.04–0.05 per run and to 0.07 when runs
  are pooled. Error k-mers grow linearly with depth, while core recovery
  saturates.
* Read N50 (8–58 kb) has no visible effect.

GEMS-0061 and GEMS-0113 lie at 3–4× the null at every sensitivity from 0.19
to 0.94.

ONT null of foreign-allele scores (`ont_null_calibration.tsv`; 7 clean
genotypes, runs and halves; leave-tested-genotypes-out):

| k | median | p95 | max | pooled max |
|---|---|---|---|---|
| 21 | 0.034 | 0.047 | 0.056 | 0.073 |
| 31 | 0.020 | 0.028 | 0.036 | 0.046 |

Pooled scores at k = 21 / k = 31, with excess over the pooled ONT null
(`foreign_scores_by_k.tsv` gives every run and line):

| query | k = 21 | excess | k = 31 | excess |
|---|---|---|---|---|
| GEMS-0061, 8 lines | 0.150–0.196 | 0.077–0.123 | 0.122–0.170 | 0.076–0.124 |
| GEMS-0113, FBLL | 0.178 | 0.105 | 0.166 | 0.120 |
| GEMS-0113, other lines | 0.078–0.108 | 0.005–0.035 | 0.062–0.094 | 0.016–0.048 |

The null falls by about 40 % from k = 21 to k = 31, but the excess does
not. The classifications agree at both k, and so do the values above the
null.

### 4.3 In-silico admixture (ONT-appropriate positive control)

Real pooled ONT sketches were mixed at known fractions. The donor's
abundances were binomially thinned to fraction f of genome depth, and each
mixture was scored with the production code (`k21/spikein_admixture.tsv`):

| f | donor lacks (4 host/donor pairs) | production class | excess dosage |
|---|---|---|---|
| 0.01 | 0.074–0.121 | supported (1 ambiguous) | ≈ 0–0.05 |
| 0.02 | 0.102–0.166 | supported or ambiguous | ≈ 0–0.06 |
| 0.05 | 0.181–0.294 | ambiguous | 0.02–0.07 |
| 0.10 | 0.304–0.473 | ambiguous | 0.08–0.11 |
| 0.20 | 0.519–0.727 | ambiguous | 0.19–0.22 |
| 0.50 | 0.922–0.976 | ambiguous | 0.44–0.50 |

Admixture leaves the host's own alleles in place (has ≥ 0.91 even at
f = 0.5). The dosage estimate tracks f. The production classifier flags
admixture from about 2–5 %. That is good sensitivity, but no specificity
between contamination and genetic difference.

### 4.4 Detection limits (`detection_limits_k21.tsv`)

The smallest whole-genome admixture of panel line X detectable from X's
private k-mers is about 0.6–1.7 % for most pooled sets. It is 1.3–2.5 % for
GEMS-0113 (sensitivity 0.80) and 1.7–4.9 % for PHBW8 (sensitivity 0.22).
FBLL (1,010 private k-mers) and GEMS-0227 (2,827) are the hardest donors to
detect.

These limits are indicative only. They assume independent k-mers, which is
anti-conservative. In GEMS-0061 and GEMS-0113, contamination by a non-panel
line at a few percent would be hidden under the larger genetic signal.

### 4.5 Normalization at low sensitivity (GEMS-0219 run 9)

The normalization is unstable at low sensitivity. Run 9 used 26 reads
(175 kb), and its sensitivity of 4.8 × 10⁻⁵ rests on 5 of 103,408 core
k-mers. Its GEMS-0219 "has" scores are 1.06–2.62 (3–9 hits per pair). They
exceed 1, which is impossible for a fraction. Half A (1 core hit) names
GEMS-0113 as "best" with all-zero scores, which is an arbitrary tie.

* The Wilson interval covers the numerator only. The uncertainty in the
  sensitivity is ignored.
* `expected_z` uses the HiFi null spread without sampling error. Run 9 has
  the highest z of all queries (78.2), so z must not be quoted as
  confidence.

The classification was nevertheless correct (insufficient data).

### 4.6 Marker definitions

* **Pairwise distinct sets** overlap across genotypes, and their composition
  differs. One generic signal therefore flags 7–8 lines at once and inflates
  FBLL (§3). The line-specific excess test (same carrier count, with vs
  without H) separates source-specific from generic signal. A true source
  casts a smaller "shadow" excess on lines that share its alleles, so read
  the largest excess, not the smaller ones.
* **Private sets** are line-specific but small, and they depend on panel
  membership.
* **Multiplicity:** k-mers at 1.5–3× the peak make up 8–14 % of informative
  k-mers. They behave as expected (dosage above 1) and do not drive the
  signal.
* **Overlapping k-mers:** at 1/2000 sampling, sampled k-mers are about 2 kb
  apart, so overlap is negligible.
* **Linkage:** k-mers are not independent (haplotype blocks; read-level
  clustering 2.4×). All binomial intervals are therefore too narrow.

### 4.7 Independence of controls

* Leave-one-run-out controls are independent biological replicates (separate
  SMRT cells) only for the six genotypes with ≥ 2 HiFi runs.
* GEMS-0061, GEMS-0113, PHBW8 and PHP02 have only split-half controls (same
  DNA). These are technical controls only, so the two lines in question
  have no independent HiFi control.
* The negative controls include genuine within-line heterogeneity
  (GEMS-0219 bc2041; GEMS-0227 bc2004, FBLL 0.096). They are a mixed null,
  not an error null.
* No ONT data set has independently known truth.

**What the controls do and do not support:**

* They support "the reads are this genotype and no other panel genotype".
* By themselves they do not support "no secondary genotype".
* The ONT null and spike-in calibration above supply that, within the
  stated limits.

## 5. GEMS-0219 run 9 (read-only verification)

* **Insufficient data: confirmed.** The run has 52 reads in total, 26 of
  them ≥ 1 kb (175 kb); 5 core k-mers were found; GEMS-0219-distinct hits
  are 3–9 per pair, and foreign hits are 0. The few hits are compatible
  with GEMS-0219 but carry no weight.
* **9b and 9c support GEMS-0219.** Both are supported at k = 21 and k = 31,
  in all halves:

  | run | has | lacks | sensitivity |
  |---|---|---|---|
  | 9b | 0.906 | 0.043 | 0.89 |
  | 9c | 0.878 | 0.034 | 0.69 |

  There is no foreign gradient or excess.
* **Caveat for 9b and 9c.** has is lower than in other clean lines, and
  recovery of the HiFi's low-dose alleles is 0.82. This reflects mild
  heterogeneity among GEMS-0219's HiFi runs (bc2041), not the ONT.
* **Not in the assembly workflow.** `config/samples.yaml` lists run 9 under
  `exclude_inputs`. That file is unchanged since commit ec639e6, and the
  working tree was clean before this work. The launch-time dry run
  (`logs/assembly_20261012_023650.dryrun/GEMS-0219.dryrun.log`) has exactly
  two `ont_filter` jobs (9b, 9c) and no reference to run 9. GEMS-0219 has
  not started yet: running now are 3AZA1, FBLL and GEMS-0227. It will read
  the same `samples.yaml` when it starts. Run 9 also has no read ≥ 50 kb, so
  it could not contribute UL reads in any case.

## 6. Genotype QC summary

"has" and "lacks" are k = 21 values for runs, with pooled values in
parentheses. "Excess" means foreign alleles above the ONT null.

| genotype | ONT runs | identity evidence | secondary-genotype evidence | limitations | recommendation |
|---|---|---|---|---|---|
| 3AZA1 | 1 | strong: has 0.998, lacks 0.050 | none detected (no gradient or excess; detection about 1 %) | single ONT run | **supported** |
| FBLL | 2 | strong: 0.96–0.99 | none | few private k-mers (1,010): FBLL is the hardest donor to detect | **supported** |
| GEMS-0061 | 3 (2 flow cells) | strong: best in all runs and halves; next-best net ≤ −0.54 | foreign alleles real, at genomic dosage, inconsistent with low-level admixture; genetic discordance between HiFi and ONT sources (16 % own alleles absent; lacks 0.15–0.20) | HiFi heterozygous or heterogeneous; no independent HiFi control; cause unknown | **supported with caveats** |
| GEMS-0113 | 2 (2 flow cells) | strong: next-best net ≤ −0.65 | FBLL-like haplotype segments plus a generic signal in the ONT source; dosage (0.45–0.8) argues against FBLL admixture; genetic discordance between HiFi and ONT sources | as GEMS-0061; pedigree unverified | **supported with caveats** |
| GEMS-0219 | 9b, 9c (9 excluded) | strong: 0.88–0.91 (0.92) | none in ONT | mild heterogeneity among HiFi runs (bc2041); run 9 insufficient | **supported**; run 9 **exclude** (done) |
| GEMS-0227 | 3 | strong: 0.91–0.97 (0.99) | none | HiFi run bc2004 differs slightly from bc2056 | **supported** |
| PHB47 | 3 | strong: 0.92–0.98 | none | — | **supported** |
| PHBW8 | 1 (pre-filtered ≥ 50 kb) | strong: 0.94; sensitivity 0.22 | none | shallow; detection 2–5 %; single HiFi run | **supported** |
| PHP02 | 1 | strong: 0.997 | none | single HiFi run | **supported** |

**No ONT run needs to be excluded on fingerprinting grounds.** Run 9 is
already excluded for insufficient data.

**Assembly inputs.** The current inputs are correctly labelled. For GEMS-0061
and GEMS-0113, however, the HiFi+ONT assemblies combine two genetically
non-identical DNA sources. Where the sources differ, UL reads carry a
different haplotype from the HiFi graph. This does not by itself invalidate
the assemblies (D1 and the 2026-10-12 decision stand), but it must be stated
in the manuscript. It is checked once the assemblies have finished, by
comparing GEMS-0061 and GEMS-0113 with the other seven lines
(`docs/post_assembly_qc_plan.md`). HiFi-only comparison assemblies are not
planned unless that QC justifies them.

## 7. Recommendations (need approval; nothing implemented)

1. Report identity as **supported for all nine genotypes**. Report
   GEMS-0061 and GEMS-0113 as showing **genetic discordance between the HiFi
   and ONT sources**. Contamination is not supported by these data but is
   not categorically excluded, and the biological provenance is not
   established. (Adopted by the project lead, 2026-10-12.)
2. Obtain seed-lot and plant provenance for the HiFi and ONT DNA of GEMS-0061
   and GEMS-0113. Where available, compare both data sets with external
   genotype data for these lines to establish which source represents the
   registered line.
3. Fingerprinting module (for a later, versioned change):
   * separate the identity call from the secondary-signal call;
   * use an ONT-matched null (leave-genotype-out) for the secondary signal;
   * replace "possible_admixture" with a reported dosage and line-specific
     excess (spike-in-calibrated), or with wording that says only that alleles
     are present;
   * refuse to rank when core hits < 1,000;
   * propagate uncertainty in the sensitivity;
   * drop or replace `expected_z`;
   * document panel-composition effects (central lines such as FBLL).

   This would make the module describable in a genome-resource paper as an
   identity test with a calibrated detection limit. The current thresholds
   were not changed.
4. A read-level spike-in (real reads, not sketches) to confirm the dosage
   calibration. **Deferred** (project lead, 2026-10-12): only if
   post-assembly QC justifies it. It writes temporary mixtures and needs
   approval.

## 8. Reproduction and resources

```bash
cd pipeline
PY=/programs/Python-3.12.7/bin/python3.12
R=../results/fingerprint/report-742a05a563; D=../results/fingerprint/diagnostics-2026-10-12
$PY scripts/g2f_fp_diagnostics.py --report $R --outdir $D --k 21      # ~3 min, 3.3 GB
$PY scripts/g2f_fp_diagnostics.py --report $R --outdir $D --k 31      # ~4 min, 5.2 GB
$PY scripts/g2f_fp_diagnostics.py --report $R --outdir $D --k 21 \
    --spikein PHB47:FBLL GEMS-0227:FBLL GEMS-0219:FBLL PHB47:PHP02    # ~3 min
SM="env PYTHONPATH=/programs/sourmash-4.9.3/lib/python3.12/site-packages $PY"
$SM scripts/g2f_fp_diagnostics.py --report $R --outdir $D --k 21 --only-reads \
    --reads-max-bases 1500000000 --reads <tar>:<LABEL> ...            # ~3 min per 1.5 Gb
$PY scripts/g2f_fp_diagnostics_summary.py --report $R --diagnostics $D
```

* Tables in the repository (`docs/fingerprinting_diagnostics_2026-10-12/`)
  are copies of the files above; the per-read tables (`reads_*.tsv`, 7.6 MB)
  and run logs stay in `results/` only and are listed with checksums in
  `MANIFEST.sha256`.
* Read scans covered GEMS-0061 run 8 (complete), GEMS-0113 run 7b (first
  1.5 Gb), PHB47 run 4c (complete) and PHP02 (first 1.5 Gb). The tars were
  streamed read-only and nothing was copied.
* All jobs ran single-threaded at `nice 19`, with peak RSS 5.2 GB and about
  35 min of CPU in total. Running assemblies were not touched.

Dosage values for excesses below about 500 k-mers are imprecise, and are
reported as NA below 100.
