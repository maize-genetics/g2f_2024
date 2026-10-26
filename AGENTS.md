# G2F Maize Genome Resource — Agent Instructions

## Project mission
Develop a scientifically defensible, reproducible genome
assembly and annotation pipeline for nine maize G2F lines.

Immediate target: assembly QC and manuscript-ready resources.
Target preprint: November 15, 2026.

Long-term goal: a reusable, species-agnostic genome resource
pipeline with modular assembly and annotation components.

## Scientific principles
- Treat PHBW8 as the initial validation sample.
- Use hifiasm v0.25.0 as the initial assembly baseline.
- Evaluate -l0 (no haplotig purging) for inbred maize.
- Do not assume all samples have identical heterozygosity.
- Do not change scientific parameters between samples silently.
- Keep HiFi-only and HiFi+ONT methods distinguishable.
- Preserve all original and historical assembly results.
- Do not interpret BUSCO duplication alone as haplotig
  contamination; maize has a history of genome duplication.
- Verify claims against actual software output and logs.
- Never invent software versions, sample metadata, or results.

## Agent execution policy
The human researcher controls all substantial computation.

You MAY:
- Inspect files and directories without modifying them.
- Read logs, scripts, metadata, and configuration files.
- Create or edit project code and documentation.
- Run lightweight, safe validation tests.
- Perform read-only environment inspections.

You MUST NOT execute without explicit approval:
- Genome assembly or other long-running analyses.
- High-CPU, high-RAM, or high-I/O operations.
- Batch jobs, background jobs, or scheduler submissions.
- Large downloads, installations, or data transfers.
- Data deletion, renaming, moving, or overwriting.
- Operations affecting files outside the project.
- Git pushes or changes to remote repositories.

For substantial jobs:
1. Generate an executable script or exact command.
2. Explain CPU, RAM, storage, and expected outputs.
3. Tell the researcher how to start and monitor the job.
4. Let the researcher execute it independently.
5. Inspect completed logs and outputs when requested.

Never watch long-running jobs through repeated tool calls.
Do not waste context polling or tailing logs.

## Reproducibility
- Use Git for reusable code, portable configuration and user
  documentation only (see Repository hygiene).
- Keep sequencing data and large outputs outside Git.
- Use one sample manifest as the authoritative input registry.
- Separate scientific parameters from resource settings.
- Record exact commands, paths, tool versions, and inputs.
- Prefer explicit named outputs and structured log files.
- Make scripts safe to rerun without overwriting results.
- Include dry-run modes wherever practical.
- Track states: planned, running, completed, failed,
  and validated. Never equate file existence with success.
- Make the pipeline executable without an AI agent.

## Development standards
- Prefer simple, auditable solutions over complexity.
- Separate reusable tools from maize-specific configuration.
- Use Bash for simple wrappers and Python for structured logic.
- Introduce a workflow manager when it provides clear value.
- Include input validation and useful error messages.
- Write small tests for new scripts, using synthetic data only.
- Document scientific decisions and their justification.
- Avoid refactoring validated analyses unnecessarily.

## Project organization
- config/: portable parameters, defaults and *.example.yaml
  templates; project-specific files are untracked
- config/local/: untracked machine and project settings
- config/workflow_stages.yaml: machine-readable stage index
- scripts/: version-controlled executable utilities
- workflow/: reproducible workflow definitions
- profiles/: optional execution profiles (e.g. BioHPC)
- tests/: lightweight tests on synthetic data
- docs/: user documentation for the reusable pipeline
- ../project_records/ (outside the repository, never in Git):
  scientific decisions, run plans, audits, incidents, QC
  reports, fingerprinting results, provenance notes

## Discovering and running the workflow
- Start with docs/workflow_guide.md: what each stage of the process does (inputs and identity → reads → k-mers →
  assembly → QC → structural validation → correction → selection → release), its entry point, and its gate.
- config/workflow_stages.yaml is the machine-readable index. Each stage has a `status`:
  - `workflow` / `tool`: code in this repository;
  - `project_script`: executed in a project, not yet packaged here;
  - `manual`: a human decision;
  - `planned`: not implemented.

  Do not present a `project_script` or `manual` stage as automated.
- Snakemake stages: always `--dry-run --printshellcmds` first, and show the researcher the plan. CLI tools: run
  `--help`, then `plan` / `check-env` / `status` where they exist, before any `run`.
- To learn what a past result came from, read its provenance JSON (schema `g2f-provenance/1`), or tabulate many with
  `scripts/g2f_provenance_report.py`. For tools with staged outputs, read `VALIDATED.json`.
- A project's history, decisions and evidence live in `../project_records/`. Read its workflow record before
  re-deriving anything.

## Validation expectations
- **A result exists only when its validation marker exists.** That is `*.validated.*.json` (workflow) or
  `VALIDATED.json` (tools), not the output file.
- **Every new script ships with a synthetic test** under tests/ and is registered in config/workflow_stages.yaml.
  `tests/test_workflow_stages.py` checks that the registry matches the code.
- **Run `bash tests/run_tests.sh`** (synthetic data, a few minutes) before committing code changes.
- **A derived file must be verified independently of the code that made it.** Derived files include a renamed,
  reordered or split FASTA. The verification is a separate pass that re-reads the source and the output. It compares
  per-record sequence hashes, IDs (each exactly once), lengths and totals.

## Release rules
- **Releases are immutable.** Never edit, regenerate or delete a file in an existing release directory.
  - Any change, including names or record order, goes into a new versioned directory.
  - Its README states exactly what changed and what did not.
- **Before and after building a derived release,** re-verify the source release (its checksum file and a tree
  snapshot) and record the result inside the new release.
- **Every release has an authoritative manifest** (one FASTA per unit, with path, sha256, record count, bp and
  source) and a checksum file. Make the files read-only after validation.
- **Identifier and order conventions** do not imply scaffolding. Contigs that may be alternative alleles are flagged
  in metadata and never silently removed. See docs/workflow_guide.md, stage 8.

## Decisions that require a human
Agents prepare evidence tables. They never make these decisions or encode them as defaults:
- **Assembly strategy for a heterozygous or mixed sample:** purge level, haplotype pair vs primary, coverage overrides,
  and which reads represent the sample.
- **Structural edits:** whether to break or split a join. The two cases are reads CONTRADICT the join, or the join is
  read-untestable at a known artefact.
- **Selection of a deliverable when QC metrics conflict**, and reporting an unexplained check failure as unexplained.
- **Promotion of a build to a release**, release approval, and anything sent to collaborators.

Keep three things separate:
- technical validation (checksums, file contracts);
- scientific QC (is the assembly good?);
- collaborator approval.

## Repository hygiene (Git policy)
GitHub distributes reusable software. It does not document
research operations or Cornell infrastructure.

Git may contain ONLY:
- Pipeline source code and workflow definitions.
- Portable configuration and *.example templates.
- Container definitions and software environments.
- Synthetic tests and fixtures (invented values and paths).
- README, installation and user documentation.
- Agent instructions (AGENTS.md, CLAUDE.md).
- Optional execution profiles without private details.

Git must NOT contain:
- Execution logs, diagnostic dumps, process histories.
- Incident or reboot reports, SSH troubleshooting, run
  plans, backup reviews, or other operational records.
- Internal IP addresses, unnecessary hostnames, personal
  usernames or paths, backup or storage architecture.
- Scientific results: assembly or QC reports, BUSCO or
  fingerprinting outputs, metrics intended for publication.
- Reads, assemblies, intermediates, or real sample
  metadata (samples.yaml, manifests, run/flow-cell IDs).
- Private AI notes or session transcripts.
- Credentials, tokens, secrets, or personal information.

Rules for agents:
- Write operational and scientific records to
  ../project_records/, never to docs/ or other tracked paths.
- Use generic names (SAMPLE_A, /path/to/...) in examples,
  tests and documentation; never real sample IDs, hosts,
  usernames or results.
- Code must not require personal paths or Cornell/BioHPC
  infrastructure; BioHPC settings stay optional
  (config/tools/biohpc.yaml, profiles/biohpc).
- Public software installation paths (e.g. /programs/...)
  are acceptable in optional BioHPC files only.
- Before every commit, review the staged diff for the items
  above. Before every push, scan every object in
  `git rev-list --objects origin/main..HEAD`, including
  commit messages; earlier versions in unpushed history
  count too. Commit author email: a GitHub noreply address.
- Never move private material out of Git by deleting it:
  copy it to ../project_records/, verify checksums, then
  `git rm --cached`.

## Manuscript requirements
Every reported numerical result must be traceable to:
- a specific sample
- an input dataset
- a software version
- a parameter configuration
- a documented execution
- a validated output

Generate summary tables from verified results rather
than manually transcribing metrics whenever possible.

## Scope boundaries
Priority 1: validated assemblies for nine G2F maize lines.
Priority 2: assembly QC, comparison, and publication tables.
Later options: EDTA, geneCAD, PlantCAD, pan-genome visualization.

Do not incorporate an unvalidated
experimental assembler into the production pipeline.

Follow Cornell BioHPC policies in addition to these rules.
