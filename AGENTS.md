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
- Use Git for all scripts, configurations, and documentation.
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
- Write small tests for new scripts.
- Document scientific decisions and their justification.
- Avoid refactoring validated analyses unnecessarily.

## Project organization
- config/: samples, parameters, and environment settings
- scripts/: version-controlled executable utilities
- workflow/: reproducible workflow definitions
- tests/: lightweight testing
- docs/: scientific decisions, benchmark data, provenance

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
