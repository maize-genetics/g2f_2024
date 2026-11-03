---
name: g3-genome-report
description: Audit, draft and review G3 Genes|Genomes|Genetics Genome Reports against verified scientific evidence. Use whenever someone is writing, checking or revising a genome-assembly or annotation manuscript for G3, comparing a manuscript with project records, drafting evidence-grounded Methods or Results, preparing a provenance table, checking word limits or the data availability statement, or preparing the AI-use disclosure, even if they only say "genome report", "check my manuscript" or "G3 submission".
argument-hint: "[audit|draft|review] [section or file]"
---

# G3 Genome Report

Three modes share one rule: **a statement enters a manuscript only with its evidence and its status.** A confident sentence with no traceable source is worse than a visible placeholder, because co-authors and reviewers cannot tell it apart from verified work.

Before starting, read:
1. `references/g3_requirements.md`: journal requirements and how current they are.
2. `references/evidence_rules.md`: status labels, reporting traps and the provenance table format.
3. Project context, if present: `manuscript_work/PROJECT_CONTEXT.md` or the path the user names. Generic journal rules stay in this skill. Sample names, paths and results stay in the project file, so the skill can be shared without unpublished material.

If the journal's page may have changed, re-read it. If it cannot be read, say which rule was not confirmed. Do not fill the gap from memory.

## Mode: audit (manuscript against evidence)

Goal: tell the authors what in the manuscript is wrong, unsupported or missing, with evidence.

1. Treat the manuscript as read-only. Record its checksum and date.
2. Find the authoritative evidence: frozen release manifest, validation markers, provenance records, then decision records. Check the release's actual contents; a README can lag behind its data.
3. Extract each quantitative claim, method claim and resource claim. Classify each as supported, outdated, incorrect, unsupported or conflicting.
4. Recompute what is cheap to recompute (for example, contig N50 from a `.fai`) instead of copying a summary.
5. List missing analyses and journal requirements.
6. Write findings ranked by consequence: claim, what the evidence says, source file, recommended correction. Name who must decide when the fix is a judgment call.

## Mode: draft (evidence-grounded writing)

Goal: text a co-author can paste in, with gaps that cannot be overlooked.

1. Update the provenance table first, then write from it.
2. Build tables from source files with a script, and check them independently. Do not transcribe metrics by hand.
3. Write plain scientific prose for a non-specialist geneticist. Report what was done, with versions and parameters, in the past tense.
4. Mark gaps inline: `[MISSING: ...]`, `[COLLABORATOR: ...]`, `[PROVISIONAL: ...]`. Never bridge them with plausible wording.
5. Keep exploratory or rejected work out unless it explains a decision or limitation. Do not present a method that was only investigated as part of the released workflow.
6. Count words for every limited section with a tool, and report the counts.
7. Give each section a one-line status so its readiness is visible at a glance.

## Mode: review (accuracy, clarity, consistency, compliance)

Run these passes in order, because a polished sentence that is false has no value:
1. **Accuracy:** do numbers, methods and claims match the provenance table?
2. **Inference:** does each conclusion follow? Look for reference-guided structure reported as validation, and identity or completeness claims stronger than the test performed.
3. **Consistency:** do the abstract, text, tables, figures and terminology agree?
4. **Clarity:** can a geneticist outside the subfield follow it? Define project-specific terms.
5. **Compliance:** use the checklist in `references/g3_requirements.md`.

Report findings ranked by consequence, with location and suggested wording.

## Boundaries

- Do not invent numbers, citations, accessions, methods, author roles or biological conclusions. Check each citation against a bibliographic source, and mark any you could not check.
- Do not make scientific decisions for the authors. Prepare the evidence and ask: assembly strategy, structural edits, selection, release, authorship.
- Do not edit a shared manuscript or frozen data unless asked. Deliver text for human review.
- Do not start expensive computation. Propose it with the command and resource estimate.
- Keep unpublished manuscripts, review comments and AI-assistance logs out of public repositories.
- Name no individual collaborators in audits or decision documents unless attribution is required; refer to roles and topics.
- Keep an AI-assistance record (application, model, dates, what it was used for), because G3 asks for disclosure.
