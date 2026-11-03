# G3 requirements for Genome Reports

Pages consulted 2026-11-02:
- https://academic.oup.com/g3journal/pages/author-guidelines
- https://academic.oup.com/g3journal/pages/data-policy
- https://academic.oup.com/g3journal/pages/About
- Official G3 Overleaf template (states "Updated: 1 April 2022")

**How far each rule was confirmed**

| label | meaning |
|---|---|
| PAGE | wording read directly from the journal page on the date above |
| AUTHOR-CHECKED | an author read the rule on the current official guidelines and reported it, on the date shown. Automated reading of the page still stops before "Preparing your manuscript" |
| AUTHOR-QUOTED | pasted from the guidelines by the authors into their working notes |
| TEMPLATE-2022 | wording in the official Overleaf template, dated 2022; lower authority than the current journal website |
| CHECK | not retrieved; read the page before submission |

## Article type

- Reports include Genome Reports, Genetic Screen Reports, Molecular Method Reports, Omics Reports, and Software and Workflow Reports. (PAGE)
- A Genome Report publishes high-quality assembled and annotated whole-genome sequence. It generally includes new genomic data and an analysis of intraspecific variation and/or comparative genomics. Its contribution must be clear relative to the existing literature. (AUTHOR-QUOTED)
- It should be a useful resource for its community, give a rationale for a general readership, and be clear to a non-specialist. (AUTHOR-QUOTED)

## Length (Reports)

| section | limit | status |
|---|---|---|
| Abstract | 150 words at most | AUTHOR-CHECKED 2026-11-02 |
| Introduction | 500 words at most | AUTHOR-CHECKED 2026-11-02 |
| Discussion | 500 words at most | AUTHOR-CHECKED 2026-11-02 |
| Materials and Methods | no word limit | AUTHOR-CHECKED 2026-11-02 |
| Article summary | 100 words | AUTHOR-CHECKED 2026-11-02 |

The current journal website takes precedence over the 2022 Overleaf template, which gives an older abstract rule: "a single paragraph of less than 250 words" that does "NOT contain citations or abbreviations" (TEMPLATE-2022). Do not apply the template's word limit. Other template wording is a suggestion unless the website confirms it.

## Figures

- Main figures require alt text. (AUTHOR-CHECKED 2026-11-02)
- Panel labels are lowercase. (AUTHOR-CHECKED 2026-11-02)

Template section order, not stated as a requirement: Abstract; Introduction; Materials and methods; Results and discussion; Data availability; Acknowledgments; Funding; Conflicts of interest; bibliography (TEMPLATE-2022). The template has a keywords field and sets no figure or table count.

Figure and table count limits, required section order, keywords rules, reference style, supplementary file rules, and CRediT, funding and author conflict statements: CHECK.

## Questions reviewers ask (AUTHOR-QUOTED)

1. Stand-alone assembly, or integrated with a previous one? If integrated, are the data described, referenced and accessible?
2. How were individuals selected and sampled? Are sample characteristics described (sex, tissue, genotype)? Does sampling limit usefulness as a reference?
3. Is quality control sufficient? Are raw data quality and scale reported? Is assembly size reasonable against expected genome size? Are N50, L50, NG50, NGA50 and percentage of Ns acceptable?
4. Are protein-coding genes, transposable elements and other repeats annotated, in enough detail to reproduce? Is a conserved-gene measure such as BUSCO reasonable?
5. Are tools up to date and appropriate? Is the approach detailed and accompanied by scripts?
6. Are other analyses (structural variation, repeats) described in detail, with scripts?

## Data and software (PAGE)

- All data and software underlying a paper must be publicly released as a condition of publication.
- "Journal policy does not allow for data to be available upon request."
- Raw sequence reads must be in a public repository. Gene sequences go to GenBank, EMBL-Bank, DDBJ or SRA. Genotypes or sequences for all individuals should be in common formats such as VCF, alignment or plain text.
- General repositories include figshare, Dryad and Zenodo; repositories that issue DOIs are encouraged. G3 figshare portal: https://gsajournals.figshare.com/submit
- New software must be in a third-party public repository and licensed; an open-source license is strongly encouraged.
- A Data Availability Statement is required at the end of the manuscript. It lists accessions or DOIs and the names and descriptions of supplemental files.
- If accessions are not yet available, state where data will be submitted. They must exist before publication, not at first submission.
- Noncompliance may lead to return without review.

The pages read do not state assembly- or annotation-specific deposition rules (for example, for genome annotation files): CHECK.

## AI use (PAGE)

- AI tools cannot be authors.
- Use of AI (for example, to help generate content or images, write code, process data, or translate) should be disclosed both in the cover letter and in the Methods or Acknowledgements.

## Other points (PAGE)

- A corresponding author must be named in the system and in the author list. ORCID iDs are encouraged.
- Peer review is single-anonymized. Articles are CC BY 4.0.
- Preprints are allowed; after acceptance, update the preprint with the published DOI.
- An editor or board member among the authors must be declared.

## Pre-submission checklist

- [ ] Abstract, Introduction, Discussion and article summary are within limits, by tool count
- [ ] Every main figure has alt text; panel labels are lowercase
- [ ] Sampling and sample characteristics described
- [ ] Raw data scale and quality reported per sample
- [ ] Size vs expected, N50, L50, NG50, % Ns (and NGA50 if reported)
- [ ] BUSCO lineage and version; consensus accuracy method
- [ ] Gene, TE and other repeat annotation, reproducibly described
- [ ] Intraspecific or comparative analysis present
- [ ] Contribution stated relative to the existing literature
- [ ] Reads, assemblies, annotations and code deposited or destinations stated; none "on request"
- [ ] Code in a third-party repository, licensed, with a persistent identifier
- [ ] Data Availability Statement at the end, listing supplemental files
- [ ] AI use disclosed in the cover letter and in the Methods or Acknowledgements
- [ ] Items marked CHECK above confirmed on the live page
