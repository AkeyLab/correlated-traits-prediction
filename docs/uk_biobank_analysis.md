# UK Biobank type 2 diabetes analysis — not in this repository

## What is missing

The paper's empirical section applies the framework to type 2 diabetes in the UK Biobank.
None of that analysis is in this repository. Specifically, the following are **not** here:

| Manuscript element | What produced it |
|---|---|
| Genome-wide association study of type 2 diabetes (11,411 cases, 132,718 controls; 37 genome-wide significant variants) | REGENIE, run on UK Biobank genotypes |
| Baseline polygenic model, AUC-ROC 0.677 | Regularized logistic regression on those 37 variants plus covariates |
| Candidate helper trait assembly (251 NMR metabolomic traits, 19 clinical and lifestyle traits) | UK Biobank phenotype extraction |
| Correlation structure of helper traits (Figure 5A–C) | Pairwise correlations among candidate helper traits |
| Redundancy-aware "tag helper trait" selection | The greedy algorithm written out in Supplementary Methods 1 |
| Incremental AUC-ROC curve (Figure 5D), reaching 0.907 | AutoGluon, sequentially adding tag helper traits |
| HbA1c reference model, AUC-ROC 0.889 | Logistic regression on covariates plus HbA1c |

The selection algorithm itself is fully specified as pseudocode in Supplementary Methods 1 of
the paper, so it is reimplementable from the text; the surrounding pipeline is not.

## Where it lives

The empirical analysis runs on the **UK Biobank Research Analysis Platform (DNAnexus)**, not
on local infrastructure. The code sits inside that controlled-access environment, which is
why it does not appear in this repository or anywhere in the authors' local research tree.

Releasing it means exporting the analysis code — and only the code — off the platform, then
confirming that nothing exported embeds participant-level values, cohort identifiers, or
intermediate per-individual files. Export of code is permitted; export of data is not.

## Why it is not here

The analysis reads **individual-level UK Biobank data**, which is controlled-access. This
repository is intentionally free of any code path that touches participant-level records, so
that it can be published openly without a data-governance review.

## Getting the data

Individual-level UK Biobank data are available to approved researchers by application:
<https://www.ukbiobank.ac.uk/enable-your-research/apply-for-access>. The analyses in the paper
were conducted under **application number 104628**.

## Software the empirical analysis used

| Tool | Purpose | URL |
|---|---|---|
| REGENIE | Genome-wide association study | <https://rgcgithub.github.io/regenie/> |
| PLINK 2.0 | Genotype quality control and handling | <https://www.cog-genomics.org/plink/2.0/> |
| AutoGluon | Automated model selection for the predictive models | <https://auto.gluon.ai/> |

## If this pipeline is added later

Genome Biology requires that analysis code be publicly available, and holds papers in
production until it is — so this is a submission blocker, not an optional extra. The code
being on the Research Analysis Platform does not exempt it; it has to be exported and
published.

If the UK Biobank pipeline is released, it should live in a `ukbiobank/` directory here with:

- extraction scripts that take **field IDs**, never participant identifiers;
- the REGENIE and PLINK invocations, with their exact parameters;
- the tag-helper-trait selection implementation matching Supplementary Methods 1;
- the model-fitting and evaluation code that produces Figure 5;
- a statement of which steps require controlled-access data and therefore cannot be run
  by a reader without their own approved application.

No participant-level data, derived individual-level files, or extracts should ever be
committed, regardless of how aggregated they appear.
