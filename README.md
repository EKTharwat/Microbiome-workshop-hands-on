# Microbiome hands-on workshop

This version contains the supplied supervisor side-by-side material. Open [the commands page](workflow-code-side-by-side.html) or [the results page](workflow-results-side-by-side.html). The complete expanded mothur runs are available as [private release downloads](RUN_ARCHIVES.md).

This repository contains 30 pseudonymized paired-read samples. Each of five study groups, `G1` through `G5`, contains three `C1` and three `C2` samples. The sample identifiers are neutral and are used consistently in filenames, metadata, analyses and plots.

The reads are a selected bacterial-panel dataset. They were screened and matched against a 16S reference panel with 100% coverage of the assembled read. Four study groups use a minimum of 99% identity; one uses 97%. The identity threshold for each sample is recorded in the metadata. Ambiguous assignments and candidate chimeras were removed. Samples were selected partly because they retained panel matches. The different identity thresholds affect which reads are retained, so do not interpret differences between study groups as purely biological. A 16S match is a putative association and may not uniquely identify a species. Diversity and differential abundance results describe these selected reads and are exploratory.

For the filtering procedure and its limitations, see [FILTERING.md](FILTERING.md).

For the tested software versions and package checklist, open [prerequisites.docx](prerequisites.docx). For step-by-step instructions, open [beginner-guide.docx](beginner-guide.docx).

## Start the workshop in RStudio

1. Install R and RStudio Desktop, then download and extract the complete workshop folder.
2. In RStudio choose **File → New Project → Existing Directory** and select the folder containing `workflow.Rmd`, `scripts`, `metadata`, and `vendor`.
3. Paste this line into the RStudio Console and press Enter. It installs packages, detects your operating system, selects mothur and prepares the data:

```r
source("scripts/prepare_packages.R")
```

Keep RStudio open while missing packages are installed and the data are prepared. The first installation requires internet access. The package installer uses a separate R process, avoiding package updates inside the active RStudio session. If RStudio restarts, run the line above again; do not run an internal line from the installer by itself.

4. Open `workflow.Rmd` in RStudio. Run its commands from the top in order using the green Run button. Each mothur block runs one mothur command and streams its output; copy and run the entire block in the RStudio Console. Other R steps remain individual commands. Earlier commands create objects used by later ones.
5. The notebook prepares the bundled data and opens one mothur session from R. Each mothur command is sent immediately; the same block waits for completion and continuously prints new output in the Console. Run all the steps in order so mothur keeps the correct current files.
6. Select **Knit** in RStudio to create a fresh `workflow.html` report with code, tables and plots after their corresponding sections. To regenerate the code-only HTML in the Console, run `source("scripts/render_code_only.R")`.

The workflow detects Windows, macOS, or Linux automatically. See [beginner-guide.docx](beginner-guide.docx) for prerequisites and the complete walkthrough. To run the entire workflow automatically in a fresh folder instead, run `source("scripts/start_workshop.R")` in the RStudio Console.

## Learn the mothur steps in RStudio

Open `workflow.Rmd` and run each command in order. Mothur commands run immediately in the shared session; the following R steps wait for each command and print its output. Every processing command names its input files explicitly. The diagnostic `get.current` command reports session state and marks command completion; it is not used to select input files. The workflow follows the MiSeq SOP sequence: assembly and screening, reference trimming and alignment, filtering and classification, and feature and diversity tables. The R sections then follow four steps: (1) quality control, with reads per sample and MetaQC (IQC, EQC, CQCg, CQCp, AQCg, AQCp) to decide which studies can be combined; (2) alpha diversity (Observed, Chao1, Shannon, Inverse Simpson) and beta diversity (Bray-Curtis, Canberra, Jclass, ThetaYC) with PCoA, PERMANOVA and PERMDISP on the included studies; (3) ANCOM-BC run separately in each included study, keeping the effect size (log fold change) and standard error of every OTU; (4) a random-effects meta-analysis of those effect sizes with metafor::rma(method = "REML"), followed by one cladogram for all included studies and the pooled effect on alpha diversity.

Each new run and each Knit creates a fresh folder under `runs/`; ensure you have enough disk space. The bundled pages are `workflow-results-side-by-side.html` (every command with its output and figures) and `workflow-code-side-by-side.html` (the same commands with explanations and a Copy button). Their mothur outputs come from the completed test run; their R outputs and figures were produced from the mothur tables in `results/mothur` and are saved in `report_outputs/by_study`. Knit the current notebook to generate a matching report from your own run. To rerun only the R analysis from the finished mothur tables, run `source("scripts/downstream_by_study.R")`.

All packages are loaded in one chunk at the start of the R part. Figure colours and theme follow Marzouk et al., npj Parkinson's Disease 2025 (https://www.nature.com/articles/s41531-025-00950-z).

Each study has only three samples per condition. The studies differ in read-filter identity threshold, so pooled and study-specific results should be interpreted with care. See [FILTERING.md](FILTERING.md) for how the supplied reads were selected.
