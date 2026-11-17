# What changed in this version

## Update: MetaQC, per-study ANCOM-BC and REML meta-analysis

The downstream analysis now follows this order:

1. **OTU/ASV tables** from mothur (unchanged).
2. **MetaQC quality control.** Before any comparison, MetaQC (github.com/donkang75/MetaQC, Kang et al. 2012) scores every study on IQC (homogeneity of the co-abundance structure), EQC, CQCg/CQCp (consistency) and AQCg/AQCp (accuracy). OTUs are CLR-transformed per study and the taxonomic groups (phylum to genus) are used as MetaQC's "pathways". A study is excluded when it fails both IQC and CQCg (score below -log10(0.05/5) = 2.0). Result: G4 fails both (IQC 1.56) and is excluded; G1, G2, G3 and G5 are kept.
3. **Alpha and beta diversity** on the four included studies.
4. **ANCOM-BC per study**, never on a merged table; the effect size (log fold change) and standard error of every OTU are kept from each study.
5. **Random-effects meta-analysis** with `metafor::rma(..., method = "REML")`. Worked examples are now Otu05 and Otu16. Otu16 Streptococcus is the only OTU significant after adjustment (pooled log fold change +1.21, 95% CI 0.34 to 2.09, adjusted p = 0.046).

Technical notes:
- MetaQC is no longer on Bioconductor. Its GitHub source (version 0.1.14, GPL-2) is bundled as `vendor/MetaQC_0.1.14.tar.gz`; `scripts/install_packages.R` installs it together with its dependencies proto, foreach and iterators.
- MetaQC's IQC normally repeats its calculation 100 times on a random 10% of the shared genes. With 10 to 17 OTUs per study that would be a single OTU, so the workflow replaces that one internal function with the same calculation done once on all shared OTUs. Everything else in MetaQC is unchanged.
- With three samples per condition no study can pass EQC, CQC or AQC, so in practice the decision rests on IQC.
- The MetaQC inputs and scores are saved in `report_outputs/by_study/` (`metaqc_scores.csv`, `metaqc_decision.csv`, `metaqc_scores.png`, and the taxonomic group file in `metaqc/`).

## Earlier changes

Open `workflow-results-side-by-side.html` (results) or `workflow-code-side-by-side.html` (commands only). Word copies of both are in `Word_Documents/`.

- **Side by side only.** The separate analysis of each study (G1 to G5) was removed. Every figure now shows the five studies together with the two conditions inside each study.
- **Working directory.** The first command of both pages and of `workflow.Rmd` sets the working directory automatically: it uses the folder of `workflow.Rmd`, or asks you to select that file if R cannot tell where it is.
- **One library cell.** All packages are loaded in a single cell; no later step calls `library()`.
- **Alpha diversity.** Chao1 added to Observed, Shannon and Inverse Simpson; drawn as violin plots with a box plot inside.
- **Beta diversity.** Bray-Curtis, Canberra, Jclass and ThetaYC. Each PCoA has an ellipse per condition and the percentage of variation on the x- and y-axis. PERMANOVA (`adonis2`) and PERMDISP are run for every distance and study; with six samples vegan enumerates every permutation, so the tests are exact. The hand-made permutation matrix was removed because supplying a matrix to `adonis2` fails with vegan 2.7-6 / permute 0.10-0 (R 4.6).
- **ANCOM-BC replaces LEfSe.** `ancombc2()` gives each OTU an effect size (log fold change) with a standard error.
- **Meta-analysis.** The chapter first explains and runs the quality control (which estimates are usable, at least three studies per OTU, comparability of the studies), then explains the random-effects model (REM): fixed versus random effects, tau2, weights, I2 and the Q test. Two OTUs are worked through step by step, with a leave-one-out check; one function then pools every OTU.
- **Cladogram.** One cladogram for all included studies, built from the ANCOM-BC output and labelled with taxon names. Tree nodes show the pooled result; one outer ring per included study show in which condition each OTU is higher in each study (salmon = C1, teal = C2, olive = not significant).
- **Figure style.** Colours and theme follow Marzouk et al., npj Parkinson's Disease 2025.

Also updated to match: `workflow.Rmd`, `scripts/downstream_by_study.R`, `scripts/install_packages.R` (ANCOMBC instead of lefser), `README.md`, `results/mothur/README.md`, and the Word documents (README, beginner guide, prerequisites, mothur results README).

Notes:
- The mothur steps and their outputs are unchanged. The mothur command `dist.shared` still calculates only Bray-Curtis and Jclass; the four distances used in the figures are calculated in R, without subsampling.
- R results were produced with R 4.5.3, vegan 2.7-5, ANCOMBC 2.14.0, metafor 5.0-1 and MetaQC 0.1.14 from the tables in `results/mothur`.
- ANCOM-BC estimates for OTUs found in a single sample of each condition (p-value reported as 1) are not carried into the meta-analysis.
