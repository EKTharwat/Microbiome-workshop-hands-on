# Downstream analysis with the five studies (G1-G5) side by side in every figure.
# Run from the workshop folder:  source("scripts/downstream_by_study.R")
# Inputs : results/mothur/final.opti_mcc.0.03.shared, its cons.taxonomy, metadata/workshop_samples.tsv
# Outputs: report_outputs/by_study/  (figures as PNG, tables as CSV)
# The same commands are explained one by one in workflow-code-side-by-side.html.

# ---- All packages, loaded once ----------------------------------------------
packages <- c("vegan", "ggplot2", "patchwork", "ANCOMBC", "metafor", "MetaQC")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Not installed: ", paste(missing_packages, collapse = ", "), ". Run source(\"scripts/prepare_packages.R\") again.")
library(vegan)
library(ggplot2)
library(patchwork)
library(ANCOMBC)
library(metafor)
suppressPackageStartupMessages(library(MetaQC))

# ---- Import the mothur tables (same objects as the Import chapter) ----------
metadata <- read.delim("metadata/workshop_samples.tsv", check.names = FALSE)
counts_table <- read.delim("results/mothur/final.opti_mcc.0.03.shared", check.names = FALSE)
stopifnot(!anyDuplicated(metadata$sample_name), setequal(counts_table$Group, metadata$sample_name))
metadata <- metadata[match(counts_table$Group, metadata$sample_name), , drop = FALSE]
rownames(metadata) <- metadata$sample_name
metadata$condition <- factor(metadata$condition, levels = c("C1", "C2"))
count_matrix <- as.matrix(counts_table[, -(1:3), drop = FALSE])
storage.mode(count_matrix) <- "numeric"
rownames(count_matrix) <- counts_table$Group
taxonomy <- read.delim("results/mothur/final.opti_mcc.0.03.0.03.cons.taxonomy", check.names = FALSE)
taxonomy <- taxonomy[match(colnames(count_matrix), taxonomy$OTU), , drop = FALSE]
taxonomy_clean <- gsub("\"", "", gsub("\\([0-9.]+\\)", "", taxonomy$Taxonomy), fixed = FALSE)
tax_matrix <- as.matrix(read.table(text = taxonomy_clean, sep = ";", fill = TRUE, quote = "", comment.char = "", stringsAsFactors = FALSE))
rownames(tax_matrix) <- taxonomy$OTU
colnames(tax_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")[seq_len(ncol(tax_matrix))]
results_dir <- "report_outputs"

# ---- Settings, colours and one table per study ----------
# List the five studies and set the colours used in every figure.
studies <- c("G1", "G2", "G3", "G4", "G5")
condition_colours <- c(C1 = "#E8847A", C2 = "#4DB3B8")
not_significant_colour <- "#C9C462"
# Set one plot theme for all figures: white panels with a black frame and the legend underneath.
theme_set(theme_bw(base_size = 12) +
  theme(panel.border = element_rect(colour = "black", linewidth = 0.8),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", colour = "black"),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom"))
# Create the output folder for this chapter.
out_dir <- file.path(results_dir, "by_study")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# Make a readable name for every OTU: its number followed by its genus.
taxon_label <- setNames(paste(rownames(tax_matrix), tax_matrix[, "Genus"]), rownames(tax_matrix))
head(taxon_label)
# Split the counts and metadata into one piece per study, with the three C1 samples first and the three C2 samples last.
study_data <- lapply(setNames(studies, studies), function(s) {
  m <- metadata[metadata$study == s, , drop = FALSE]
  m <- m[order(m$condition), , drop = FALSE]
  stopifnot(nrow(m) == 6, all(table(m$condition) == 3))
  mat <- count_matrix[rownames(m), , drop = FALSE]
  list(m = m, mat = mat, relative = mat / rowSums(mat))
})
# Set the random seed so that every run gives the same numbers.
set.seed(20261114)

# ---- Quality control: reads per sample ----------
# Count the reads of every sample.
qc <- do.call(rbind, lapply(studies, function(s) {
  d <- study_data[[s]]
  data.frame(study = s, sample_name = rownames(d$m), condition = d$m$condition, reads = rowSums(d$mat))
}))
write.csv(qc, file.path(out_dir, "qc_by_study.csv"), row.names = FALSE)
# Draw the reads per sample, one panel per study.
qc_plot <- ggplot(qc, aes(sample_name, reads, fill = condition)) +
  geom_col(colour = "black", linewidth = 0.3) +
  facet_wrap(~study, nrow = 1, scales = "free") +
  scale_fill_manual(values = condition_colours) +
  labs(title = "Reads per sample after mothur", x = "Sample", y = "Reads", fill = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(out_dir, "qc_by_study.png"), qc_plot, width = 14, height = 4.2, dpi = 180)
print(qc_plot)

# ---- MetaQC: can the studies be combined? ----------
# Prepare the MetaQC input. For every study, log-ratio transform the counts (CLR) and keep the OTUs seen in at least two of its six samples.
metaqc_input <- lapply(study_data, function(d) {
  present <- d$mat[, colSums(d$mat > 0) > 0, drop = FALSE]
  clr <- log(present + 1)
  clr <- clr - rowMeans(clr)
  keep <- colSums(present > 0) >= 2 & apply(clr, 2, var) > 0
  list(x = t(clr[, keep, drop = FALSE]), y = as.character(d$m$condition))
})
sapply(metaqc_input, function(d) nrow(d$x))
# Write the taxonomic groups in the GMT format that MetaQC reads for its pathway-based measures (EQC, CQCp, AQCp).
metaqc_dir <- file.path(out_dir, "metaqc")
dir.create(metaqc_dir, recursive = TRUE, showWarnings = FALSE)
taxon_sets <- unlist(lapply(c("Phylum", "Class", "Order", "Family", "Genus"), function(rank)
  split(rownames(tax_matrix), paste(rank, tax_matrix[, rank], sep = "_"))), recursive = FALSE)
taxon_sets <- taxon_sets[lengths(taxon_sets) >= 2 & !duplicated(taxon_sets)]
writeLines(vapply(names(taxon_sets), function(n) paste(c(n, "taxonomy", taxon_sets[[n]]), collapse = "\t"), character(1)),
           file.path(metaqc_dir, "taxon_sets.gmt"))
lengths(taxon_sets)
# Adapt the IQC calculation to small tables.
iqc_all_shared <- function(.) {
  pairs <- combn(length(.$.DListF), 2)
  distance <- apply(pairs, 2, function(ii) {
    shared <- intersect(rownames(.$.DListF[[ii[1]]]$x), rownames(.$.DListF[[ii[2]]]$x))
    correlations <- sapply(1:2, function(j) as.vector(as.dist(cor(t(.$.DListF[[ii[j]]]$x[shared, , drop = FALSE])))))
    correlations <- na.omit(correlations)
    (1 - cor(correlations[, 1], correlations[, 2], method = "spearman")) / 2
  })
  m <- matrix(0, length(.$.DListF), length(.$.DListF), dimnames = list(names(.$.DListF), names(.$.DListF)))
  m[t(pairs)] <- distance
  m[t(pairs[2:1, ])] <- distance
  .$.DistOfStudies <- as.dist(m)
}
# Run MetaQC: create the MetaQC object, plug in the IQC function and compute all six measures with 10,000 permutations.
foreach::registerDoSEQ()
run_metaqc <- function() {
  old_dir <- setwd(metaqc_dir)
  on.exit(setwd(old_dir))
  set.seed(20261114)
  object <- MetaQC(metaqc_input, "taxon_sets.gmt", filterGenes = FALSE, minNumGenes = 2, verbose = FALSE)
  object$.CalcDistOfStudies <- iqc_all_shared
  suppressWarnings(runQC(object, B = 1e4, fileForCQCp = "taxon_sets.gmt"))
  object
}
metaqc <- run_metaqc()
metaqc$Print()
# Print how different the studies are from each other (the IQC distances).
round(metaqc$.DistOfStudies, 2)
# Draw the six MetaQC scores of every study with the pass line.
metaqc_threshold <- -log10(0.05 / length(studies))
metaqc_scores <- metaqc$.Scores
metaqc_long <- data.frame(study = rep(rownames(metaqc_scores), ncol(metaqc_scores)),
                          measure = factor(rep(colnames(metaqc_scores), each = nrow(metaqc_scores)), levels = colnames(metaqc_scores)),
                          score = unlist(metaqc_scores))
metaqc_long$result <- ifelse(metaqc_long$score >= metaqc_threshold, "Passes", "Does not pass")
write.csv(cbind(study = rownames(metaqc_scores), round(metaqc_scores, 3)), file.path(out_dir, "metaqc_scores.csv"), row.names = FALSE)
metaqc_plot <- ggplot(metaqc_long, aes(study, score, fill = result)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.7) +
  geom_hline(yintercept = metaqc_threshold, linetype = "dashed") +
  geom_text(aes(label = sprintf("%.2f", pmax(score, 0))), vjust = -0.4, size = 3) +
  facet_wrap(~measure, nrow = 1) +
  scale_fill_manual(values = c("Passes" = condition_colours[["C2"]], "Does not pass" = not_significant_colour)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "MetaQC quality scores per study", subtitle = sprintf("Score = -log10(p); dashed line = -log10(0.05/%d) = %.1f", length(studies), metaqc_threshold),
       x = "Study", y = "Score", fill = NULL)
ggsave(file.path(out_dir, "metaqc_scores.png"), metaqc_plot, width = 14, height = 4.6, dpi = 180)
print(metaqc_plot)
# Apply the decision rule: exclude a study that fails both IQC and CQCg, and keep the others for every later step.
metaqc_decision <- data.frame(study = rownames(metaqc_scores), IQC = round(metaqc_scores$IQC, 2), CQCg = round(metaqc_scores$CQCg, 2),
                              passes_IQC = metaqc_scores$IQC >= metaqc_threshold, passes_CQCg = metaqc_scores$CQCg >= metaqc_threshold)
metaqc_decision$included <- metaqc_decision$passes_IQC | metaqc_decision$passes_CQCg
metaqc_decision <- metaqc_decision[order(metaqc_decision$study), ]
write.csv(metaqc_decision, file.path(out_dir, "metaqc_decision.csv"), row.names = FALSE)
print(metaqc_decision, row.names = FALSE)
excluded_studies <- metaqc_decision$study[!metaqc_decision$included]
studies <- setdiff(studies, excluded_studies)
study_data <- study_data[studies]
studies

# ---- Alpha diversity: Observed, Chao1, Shannon and Inverse Simpson ----------
# Calculate the four alpha-diversity indices for every sample.
alpha <- do.call(rbind, lapply(studies, function(s) {
  d <- study_data[[s]]
  data.frame(study = s, sample_name = rownames(d$mat), condition = d$m$condition,
             observed = rowSums(d$mat > 0),
             chao1 = estimateR(d$mat)["S.chao1", ],
             shannon = diversity(d$mat, "shannon"),
             inverse_simpson = diversity(d$mat, "invsimpson"))
}))
write.csv(alpha, file.path(out_dir, "alpha_by_study.csv"), row.names = FALSE)
print(alpha, row.names = FALSE, digits = 3)
# Write the exact permutation test: compare the observed C2-minus-C1 difference with all 20 ways of splitting six samples into two groups of three.
allocations <- combn(6, 3)
exact_test <- function(x) {
  observed <- mean(x[4:6]) - mean(x[1:3])
  permuted <- apply(allocations, 2, function(i) mean(x[-i]) - mean(x[i]))
  c(difference = observed, p = mean(abs(permuted) >= abs(observed) - 1e-12))
}
# Test every index in every study and adjust the p-values within each study.
alpha_metrics <- c(observed = "Observed OTUs", chao1 = "Chao1 richness", shannon = "Shannon", inverse_simpson = "Inverse Simpson")
alpha_tests <- do.call(rbind, lapply(studies, function(s) do.call(rbind, lapply(names(alpha_metrics), function(metric) {
  r <- exact_test(alpha[alpha$study == s, metric])
  data.frame(study = s, metric = metric, difference_C2_minus_C1 = r[["difference"]], p_value = r[["p"]])
}))))
alpha_tests$adjusted_p <- ave(alpha_tests$p_value, alpha_tests$study, FUN = function(p) p.adjust(p, "BH"))
write.csv(alpha_tests, file.path(out_dir, "alpha_tests_by_study.csv"), row.names = FALSE)
print(alpha_tests, row.names = FALSE, digits = 3)
# Draw the four indices as violin plots: studies along the x-axis, the two conditions next to each other, and the exact p-value above each study.
alpha_long <- do.call(rbind, lapply(names(alpha_metrics), function(metric)
  data.frame(alpha[, c("study", "sample_name", "condition")], metric = alpha_metrics[[metric]], value = alpha[[metric]])))
alpha_long$metric <- factor(alpha_long$metric, levels = alpha_metrics)
alpha_labels <- do.call(rbind, lapply(names(alpha_metrics), function(metric) {
  top <- aggregate(value ~ study, alpha_long[alpha_long$metric == alpha_metrics[[metric]], ], max)
  tests <- alpha_tests[alpha_tests$metric == metric, ]
  data.frame(study = top$study, metric = alpha_metrics[[metric]], value = top$value,
             label = sprintf("p = %.1f", tests$p_value[match(top$study, tests$study)]))
}))
alpha_labels$metric <- factor(alpha_labels$metric, levels = alpha_metrics)
alpha_plot <- ggplot(alpha_long, aes(study, value, fill = condition)) +
  geom_violin(trim = TRUE, scale = "width", width = 0.7, colour = "black", linewidth = 0.4, position = position_dodge(0.8)) +
  geom_boxplot(width = 0.12, colour = "black", linewidth = 0.3, outlier.shape = NA, position = position_dodge(0.8), show.legend = FALSE) +
  geom_point(shape = 21, size = 1.6, colour = "black", position = position_dodge(0.8), show.legend = FALSE) +
  geom_text(data = alpha_labels, aes(study, Inf, label = label), inherit.aes = FALSE, vjust = 1.6, size = 3.2) +
  facet_wrap(~metric, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = condition_colours) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(title = "Alpha diversity by study and condition", x = "Study", y = "Index value", fill = "Condition")
ggsave(file.path(out_dir, "alpha_by_study.png"), alpha_plot, width = 12, height = 8, dpi = 180)
print(alpha_plot)

# ---- Beta diversity: Bray-Curtis, Canberra, Jclass and ThetaYC ----------
# Write the ThetaYC distance, which vegan does not include.
thetayc <- function(relative) {
  cross <- relative %*% t(relative)
  squares <- diag(cross)
  as.dist(1 - cross / (outer(squares, squares, "+") - cross))
}
# List the four distances. Each entry is a small function that takes one study and returns its distance matrix.
distance_functions <- list(
  "Bray-Curtis" = function(d) vegdist(d$relative, method = "bray"),
  "Canberra"    = function(d) vegdist(d$relative, method = "canberra"),
  "Jclass"      = function(d) vegdist(d$mat, method = "jaccard", binary = TRUE),
  "ThetaYC"     = function(d) thetayc(d$relative))
# For every distance and every study: run the PCoA, the PERMANOVA and the PERMDISP test.
beta <- lapply(names(distance_functions), function(distance_name) lapply(studies, function(s) {
  d <- study_data[[s]]
  distance <- distance_functions[[distance_name]](d)
  pcoa <- cmdscale(distance, k = 2, eig = TRUE)
  explained <- 100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0])
  permanova <- suppressMessages(adonis2(distance ~ condition, data = d$m, permutations = 999))
  permdisp <- suppressMessages(vegan::permutest(suppressWarnings(betadisper(distance, d$m$condition)), permutations = 999))
  list(points = data.frame(distance = distance_name, study = s, sample_name = rownames(d$m), condition = d$m$condition,
                           PCoA1 = pcoa$points[, 1], PCoA2 = pcoa$points[, 2]),
       tests = data.frame(distance = distance_name, study = s, axis1_percent = explained[1], axis2_percent = explained[2],
                          R2 = permanova$R2[1], F = permanova$F[1], permanova_p = permanova$`Pr(>F)`[1],
                          permdisp_p = permdisp$tab$`Pr(>F)`[1]))
}))
pcoa_points <- do.call(rbind, lapply(beta, function(x) do.call(rbind, lapply(x, `[[`, "points"))))
beta_tests <- do.call(rbind, lapply(beta, function(x) do.call(rbind, lapply(x, `[[`, "tests"))))
write.csv(pcoa_points, file.path(out_dir, "pcoa_by_study.csv"), row.names = FALSE)
write.csv(beta_tests, file.path(out_dir, "permanova_by_study.csv"), row.names = FALSE)
# Print the PERMANOVA and PERMDISP results with the percentage of variation on the two PCoA axes.
print(beta_tests, row.names = FALSE, digits = 3)
# Write a function that returns the outline of an ellipse around the samples of one condition.
ellipse_points <- function(points, level = 0.68) {
  xy <- as.matrix(points[, c("PCoA1", "PCoA2")])
  spread <- eigen(cov(xy), symmetric = TRUE)
  circle <- rbind(cos(seq(0, 2 * pi, length.out = 100)), sin(seq(0, 2 * pi, length.out = 100)))
  outline <- t(colMeans(xy) + spread$vectors %*% (sqrt(pmax(spread$values, 0) * qchisq(level, 2)) * circle))
  data.frame(condition = points$condition[1], PCoA1 = outline[, 1], PCoA2 = outline[, 2])
}
# Write a function that draws the PCoA of one distance, one panel per study, with the percentage of variation written on the x- and y-axis.
pcoa_figure <- function(distance_name) {
  panels <- lapply(studies, function(s) {
    points <- pcoa_points[pcoa_points$distance == distance_name & pcoa_points$study == s, ]
    tests <- beta_tests[beta_tests$distance == distance_name & beta_tests$study == s, ]
    ellipses <- do.call(rbind, lapply(split(points, points$condition), ellipse_points))
    ggplot(points, aes(PCoA1, PCoA2, colour = condition, fill = condition)) +
      geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
      geom_polygon(data = ellipses, alpha = 0.2, linewidth = 0.5) +
      geom_point(size = 3) +
      geom_text(aes(label = sample_name), vjust = -1, size = 3, show.legend = FALSE) +
      scale_fill_manual(values = condition_colours) +
      scale_colour_manual(values = condition_colours) +
      scale_x_continuous(expand = expansion(mult = 0.2)) + scale_y_continuous(expand = expansion(mult = 0.2)) +
      labs(title = s,
           subtitle = sprintf("PERMANOVA R2 = %.2f, p = %.2f\nPERMDISP p = %.2f", tests$R2, tests$permanova_p, tests$permdisp_p),
           x = sprintf("PCoA1 (%.1f%%)", tests$axis1_percent), y = sprintf("PCoA2 (%.1f%%)", tests$axis2_percent),
           colour = "Condition", fill = "Condition") +
      theme(plot.subtitle = element_text(size = 8.5), panel.grid = element_blank())
  })
  wrap_plots(panels, nrow = 1, guides = "collect") +
    plot_annotation(title = paste(distance_name, "PCoA")) & theme(legend.position = "bottom")
}
# Draw and save the Bray-Curtis PCoA.
pcoa_bray <- pcoa_figure("Bray-Curtis")
ggsave(file.path(out_dir, "pcoa_braycurtis_by_study.png"), pcoa_bray, width = 17, height = 4.8, dpi = 180)
print(pcoa_bray)
# Draw and save the Canberra PCoA.
pcoa_canberra <- pcoa_figure("Canberra")
ggsave(file.path(out_dir, "pcoa_canberra_by_study.png"), pcoa_canberra, width = 17, height = 4.8, dpi = 180)
print(pcoa_canberra)
# Draw and save the Jclass PCoA.
pcoa_jclass <- pcoa_figure("Jclass")
ggsave(file.path(out_dir, "pcoa_jclass_by_study.png"), pcoa_jclass, width = 17, height = 4.8, dpi = 180)
print(pcoa_jclass)
# Draw and save the ThetaYC PCoA.
pcoa_thetayc <- pcoa_figure("ThetaYC")
ggsave(file.path(out_dir, "pcoa_thetayc_by_study.png"), pcoa_thetayc, width = 17, height = 4.8, dpi = 180)
print(pcoa_thetayc)
# Summarise PERMANOVA in one picture: the share of variation explained by condition (R2) for every distance and study, with the p-value underneath.
beta_tests$distance <- factor(beta_tests$distance, levels = rev(names(distance_functions)))
permanova_plot <- ggplot(beta_tests, aes(study, distance, fill = R2)) +
  geom_tile(colour = "black") +
  geom_text(aes(label = sprintf("R2 = %.2f\np = %.2f", R2, permanova_p)), size = 3.4) +
  scale_fill_gradient(low = "white", high = condition_colours[["C2"]], limits = c(0, 1)) +
  labs(title = "PERMANOVA: variation explained by condition", x = "Study", y = NULL, fill = "R2") +
  theme(panel.grid = element_blank(), legend.position = "right")
ggsave(file.path(out_dir, "permanova_by_study.png"), permanova_plot, width = 9, height = 3.8, dpi = 180)
print(permanova_plot)

# ---- ANCOM-BC per study: effect size and standard error for every OTU ----------
# Run ANCOM-BC in every study and collect, for every OTU, the effect size (log fold change, C2 versus C1), its standard error, the p-value and the adjusted p-value.
ancom <- do.call(rbind, lapply(studies, function(s) {
  d <- study_data[[s]]
  keep <- colSums(d$mat > 0) >= 2
  set.seed(20261114)
  fit <- suppressWarnings(suppressMessages(ancombc2(data = t(d$mat[, keep, drop = FALSE]), taxa_are_rows = TRUE,
           meta_data = d$m, fix_formula = "condition", group = "condition", p_adj_method = "BH",
           prv_cut = 0, lib_cut = 0, struc_zero = FALSE, pseudo_sens = FALSE, verbose = FALSE)))
  r <- fit$res
  data.frame(study = s, feature = r$taxon, lfc_C2_vs_C1 = r$lfc_conditionC2, se = r$se_conditionC2,
             W = r$W_conditionC2, p_value = r$p_conditionC2, q_value = r$q_conditionC2)
}))
ancom$taxon <- unname(taxon_label[ancom$feature])
ancom$usable <- is.finite(ancom$lfc_C2_vs_C1) & is.finite(ancom$se) & ancom$se > 0 & ancom$p_value < 1
# Label each estimate: higher in C1, higher in C2, or not significant.
significance_level <- 0.05
ancom$direction <- ifelse(!ancom$usable | ancom$p_value >= significance_level, "Not significant",
                          ifelse(ancom$lfc_C2_vs_C1 > 0, "Higher in C2", "Higher in C1"))
ancom$direction <- factor(ancom$direction, levels = c("Higher in C1", "Higher in C2", "Not significant"))
direction_colours <- c("Higher in C1" = condition_colours[["C1"]], "Higher in C2" = condition_colours[["C2"]], "Not significant" = not_significant_colour)
write.csv(ancom, file.path(out_dir, "ancombc_by_study.csv"), row.names = FALSE)
# Print the effect size of every OTU in every study.
print(ancom[, c("study", "taxon", "lfc_C2_vs_C1", "se", "p_value", "q_value", "usable", "direction")], row.names = FALSE, digits = 3)
# Count the estimates and how many are significant.
table(study = ancom$study, result = ancom$direction)
c(estimates = nrow(ancom), usable = sum(ancom$usable), p_below_0.05 = sum(ancom$usable & ancom$p_value < 0.05), q_below_0.05 = sum(ancom$usable & ancom$q_value < 0.05))
# Draw the effect size of every OTU with its 95% confidence interval, one panel per study.
ancom_drawn <- ancom[ancom$usable, ]
ancom_plot <- ggplot(ancom_drawn, aes(lfc_C2_vs_C1, reorder(paste(study, taxon, sep = "|"), lfc_C2_vs_C1), fill = direction)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_errorbar(aes(xmin = lfc_C2_vs_C1 - 1.96 * se, xmax = lfc_C2_vs_C1 + 1.96 * se), width = 0, orientation = "y") +
  geom_point(shape = 21, size = 3, colour = "black") +
  facet_wrap(~study, nrow = 1, scales = "free_y") +
  scale_y_discrete(labels = function(x) sub("^[^|]*\\|", "", x)) +
  scale_fill_manual(values = direction_colours, drop = FALSE) +
  labs(title = "ANCOM-BC effect sizes: C2 versus C1", x = "Log fold change (natural log), with 95% confidence interval", y = NULL, fill = NULL) +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(out_dir, "ancombc_by_study.png"), ancom_plot, width = 18, height = 5.6, dpi = 180)
print(ancom_plot)

# ---- Random-effects meta-analysis (metafor, REML) for every OTU ----------
# Before pooling: for every OTU, list which studies give a usable estimate and whether the OTU can be pooled.
qc_meta <- expand.grid(feature = unique(ancom$feature), study = studies, stringsAsFactors = FALSE)
found <- match(paste(qc_meta$feature, qc_meta$study), paste(ancom$feature, ancom$study))
qc_meta$status <- ifelse(is.na(found), "Not detected", ifelse(ancom$usable[found], "Usable", "Not testable"))
qc_meta$taxon <- unname(taxon_label[qc_meta$feature])
qc_summary <- as.data.frame.matrix(table(qc_meta$taxon, factor(qc_meta$status, levels = c("Usable", "Not testable", "Not detected"))))
qc_summary$pooled <- ifelse(qc_summary$Usable >= 3, "yes", "no")
write.csv(qc_meta, file.path(out_dir, "meta_qc.csv"), row.names = FALSE)
print(qc_summary)
# Before pooling: draw the same information as a map of OTUs by studies.
qc_meta$status <- factor(qc_meta$status, levels = c("Usable", "Not testable", "Not detected"))
usable_count <- tapply(qc_meta$status == "Usable", qc_meta$taxon, sum)
qc_meta$taxon <- factor(qc_meta$taxon, levels = names(sort(usable_count)))
qc_meta_plot <- ggplot(qc_meta, aes(study, taxon, fill = status)) +
  geom_tile(colour = "black") +
  geom_hline(yintercept = sum(usable_count < 3) + 0.5, linewidth = 1.2) +
  scale_fill_manual(values = c("Usable" = condition_colours[["C2"]], "Not testable" = not_significant_colour, "Not detected" = "white")) +
  labs(title = "Estimates available for pooling", subtitle = "OTUs above the thick line have a usable estimate in at least three studies and are pooled",
       x = "Study", y = NULL, fill = "ANCOM-BC estimate") +
  theme(panel.grid = element_blank(), plot.subtitle = element_text(size = 9))
ggsave(file.path(out_dir, "meta_qc.png"), qc_meta_plot, width = 9, height = 6.5, dpi = 180)
print(qc_meta_plot)
# Example 1, Otu05: take its ANCOM-BC effect sizes from the included studies.
example_1 <- ancom[ancom$feature == "Otu05" & ancom$usable, ]
print(example_1[, c("study", "taxon", "lfc_C2_vs_C1", "se", "p_value")], row.names = FALSE, digits = 3)
# Example 1, Otu05: fit the random-effects model and print the pooled effect size.
model_1 <- rma(yi = lfc_C2_vs_C1, sei = se, data = example_1, method = "REML")
print(model_1)
# Example 1, Otu05: draw the forest plot and save it. Each square is one study; the diamond is the pooled effect size.
forest_1 <- function() forest(model_1, slab = example_1$study, header = c("Study", "Log fold change [95% CI]"),
  xlab = "Log fold change, C2 versus C1", main = taxon_label[["Otu05"]], colout = "black",
  col = condition_colours[[ifelse(coef(model_1) > 0, "C2", "C1")]], border = "black")
forest_1()
png(file.path(out_dir, "meta_example_Otu05.png"), width = 1500, height = 850, res = 180)
forest_1()
invisible(dev.off())
# Example 1, Otu05: quality control of the pooled result. Refit the model once per study, each time leaving that study out.
leave_one_out <- leave1out(model_1)
print(data.frame(left_out = example_1$study, estimate = leave_one_out$estimate, ci_lower = leave_one_out$ci.lb, ci_upper = leave_one_out$ci.ub,
                 p_value = leave_one_out$pval, I_squared = leave_one_out$I2), row.names = FALSE, digits = 3)
# Example 2, Otu16: take its ANCOM-BC effect sizes from the included studies.
example_2 <- ancom[ancom$feature == "Otu16" & ancom$usable, ]
print(example_2[, c("study", "taxon", "lfc_C2_vs_C1", "se", "p_value")], row.names = FALSE, digits = 3)
# Example 2, Otu16: fit the random-effects model and print the pooled effect size.
model_2 <- rma(yi = lfc_C2_vs_C1, sei = se, data = example_2, method = "REML")
print(model_2)
# Example 2, Otu16: draw the forest plot and save it.
forest_2 <- function() forest(model_2, slab = example_2$study, header = c("Study", "Log fold change [95% CI]"),
  xlab = "Log fold change, C2 versus C1", main = taxon_label[["Otu16"]], colout = "black",
  col = condition_colours[[ifelse(coef(model_2) > 0, "C2", "C1")]], border = "black")
forest_2()
png(file.path(out_dir, "meta_example_Otu16.png"), width = 1500, height = 850, res = 180)
forest_2()
invisible(dev.off())
# All OTUs: write one function that does for any OTU what the two examples did by hand.
pool_otu <- function(feature) {
  d <- ancom[ancom$feature == feature & ancom$usable, ]
  not_estimable <- data.frame(feature = feature, studies = nrow(d), estimate = NA, se = NA, ci_lower = NA, ci_upper = NA,
                              p_value = NA, tau_squared = NA, I_squared = NA, Q_p_value = NA)
  if (nrow(d) < 3) return(not_estimable)
  model <- tryCatch(rma(yi = lfc_C2_vs_C1, sei = se, data = d, method = "REML"), error = function(e) NULL)
  if (is.null(model)) return(not_estimable)
  data.frame(feature = feature, studies = model$k, estimate = as.numeric(coef(model)), se = model$se, ci_lower = model$ci.lb,
             ci_upper = model$ci.ub, p_value = model$pval, tau_squared = model$tau2, I_squared = model$I2, Q_p_value = model$QEp)
}
# All OTUs: apply the function to every OTU, adjust the p-values for multiple testing and save the table.
meta_otu <- do.call(rbind, lapply(unique(ancom$feature), pool_otu))
meta_otu$adjusted_p <- p.adjust(meta_otu$p_value, "BH")
meta_otu$taxon <- unname(taxon_label[meta_otu$feature])
meta_otu$direction <- ifelse(is.na(meta_otu$p_value) | meta_otu$p_value >= significance_level, "Not significant",
                             ifelse(meta_otu$estimate > 0, "Higher in C2", "Higher in C1"))
meta_otu$direction <- factor(meta_otu$direction, levels = names(direction_colours))
meta_otu <- meta_otu[order(meta_otu$p_value), ]
write.csv(meta_otu, file.path(out_dir, "meta_otu_ancombc.csv"), row.names = FALSE)
print(meta_otu[, c("taxon", "studies", "estimate", "ci_lower", "ci_upper", "p_value", "adjusted_p", "tau_squared", "I_squared", "Q_p_value")], row.names = FALSE, digits = 3)
# All OTUs: draw the pooled effect size of every OTU with its 95% confidence interval.
pooled <- meta_otu[!is.na(meta_otu$estimate), ]
pooled$label <- sprintf("%s  (k = %d)%s", pooled$taxon, pooled$studies, ifelse(pooled$adjusted_p < 0.05, " *", ""))
meta_plot <- ggplot(pooled, aes(estimate, reorder(label, estimate), fill = direction)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0, orientation = "y") +
  geom_point(shape = 23, size = 4, colour = "black") +
  scale_fill_manual(values = direction_colours, drop = FALSE) +
  labs(title = "Pooled effect size per OTU", subtitle = "Random-effects meta-analysis of ANCOM-BC log fold changes; k = number of studies; * adjusted p < 0.05",
       x = "Pooled log fold change, C2 versus C1, with 95% confidence interval", y = NULL, fill = NULL) +
  theme(plot.subtitle = element_text(size = 9))
ggsave(file.path(out_dir, "meta_otu_ancombc.png"), meta_plot, width = 10, height = 0.34 * nrow(pooled) + 2.4, dpi = 180)
print(meta_plot)

# ---- One cladogram for all included studies, coloured from the ANCOM-BC results ----------
# Turn the taxonomy of the OTUs into tree coordinates: one row per node with its rank, its position around the circle and its parent.
features <- unique(ancom$feature)
ranks <- cbind(tax_matrix[features, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), drop = FALSE], OTU = unname(taxon_label[features]))
ranks <- ranks[do.call(order, as.data.frame(ranks)), , drop = FALSE]
paths <- sapply(seq_len(ncol(ranks)), function(level) apply(ranks[, seq_len(level), drop = FALSE], 1, paste, collapse = ";"))
leaf_position <- seq_len(nrow(ranks))
nodes <- do.call(rbind, lapply(seq_len(ncol(ranks)), function(level) {
  first <- match(unique(paths[, level]), paths[, level])
  data.frame(path = paths[first, level], level = level - 1, name = ranks[first, level],
             position = as.numeric(tapply(leaf_position, factor(paths[, level], levels = unique(paths[, level])), function(p) mean(range(p)))),
             parent = if (level == 1) NA else paths[first, level - 1])
}))
nodes$feature <- ifelse(nodes$level == 6, rownames(ranks)[match(nodes$path, paths[, 7])], NA)
nodes$parent_position <- nodes$position[match(nodes$parent, nodes$path)]
# Colour the tree nodes from the pooled ANCOM-BC effect sizes. A higher-rank node takes a colour when all the significant OTUs below it point the same way.
nodes$direction <- as.character(meta_otu$direction[match(nodes$feature, meta_otu$feature)])
for (i in which(is.na(nodes$feature))) {
  below <- rownames(ranks)[paths[, nodes$level[i] + 1] == nodes$path[i]]
  hits <- setdiff(unique(as.character(meta_otu$direction[match(below, meta_otu$feature)])), "Not significant")
  nodes$direction[i] <- if (length(hits) == 1) hits else "Not significant"
}
nodes$direction <- factor(nodes$direction, levels = names(direction_colours))
# Make the ring tiles: one tile per OTU and study, coloured by the condition in which the OTU is higher.
tiles <- ancom[ancom$usable, ]
tiles$position <- match(tiles$feature, rownames(ranks))
tiles$ring <- 8.9 + 0.55 * match(tiles$study, studies)
tiles$higher_in <- ifelse(tiles$lfc_C2_vs_C1 > 0, "C2", "C1")
tiles$significant <- tiles$p_value < significance_level
# Draw the cladogram and save it.
n <- nrow(ranks)
tips <- nodes[!is.na(nodes$feature), ]
angle <- 90 - 360 * (tips$position - 0.5) / n
tips$hjust <- ifelse(angle < -90, 1, 0)
tips$angle <- ifelse(angle < -90, angle + 180, angle)
inner <- nodes[nodes$level %in% c(1, 4), ]
sectors <- tips[tips$direction != "Not significant", ]
branches <- nodes[!is.na(nodes$parent), ]
ring_labels <- data.frame(position = 0, ring = 8.9 + 0.55 * seq_along(studies), label = studies)
cladogram_plot <- ggplot(nodes) +
  geom_rect(data = sectors, aes(xmin = position - 0.5, xmax = position + 0.5, ymin = 2, ymax = 8.4, fill = direction), alpha = 0.25, show.legend = FALSE) +
  geom_segment(data = branches, aes(x = position, xend = position, y = level + 1, yend = level + 2), linewidth = 0.5) +
  geom_segment(data = branches, aes(x = parent_position, xend = position, y = level + 1, yend = level + 1), linewidth = 0.5) +
  geom_point(aes(position, level + 2, fill = direction, size = is.na(feature)), shape = 21, colour = "black", show.legend = c(fill = TRUE, size = FALSE)) +
  geom_tile(data = tiles, aes(position, ring, fill = paste("Higher in", higher_in), alpha = significant), width = 0.9, height = 0.5, colour = "black", linewidth = 0.2, show.legend = FALSE) +
  geom_text(data = ring_labels, aes(position, ring, label = label), size = 2.4, fontface = "bold") +
  geom_label(data = inner, aes(position, level + 2, label = name, fontface = ifelse(level == 1, "bold", "plain")), size = 3,
             vjust = -0.6, fill = "white", alpha = 0.7, linewidth = 0, label.padding = unit(0.08, "lines")) +
  geom_text(data = tips, aes(position, 12.1, label = name, angle = angle, hjust = hjust, fontface = ifelse(direction == "Not significant", "plain", "bold")), size = 3.6) +
  scale_fill_manual(values = direction_colours, drop = FALSE, name = NULL,
                    labels = c("Higher in C1", "Higher in C2", "Not significant")) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3)) +
  scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 5)) +
  scale_x_continuous(limits = c(-0.5, n + 0.5)) + scale_y_continuous(limits = c(0, 18)) +
  coord_polar(theta = "x", clip = "off") +
  labs(title = "Cladogram of all included studies, coloured by condition",
       subtitle = paste0("Tree nodes: pooled ANCOM-BC result. Rings ", paste(studies, collapse = ", "), " (inside to outside): condition in which the OTU is higher in each study (strong colour = p < 0.05)")) +
  theme_void(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 10), legend.position = "bottom",
        plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(10, 30, 10, 30)) +
  guides(fill = guide_legend(override.aes = list(size = 5)))
ggsave(file.path(out_dir, "ancombc_cladogram.png"), cladogram_plot, width = 12, height = 12.5, dpi = 180)
print(cladogram_plot)

# ---- Effect size of alpha diversity across the included studies ----------
# Calculate the C2-minus-C1 difference and its variance for every index and study.
alpha_effects <- do.call(rbind, lapply(names(alpha_metrics), function(metric) do.call(rbind, lapply(studies, function(s) {
  x <- alpha[alpha$study == s, metric]
  data.frame(metric = metric, study = s, effect = mean(x[4:6]) - mean(x[1:3]), variance = var(x[1:3]) / 3 + var(x[4:6]) / 3)
}))))
print(alpha_effects, row.names = FALSE, digits = 3)
# Pool the included studies for every index with a random-effects model.
alpha_models <- lapply(setNames(names(alpha_metrics), names(alpha_metrics)), function(metric)
  rma(yi = effect, vi = variance, data = alpha_effects[alpha_effects$metric == metric & alpha_effects$variance > 0, ], method = "REML"))
meta_alpha <- do.call(rbind, lapply(names(alpha_models), function(metric) { model <- alpha_models[[metric]]
  data.frame(metric = metric, studies = model$k, estimate = as.numeric(coef(model)), ci_lower = model$ci.lb, ci_upper = model$ci.ub,
             p_value = model$pval, tau_squared = model$tau2, I_squared = model$I2) }))
write.csv(meta_alpha, file.path(out_dir, "meta_alpha.csv"), row.names = FALSE)
print(meta_alpha, row.names = FALSE, digits = 3)
# Draw the forest plots of the four indices and save them.
forest_alpha <- function() {
  old_par <- par(mfrow = c(2, 2), mar = c(4.5, 2, 3, 2))
  for (metric in names(alpha_models)) {
    model <- alpha_models[[metric]]
    forest(model, slab = alpha_effects$study[alpha_effects$metric == metric & alpha_effects$variance > 0], header = "Study",
           xlab = "Difference, C2 minus C1", main = alpha_metrics[[metric]], colout = "black",
           col = condition_colours[[ifelse(coef(model) > 0, "C2", "C1")]], border = "black")
  }
  par(old_par)
}
forest_alpha()
png(file.path(out_dir, "meta_alpha_forest.png"), width = 2400, height = 1500, res = 180)
forest_alpha()
invisible(dev.off())

message("Finished. Figures and tables are in ", out_dir)
