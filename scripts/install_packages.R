# This helper runs in a separate R process, so packages can be updated safely.
cran_packages <- c("BiocManager", "vegan", "ggplot2", "patchwork", "knitr", "rmarkdown", "digest", "jsonlite", "metafor", "ape", "processx", "proto", "foreach", "iterators")
installed <- rownames(installed.packages())
missing_cran <- setdiff(cran_packages, installed)
minimum_versions <- c(ggplot2 = "4.0.0", patchwork = "1.3.0")
for (package in names(minimum_versions)) {
  if (package %in% installed && utils::packageVersion(package) < minimum_versions[[package]])
    missing_cran <- union(missing_cran, package)
}
if (length(missing_cran)) {
  message("Installing required R packages: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}
installed <- rownames(installed.packages())
missing_cran <- setdiff(cran_packages, installed)
if (length(missing_cran)) stop("R package installation did not complete: ", paste(missing_cran, collapse = ", "))
bioc_packages <- c("phyloseq", "ANCOMBC")
missing_bioc <- setdiff(bioc_packages, installed)
if (length(missing_bioc)) {
  options(BioC_mirror = "https://bioconductor.posit.co")
  # The mirror also hosts the version map; do not require an existing BiocVersion package.
  Sys.setenv(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = "TRUE",
             BIOCONDUCTOR_CONFIG_FILE = "https://bioconductor.posit.co/config.yaml")
  tryCatch(BiocManager::install(missing_bioc, ask = FALSE, update = FALSE),
    error = function(e) stop("Bioconductor setup failed: ", conditionMessage(e),
      "\nCheck your internet connection and retry source(\"scripts/prepare_packages.R\")."))
}
# MetaQC is no longer distributed by Bioconductor. A copy of its source from github.com/donkang75/MetaQC
# (version 0.1.14, GPL-2) is bundled in vendor/, so it installs without internet access. It is pure R: no compiler needed.
if (!"MetaQC" %in% rownames(installed.packages())) {
  message("Installing MetaQC from vendor/MetaQC_0.1.14.tar.gz")
  install.packages("vendor/MetaQC_0.1.14.tar.gz", repos = NULL, type = "source")
}
required <- c(cran_packages, bioc_packages, "MetaQC")
missing <- setdiff(required, rownames(installed.packages()))
if (length(missing)) stop("Package installation did not complete: ", paste(missing, collapse = ", "))
for (package in names(minimum_versions)) {
  if (utils::packageVersion(package) < minimum_versions[[package]])
    stop(package, " could not be updated. Restart RStudio and run the starter again.")
}
# Try to load every package. If one fails because a dependency is missing, install that dependency and try again.
load_error <- function(package) tryCatch({ loadNamespace(package); "" }, error = function(e) conditionMessage(e))
for (attempt in 1:6) {
  errors <- vapply(required, load_error, character(1))
  errors <- errors[nzchar(errors)]
  if (!length(errors)) break
  found <- regmatches(errors, regexec("there is no package called [^[:alnum:]]*([[:alnum:].]+)", errors))
  missing_dependencies <- unique(unlist(lapply(found, function(x) if (length(x) == 2) x[2])))
  if (!length(missing_dependencies)) {
    # Not a missing dependency: the installed copy itself is broken or too old (for example an old ANCOMBC that needs
    # a function the current CVXR no longer has). Reinstall the failing Bioconductor packages once, then check again.
    broken <- intersect(names(errors), bioc_packages)
    if (attempt > 1 || !length(broken)) break
    message("Reinstalling packages that could not be loaded: ", paste(broken, collapse = ", "))
    options(BioC_mirror = "https://bioconductor.posit.co")
    Sys.setenv(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = "TRUE", BIOCONDUCTOR_CONFIG_FILE = "https://bioconductor.posit.co/config.yaml")
    for (library_path in .libPaths()) try(suppressWarnings(remove.packages(broken, lib = library_path)), silent = TRUE)
    try(BiocManager::install(broken, ask = FALSE, update = FALSE, force = TRUE))
    next
  }
  message("Installing missing dependencies: ", paste(missing_dependencies, collapse = ", "))
  options(BioC_mirror = "https://bioconductor.posit.co")
  Sys.setenv(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = "TRUE", BIOCONDUCTOR_CONFIG_FILE = "https://bioconductor.posit.co/config.yaml")
  try(BiocManager::install(missing_dependencies, ask = FALSE, update = FALSE))
}
if (length(errors)) stop("These installed packages could not be loaded:\n", paste0("  ", names(errors), ": ", errors, collapse = "\n"),
                         "\nSend this message to the workshop organiser, or reinstall the package named in it.")
message("Required packages are installed and can be loaded.")
