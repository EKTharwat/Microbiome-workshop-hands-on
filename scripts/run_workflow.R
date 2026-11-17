# Render the step-by-step workflow in a separate R process.
if (!requireNamespace("rmarkdown", quietly = TRUE)) stop("Install the rmarkdown package first")
rmarkdown::render("workflow.Rmd", output_file = "workflow-results.html", quiet = FALSE)
