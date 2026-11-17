# Completed mothur runs

The supplied material includes two complete expanded run folders. They contain multi-hundred-megabyte alignments and reference files that exceed GitHub’s normal file limit. Both run folders are preserved as ZIP downloads in the private repository release:

[Download the completed runs](https://github.com/EKTharwat/Microbiome-workshop-hands-on/releases/tag/supervisor-material-2026-11-17)

Download `teaching-20261115-212513.zip` and/or `teaching-20261116-183841.zip` and `SHA256SUMS` from the release. Each ZIP includes its `runs/` path. In the RStudio Console, set the working directory to the repository root and extract the downloaded ZIP:

```r
utils::unzip(file.choose(), exdir = getwd())
```

The raw paired-read archive, references, final mothur tables, downstream results, HTML reports and Word documents are included directly in the repository. The expanded runs are optional for starting a new workshop analysis.

SHA-256 checksums for the run downloads are included in `completed-runs-sha256.txt` and as a release asset.
