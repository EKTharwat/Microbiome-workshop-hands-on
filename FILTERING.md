# How the workshop reads were selected

The workshop contains 30 samples, with three samples per condition in each of five study groups. Samples were chosen partly because their reads matched the workshop's bacterial reference panel; they are not a random representation of the original studies.

Paired reads were assembled and screened for ambiguous bases, length (400–550 bases), and homopolymers (maximum eight). Unique assembled sequences were compared with the bundled full-length 16S reference sequences for the selected bacterial panel using VSEARCH on both strands. The search required 100% coverage of the assembled query sequence, not 100% coverage of the full-length reference. Four groups used a minimum identity of 99%; the remaining group used 97%. The threshold for every sample is recorded in `metadata/workshop_samples.tsv`.

Matches with near-equivalent assignments to different panel taxa (within 0.1 percentage point of the best identity) were excluded. Candidate sequences were then checked for chimeras with mothur's `chimera.vsearch` and detected candidate chimeras were removed. The matching original paired FASTQ records were recovered and their mate identifiers checked. Reads retained in the previous workshop dataset were preserved when its two samples in the 97% group were rechecked; all of those reads also passed that recheck.

The archive `vendor/workshop_reads.zip` contains the resulting neutral-named paired FASTQ files, assembled FASTA sequences, putative panel assignments, and filtering summaries. The workflow in `workflow.Rmd` starts with these filtered reads and performs the subsequent mothur and R analyses. The original unfiltered data are not included.

Identity to a 16S reference provides a putative taxonomic match, not certain species identification. Different identity thresholds and selection based on retained reads affect comparisons; the analyses are exploratory and between-group differences should not be attributed solely to biology.
