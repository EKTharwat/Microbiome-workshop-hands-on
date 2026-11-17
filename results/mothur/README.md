# Mothur outputs

These files are from the completed run on the 30 anonymized workshop samples. The full analysis and code are in `workflow-results-side-by-side.html`; the side-by-side tables and figures are in `report_outputs/by_study/` at the repository root.

The shared tables contain feature counts by sample. The consensus taxonomy files provide the corresponding classifications; `final.fasta` and `final.count_table` contain processed sequences and their sample counts. `workshop.batch` records the commands, and `mothur-console.log` records the run.

Mothur warned that 13 unique sequences produced alignments that eliminated too many bases; their identifiers are in `reads.trim.contigs.good.unique.flip.accnos`. Interpret these selected panel data and downstream results as exploratory.
