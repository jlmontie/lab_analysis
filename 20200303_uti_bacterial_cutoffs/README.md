# UTI Bacterial Cutoffs
Evaluating the cutoffs for bacterial organisms in the UTI test.

Utilizing two data sources
1. Synergy patient samples. `synergy_cli
2. ARUP samples run in SLC lab

Analysis steps
1. with `get_synergy_summary_data.py` script, iterate through rna bacterial
summary paths in `synergy_clinical_sample_info.txt`. Output the coverages
of all organisms to `synergy_coverages.csv`.
2. with `merge_synergy_md_review.py` script, merge MD review detection
status from the review portal dump
`explify-synergy-export-200303.xlsx` and output to
`synergy_coverages_with_review.csv`.

The coverages for organisms are compared in Plotly Dash webapp.
