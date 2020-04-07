SCRIPT=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/merge_coverage_and_md_review.py

COVERAGES=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/arup_coverages.csv
REVIEW=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/explify-synergy-dev-export-200303.xlsx
OUTPUT=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/arup_coverages_with_review.csv

# COVERAGES=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/arup_coverages_fungpar.csv
# REVIEW=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/explify-synergy-dev-export-200303.xlsx
# OUTPUT=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/arup_coverages_fungpar_with_review.csv

python $SCRIPT $COVERAGES $REVIEW $OUTPUT