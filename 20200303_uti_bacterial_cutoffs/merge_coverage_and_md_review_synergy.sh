SCRIPT=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/merge_coverage_and_md_review_synergy.py

COVERAGES1=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/synergy_coverages.csv
REVIEW1=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/explify-synergy-export-200303.xlsx
OUTPUT1=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/synergy_coverages_with_review.csv

COVERAGES2=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/synergy_coverages_fungpar.csv
REVIEW2=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/explify-synergy-export-200303.xlsx
OUTPUT2=/home/jmontgomery/lab_analysis/20200303_uti_bacterial_cutoffs/synergy_coverages_fungpar_with_review.csv

python $SCRIPT $COVERAGES1 $REVIEW1 $OUTPUT1
python $SCRIPT $COVERAGES2 $REVIEW2 $OUTPUT2