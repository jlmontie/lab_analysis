import pandas as pd
info = pd.read_csv('synergy_clinical_sample_info.txt', sep='\t')
runids = info['seq_id'].unique().tolist()
print(",".join(runids))