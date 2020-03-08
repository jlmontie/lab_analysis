import pandas as pd

arup = pd.read_csv('arup_coverages_with_review.csv')
synergy = pd.read_csv('synergy_coverages_with_review.csv')
arup['source'] = 'ARUP'
synergy['source'] = 'Synergy'
concat = pd.concat([arup, synergy])
concat = concat.drop(columns='reporting_id')
print(f"arup.shape: {arup.shape}")
print(f"synergy.shape: {synergy.shape}")
print(f"concat.shape: {concat.shape}")
print(f"concat review counts:\n{concat['review'].value_counts()}")
concat.to_csv('coverages_all.csv', index=False)
