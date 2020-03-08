import pandas as pd

info = pd.read_csv('Summary_Accuracy_200226.csv', dtype={'Batch': str})
batchids = info['Batch'].unique().tolist()
print(",".join(batchids))
