import os

import pandas as pd

output_ls = []
for output in os.listdir('output'):
    df = pd.read_csv(os.path.join('output', output))
    output_ls.append(df)
pd.concat(output_ls).to_csv('merged_output.csv', index=False)
