import pandas as pd
import plotly.graph_objects as go
import plotly.express as px

data = pd.read_csv('191209-2-1_results.csv')
info = pd.read_csv('sample_info.csv')
merged = data.merge(info, on='Accession')

fig_t7 = px.boxplot(merged, x='sample')