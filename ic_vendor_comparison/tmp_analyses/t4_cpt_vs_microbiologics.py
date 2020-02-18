import pandas as pd
import plotly.express as px

df = pd.read_excel('../ecoli_results_annotated_16s_read_count_for_RNA_phage.xlsx')
df = df[df['phage'] == 'T4']
fig = px.box(df, x='source', y='t4_normalized_reads', log_y=True)
fig.write_html('t4_cpt_vs_microbiologics.html')