import plotly.express as px
import pandas as pd

ms2 = pd.read_csv('../data/MS2_Ecoli_Coverage.csv')
qbeta = pd.read_csv('../data/Qbeta_Ecoli_Coverage.csv')

fig_ms2 = px.line(ms2, x='Position', y='Coverage', log_y=True)
fig_ms2.update_layout(title='MS2 <em>E. coli</em> Coverage')
fig_ms2.write_html('ms2_ecoli_coverage.html')

fig_qbeta = px.line(qbeta, x='Position', y='Coverage', log_y=True)
fig_qbeta.update_layout(title='Qbeta <em>E. coli</em> Coverage')
fig_qbeta.write_html('qbeta_ecoli_coverage.html')
