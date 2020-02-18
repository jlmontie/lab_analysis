import pandas as pd
import plotly.express as px

df = pd.read_excel('ecoli_results_annotated.xlsx')
print(df['source'].unique())
fig = px.box(df, y='phage_count_ratio', x='phage', color='source', points='all',
             labels={'phage_count_ratio': 'Normalized Ecoli reads / Normalized phage reads'},
             log_y=True, hover_data=['phage_count_ratio', 'phage', 'source', 'a549', 'accession'],
             category_orders={'source': ['Microbiologics', 'ViraPur', 'Zeptometrix', 'ATCC', 'CPT', 'CPT-2nd batch', 'CPT-2nd batch (repeat)'],
                              'phage': ['MS2', 'Qbeta', 'T7', 'T4']})
fig.write_html('comparison_analysis_plot.html')

df_genomic = df[~df['genome_coverage'].isna() & ~(df['genome_coverage'] == 0)]
fig2 = px.box(df_genomic, y='genome_coverage', x='phage', color='source', points='all',
              labels={'phage_count_ratio': 'Normalized Ecoli reads / Normalized phage reads'},
              log_y=True, hover_data=['phage_count_ratio', 'phage', 'source', 'a549', 'accession'])
fig2.write_html('comparison_analysis_plot_genome_coverage.html')

fig3 = px.box(df, y='16s_coverage', x='phage', color='source', points='all',
              labels={'phage_count_ratio': 'Normalized Ecoli reads / Normalized phage reads'},
              log_y=True, hover_data=['phage_count_ratio', 'phage', 'source', 'a549', 'accession'])
fig3.write_html('comparison_analysis_plot_16s_coverage.html')

df_fill_small = df.fillna(value=1e-6)
fig4 = px.box(df_fill_small, y='phage_count_ratio', x='phage', color='source', points='all',
             labels={'phage_count_ratio': 'Normalized Ecoli reads / Normalized phage reads'},
             log_y=True, hover_data=['phage_count_ratio', 'phage', 'source', 'a549', 'accession'])
fig4.write_html('comparison_analysis_plot_with_zeros.html')
df_fill_zeros = df.fillna(value=0)
fig5 = px.box(df_fill_zeros, y='16s_coverage', x='phage', color='source', points='all',
              labels={'phage_count_ratio': 'Normalized Ecoli reads / Normalized phage reads'},
              hover_data=['phage_count_ratio', 'phage', 'source', 'a549', 'accession'])
fig5.write_html('comparison_analysis_plot_16s_coverage_with_zeros.html')