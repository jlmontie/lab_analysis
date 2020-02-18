import plotly.graph_objects as go
import pandas as pd
import numpy as np

df_dict = pd.read_excel('heatmap_data.xlsx', sheet_name=None, header=[0, 1], index_col=0)
for name, df in df_dict.items():
    accessions = list(df.columns.get_level_values(1))
    df.columns = df.columns.droplevel(level=1)
    df = df * 100
    customdata = [accessions for i in range(df.shape[0])]
    fig = go.Figure()
    fig.add_trace(
        go.Heatmap(
            z=df.values,
            y=df.index,
            hovertemplate='Accession: %{customdata}<br>Metric: %{y}<br>Value: %{z}<extra></extra>',
            customdata=customdata
        )
    )
    fig.update_layout(title=name.capitalize())
    fig.update_yaxes(autorange='reversed')
    fig.update_xaxes(tickmode='array', tickvals=np.arange(df.shape[1])[::2] + 0.5, ticktext=df.columns.values[::2])
    fig.write_html(f'plots/{name}.html')
