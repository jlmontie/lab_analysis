import json
import urllib
import os

import pandas as pd
import numpy as np
import dash
import dash_table
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
import plotly.express as px
from dash.dependencies import Input, Output


def generate_fig(threshold, df, col, cutoff=None):
    df_review_vals = df['review'].unique().tolist()
    fig = go.Figure()
    if not 'positive' in df_review_vals:
        fig.add_trace(go.Violin(y=df[col], x=[1] * len(df), box_visible=True))
    else:
        fig.add_trace(
            go.Violin(y=df.loc[~(df['review'] == 'positive'), col],
            x=[1] * len(df), side='negative', line_color='blue',
            pointpos=-1.5)
        )
        fig.add_trace(
            go.Violin(y=df.loc[df['review'] == 'positive', col],
            x=[1] * len(df), side='positive', line_color='orange',
            pointpos=1.5)
        )
    fig.update_traces(points='all')
    fig.add_trace(
        go.Scatter(
            y=[threshold] * 3,
            mode='lines',
            line_color='mediumseagreen'
        )
    )
    if cutoff is not None:
        fig.add_trace(
            go.Scatter(
                y=[cutoff] * 3,
                mode='lines',
                line_color='indianred'
            )
        )
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        showlegend=False,
        font=dict(size=14)
    )
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(range=[-0.1, 1.1])
    return fig


def generate_scatter(df, thresholds, cols, cutoffs, axes_titles):
    df_review_vals = df['review'].unique().tolist()
    fig = go.Figure()
    if not 'positive' in df_review_vals:
        fig.add_trace(go.Scatter(x=df[cols[0]], y=df[cols[1]]))
    else:
        fig.add_trace(
            go.Scatter(x=df.loc[~(df['review'] == 'positive'), cols[0]],
            y=df.loc[~(df['review'] == 'positive'), cols[1]], mode='markers',
            marker=dict(color='blue'))
        )
        fig.add_trace(
            go.Scatter(x=df.loc[df['review'] == 'positive', cols[0]],
            y=df.loc[df['review'] == 'positive', cols[1]], mode='markers',
            marker=dict(color='orange'))
        )
    fig.add_trace(
        go.Scatter(
            x=[thresholds[0], thresholds[0]],
            y=[-0.2, 1.2],
            mode='lines',
            line_color='mediumseagreen'
        )
    )
    fig.add_trace(
        go.Scatter(
            y=[thresholds[1], thresholds[1]],
            x=[-0.2, 1.2],
            mode='lines',
            line_color='mediumseagreen'
        )
    )
    fig.add_trace(
        go.Scatter(
            x=[cutoffs[0], cutoffs[0]],
            y=[-0.2, 1.2],
            mode='lines',
            line_color='indianred'
        )
    )
    if len(cutoffs) > 1:
        fig.add_trace(
            go.Scatter(
                y=[cutoffs[1], cutoffs[1]],
                x=[-0.2, 1.2],
                mode='lines',
                line_color='indianred'
            )
        )
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        showlegend=False,
        font=dict(size=14)
    )
    fig.update_xaxes(range=[-0.1, 1.1], title=axes_titles[0])
    fig.update_yaxes(range=[-0.1, 1.1], title=axes_titles[1])
    return fig


def generate_empty_fig():
    fig = go.Figure()
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        showlegend=False,
        font=dict(size=14)
    )
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(range=[-0.1, 1.1])
    return fig


def get_table(df_in, col, selected_threshold, cutoff=None):
    sources = df_in['source'].sort_values().unique().tolist()
    n = [np.sum(~df_in.loc[df_in['source'] == source, col].isna())
        for source in sources]
    n_above_selected = [np.sum(df_in.loc[df_in['source'] == source, col]
        >= selected_threshold) for source in sources]
    if cutoff is not None:
        n_above = [np.sum(df_in.loc[df_in['source'] == source, col] >= cutoff)
            for source in sources]
        df_return = pd.DataFrame({
            'Source': sources,
            'N': n,
            '>= Set Threshold': n_above_selected,
            'Cutoff': cutoff,
            '>= Cutoff': n_above
        })
    else:
        n_above = [np.nan] * len(sources)
        cutoff = [np.nan] * len(sources)
        df_return = pd.DataFrame({
            'Source': sources,
            'N': n,
            '>= Set Threshold': n_above_selected,
        })
    return df_return


def get_table_scatter(df_in, cols, selected_thresholds, cutoffs):
    sources = df_in['source'].sort_values().unique().tolist()
    n = [np.sum(np.all(~df_in.loc[df_in['source'] == source, cols].isna(),
            axis=1))
        for source in sources]
    n_above_selected = []
    n_above_cutoffs = []
    for source in sources:
        above_threshold_0 = df_in.loc[df_in['source'] == source,
                cols[0]] >= cutoffs[0]
        above_threshold_1 = df_in.loc[df_in['source'] == source,
                cols[1]] >= cutoffs[1]
        print(f"source col_0\n{df_in.loc[df_in['source'] == source, cols[0]]}")
        print(f"threshold 0\n{cutoffs[0]}")
        print(f"above threshold 0:\n{above_threshold_0}")
        print(f"source col_1\n{df_in.loc[df_in['source'] == source, cols[1]]}")
        print(f"threshold 1\n{cutoffs[1]}")
        print(f"above threshold 1:\n{above_threshold_1}")
        above = np.sum(np.all([
            df_in.loc[df_in['source'] == source,
                cols[0]] >= selected_thresholds[0],
            df_in.loc[df_in['source'] == source,
                cols[1]] >= selected_thresholds[1]
            ], axis=0)
        )
        n_above_selected.append(above)
        # if len(cutoffs) == 1:
        #     above_cutoff = np.sum(
        #         df_in.loc[df_in['source'] == source, cols[1]] >= cutoffs[0]
        #     )
        #     n_above_cutoffs.append(above_cutoff)
        # else:
        above_cutoff = np.sum(np.all([
            df_in.loc[df_in['source'] == source,
                cols[0]] >= cutoffs[0],
            df_in.loc[df_in['source'] == source,
                cols[1]] >= cutoffs[1]
            ], axis=0)
        )
        n_above_cutoffs.append(above_cutoff)
    df_return = pd.DataFrame({
        'Source': sources,
        'N': n,
        '>= Set Thresholds': n_above_selected,
        '>= Cutoffs': n_above_cutoffs
    })
    return df_return


def get_empty_table(core=False):
    if not core:
        df_return = pd.DataFrame({
            'Source': [np.nan],
            'N': [np.nan],
            # 'Threshold': [np.nan],
            '>= Set Threshold': [np.nan],
            'Cutoff': [np.nan],
            '>= Cutoff': [np.nan]
        })
    else:
        df_return = pd.DataFrame({
            'Source': [np.nan],
            'N': [np.nan],
            # 'Set Threshold': [np.nan],
            '>= Set Threshold': [np.nan],
        })
    return df_return


def get_empty_scatter_table():
    df_return = pd.DataFrame({
        'Source': [np.nan],
        'N': [np.nan],
        '>= Set Thresholds': [np.nan],
        '>= Cutoffs': [np.nan]
    })
    return df_return

df = pd.read_csv('data/562.csv',
    dtype={'accession': str, 'seq_sple': str, 'batch_id': str, 'seq_sple': str,
        'organism': str, 'taxid': int, 'review': str, 'above_cutoff': str,
        'source': str})

cutoffs = {}
with open('cutoffs.txt') as cutoff_file:
    for line in cutoff_file:
        obj = json.loads(line)
        for taxid in obj['taxids']:
            if not 'bacterial' in obj['subclass']:
                continue
            if 'rna_specific' in obj:
                rna_cutoff = obj['rna_specific']
            else:
                rna_cutoff = 0
            if 'dna_specific' in obj:
                dna_cutoff = obj['dna_specific']
            else:
                dna_cutoff = 0
            cutoffs.update({
                taxid:
                    {
                        'dna': dna_cutoff,
                        'rna': rna_cutoff
                    }
            })

uti_org_taxa = []
with open('org_taxids_uti.txt') as orgfile:
    orgfile.readline()
    for line in orgfile:
        data = line.strip().split('\t')
        uti_org_taxa.append((data[0], data[1]))

uti_options = [
    {'label': entry[0], 'value': entry[1]} for entry in uti_org_taxa
]

top_8_orgs = [
    'Escherichia coli', 'Klebsiella pneumoniae', 'Pseudomonas aeruginosa',
    'Proteus mirabilis', 'Enterococcus faecalis', 'Enterococcus faecium',
    'Staphylococcus aureus', 'Staphylococcus saprophyticus'
]
top_8_options = [
    {'label': entry[0], 'value': entry[1]} for entry in uti_org_taxa
    if entry[0] in top_8_orgs
]

all_org_taxa = []
with open('org_taxids_all.txt') as orgfile2:
    orgfile2.readline()
    for line in orgfile2:
        data = line.strip().split('\t')
        all_org_taxa.append((data[0], data[1]))

all_options = [
    {'label': entry[0], 'value': entry[1]} for entry in all_org_taxa
]

rna_table = get_table(df[df['source'] == 'ARUP'], 'rna_coverage', 0.95,
    cutoff=cutoffs[562]['rna'])
dna_table = get_table(df[df['source'] == 'ARUP'], 'dna_coverage', 0.5,
    cutoff=cutoffs[562]['dna'])
core_table = get_table(df[df['source'] == 'ARUP'], 'core_coverage', 0.5)
rna_v_dna_table = get_table_scatter(df[df['source'] == 'ARUP'],
    ['dna_coverage', 'rna_coverage'], [0.5, 0.95],
    [cutoffs[562]['dna'], cutoffs[562]['rna']])
rna_v_core_table = get_table_scatter(df[df['source'] == 'ARUP'],
    ['core_coverage', 'rna_coverage'], [0.5, 0.95],
    [0, cutoffs[562]['rna']])

navbar = dbc.NavbarSimple(
    # children=[
    #     dbc.NavItem(dbc.NavLink("Link", href="#")),
    #     dbc.DropdownMenu(
    #         nav=True,
    #         in_navbar=True,
    #         label="Menu",
    #         children=[
    #             dbc.DropdownMenuItem("Entry 1"),
    #             dbc.DropdownMenuItem("Entry 2"),
    #             dbc.DropdownMenuItem(divider=True),
    #             dbc.DropdownMenuItem("Entry 3"),
    #         ],
    #     ),
    # ],
    brand="Organism Coverages",
    brand_href="#",
    sticky="top",
)

body = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H5("Samples"),
                        dcc.Dropdown(
                            id='samples',
                            options=[
                                {'label': 'ARUP', 'value': 'ARUP'},
                                {'label': 'Synergy', 'value': 'Synergy'}
                            ],
                            value='ARUP',
                            multi=True
                        )
                    ]
                ),
                dbc.Col(
                    [
                        html.H5("Organism Pool"),
                        dcc.Dropdown(
                            id='org_pool',
                            options=[
                                {'label': 'Top 8 Uropathogens',
                                'value': 'top8'},
                                {'label': 'UTI Panel', 'value': 'uti'},
                                {'label': 'All', 'value': 'all'}
                            ],
                            value='top8'
                        )
                    ]
                ),
                dbc.Col(
                    [
                        html.H5("Organism"),
                        dcc.Dropdown(
                            id='organism',
                            options=top_8_options,
                            value=562
                        )
                    ]
                )
            ], style={'padding-top': 0}
        ),
        dbc.Row(
            [
                html.A(
                    'Download Data',
                    id='download-link',
                    download="coverage_data.csv",
                    href="",
                    target="_blank"
                ),
            ], style={'padding-top': 10}
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H4("RNA"),
                        dcc.Graph(
                            id='rna_fig',
                            figure=generate_fig(0.95, df, 'rna_coverage',
                                cutoff=cutoffs[562]['rna'])
                        ),
                    ]
                ),
                dbc.Col(
                    [
                        html.H4("DNA"),
                        dcc.Graph(
                            id='dna_fig',
                            figure=generate_fig(0.5, df, 'dna_coverage',
                                cutoff=cutoffs[562]['dna'])
                        ),
                    ]
                ),
                dbc.Col(
                    [
                        html.H4("Core"),
                        dcc.Graph(
                            id='core_fig',
                            figure=generate_fig(0.5, df, 'core_coverage')
                        ),
                    ]
                ),
            ], style={'padding-top': 20}
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        # dcc.Slider(
                        #     id='rna_slider',
                        #     min=0,
                        #     max=1,
                        #     step=0.0001,
                        #     value=0.95
                        # ),
                        html.H6("RNA Threshold"),
                        dcc.Input(
                            id='rna_input',
                            type='number',
                            min=0,
                            max=1,
                            step=0.0001,
                            value=0.95
                        )
                    ]
                ),
                dbc.Col(
                    [
                        # dcc.Slider(
                        #     id='dna_slider',
                        #     min=0,
                        #     max=1,
                        #     step=0.0001,
                        #     value=0.50
                        # ),
                        html.H6("DNA Threshold"),
                        dcc.Input(
                            id='dna_input',
                            type='number',
                            min=0,
                            max=1,
                            step=0.0001,
                            value=0.50
                        )
                    ]
                ),
                dbc.Col(
                    [
                        # dcc.Slider(
                        #     id='core_slider',
                        #     min=0,
                        #     max=1,
                        #     step=0.0001,
                        #     value=0.5
                        # ),
                        html.H6("Core Threshold"),
                        dcc.Input(
                            id='core_input',
                            type='number',
                            min=0,
                            max=1,
                            step=0.0001,
                            value=0.5
                        )
                    ]
                )
            ], style={'padding-top': 10}
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dash_table.DataTable(
                            id='rna_table',
                            columns=[{"name": i, "id": i} for i in
                                rna_table.columns],
                            data=rna_table.to_dict('records'),
                            style_cell={
                                'height': 'auto',
                                # all three widths are needed
                                'minWidth': '30px', 'width': '30px',
                                'maxWidth': '30px',
                                'whiteSpace': 'normal'
                            }
                        )
                    ], style={'padding-left': 20, 'padding-right': 20}
                ),
                dbc.Col(
                    [
                        dash_table.DataTable(
                            id='dna_table',
                            columns=[{"name": i, "id": i} for i in
                                dna_table.columns],
                            data=dna_table.to_dict('records'),
                            style_cell={
                                'height': 'auto',
                                # all three widths are needed
                                'minWidth': '30px', 'width': '30px',
                                'maxWidth': '30px',
                                'whiteSpace': 'normal'
                            }
                        )
                    ], style={'padding-left': 20, 'padding-right': 20}
                ),
                dbc.Col(
                    [
                        dash_table.DataTable(
                            id='core_table',
                            columns=[{"name": i, "id": i} for i in
                                core_table.columns],
                            data=core_table.to_dict('records'),
                            style_cell={
                                'height': 'auto',
                                # all three widths are needed
                                'minWidth': '28px', 'width': '28px',
                                'maxWidth': '28px',
                                'whiteSpace': 'normal'
                            }
                        )
                    ], style={'padding-left': 30, 'padding-right': 30}
                )
            ], style={'padding-top': 20}
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.H4("RNA vs DNA Coverage"),
                        dcc.Graph(
                            id='rna_v_dna',
                            figure=generate_scatter(df, [0.5, 0.95],
                                ['dna_coverage', 'rna_coverage'],
                                [cutoffs[562]['dna'],
                                    cutoffs[562]['rna']],
                                ['DNA Coverage', 'RNA Coverage'])
                        ),
                    ]
                ),
                dbc.Col(
                    [
                        html.H4("RNA vs Core Coverage"),
                        dcc.Graph(
                            id='rna_v_core',
                            figure=generate_scatter(df, [0.5, 0.95],
                                ['core_coverage', 'rna_coverage'],
                                [np.nan, cutoffs[562]['rna']],
                                ['Core Coverage', 'RNA Coverage'])
                        ),
                    ]
                )
            ], style={'padding-top': 20}
        ),
        dbc.Row(
            [
                dbc.Col(
                    [
                        dash_table.DataTable(
                            id='rna_v_dna_table',
                            columns=[{"name": i, "id": i} for i in
                                rna_v_dna_table.columns],
                            data=rna_v_dna_table.to_dict('records')
                        )
                    ], style={'padding-left': 20, 'padding-right': 20}
                ),
                dbc.Col(
                    [
                        dash_table.DataTable(
                            id='rna_v_core_table',
                            columns=[{"name": i, "id": i} for i in
                                rna_v_core_table.columns],
                            data=rna_v_core_table.to_dict('records')
                        )
                    ], style={'padding-left': 20, 'padding-right': 20}
                )
            ], style={'padding-top': 20, 'padding-bottom': 20}
        ),
        dcc.Store(
            id='store',
            data=df[df['source'] == 'ARUP'].to_dict('records')
        )
    ],
    className="mt-4",
)
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.FLATLY])
app.layout = html.Div([navbar, body])


# @app.callback(
#     Output('core_input', 'value'),
#     [Input('core_slider', 'value')]
# )
# def update_core_input(value):
#     return value

# @app.callback(
#     Output('dna_input', 'value'),
#     [Input('dna_slider', 'value')]
# )
# def update_dna_input(value):
#     return value

# @app.callback(
#     Output('rna_input', 'value'),
#     [Input('rna_slider', 'value')]
# )
# def update_rna_input(value):
#     return value

@app.callback(
    [Output('rna_fig', 'figure'),
    Output('dna_fig', 'figure'),
    Output('core_fig', 'figure'),
    Output('rna_v_dna', 'figure'),
    Output('rna_v_core', 'figure'),
    Output('rna_table', 'data'),
    Output('dna_table', 'data'),
    Output('core_table', 'data'),
    Output('rna_v_dna_table', 'data'),
    Output('rna_v_core_table', 'data'),
    Output('store', 'data')],
    [Input('rna_input', 'value'),
    Input('dna_input', 'value'),
    Input('core_input', 'value'),
    Input('samples', 'value'),
    Input('organism', 'value')]
)
def update_plots_and_tables(rna_input, dna_input, core_input, samples,
        organism):
    df_filepath = f"data/{organism}.csv"
    if os.path.isfile(df_filepath):
        df_filtered = pd.read_csv(f"data/{organism}.csv")
        if not isinstance(organism, int):
            organism = int(organism)
        if not isinstance(samples, list):
            samples = [samples]
        if samples:
            df_filtered = df_filtered[df_filtered['source'].isin(samples)]
        if organism in cutoffs:
            if 'rna' in cutoffs[organism]:
                rna_cutoff = cutoffs[organism]['rna']
            else:
                rna_cutoff = None
            if 'dna' in cutoffs[organism]:
                dna_cutoff = cutoffs[organism]['dna']
            else:
                dna_cutoff = None
        else:
            rna_cutoff = None
            dna_cutoff = None
        rna_fig = generate_fig(rna_input, df_filtered, 'rna_coverage',
            cutoff=rna_cutoff)
        dna_fig = generate_fig(dna_input, df_filtered, 'dna_coverage',
            cutoff=dna_cutoff)
        core_fig = generate_fig(core_input, df_filtered, 'core_coverage')
        rna_v_dna_fig = generate_scatter(df_filtered,
            [dna_input, rna_input], ['dna_coverage', 'rna_coverage'],
            [dna_cutoff, rna_cutoff], ['DNA Coverage', 'RNA Coverage'])
        rna_v_core_fig = generate_scatter(df_filtered,
            [core_input, rna_input], ['core_coverage', 'rna_coverage'],
            [np.nan, rna_cutoff], ['Core Coverage', 'RNA Coverage'])
        # Tables
        rna_table_df = get_table(df_filtered, 'rna_coverage', rna_input,
            cutoff=cutoffs[organism]['rna'])
        dna_table_df = get_table(df_filtered, 'dna_coverage', dna_input,
            cutoff=cutoffs[organism]['dna'])
        core_table_df = get_table(df_filtered, 'core_coverage', core_input)
        rna_v_dna_table = get_table_scatter(df_filtered,
            ['dna_coverage', 'rna_coverage'], [dna_input, rna_input],
            [cutoffs[organism]['dna'], cutoffs[organism]['rna']])
        rna_v_core_table = get_table_scatter(df_filtered,
            ['core_coverage', 'rna_coverage'], [core_input, rna_input],
            [0, cutoffs[organism]['rna']])
    else:
        rna_fig = generate_empty_fig()
        dna_fig = generate_empty_fig()
        core_fig = generate_empty_fig()
        rna_v_dna_fig = generate_empty_fig()
        rna_v_core_fig = generate_empty_fig()
        rna_table_df = get_empty_table()
        dna_table_df = get_empty_table()
        core_table_df = get_empty_table(core=True)
        rna_v_dna_table = get_empty_scatter_table()
        rna_v_core_table = get_empty_scatter_table()
        df_filtered = pd.DataFrame()
    return (rna_fig, dna_fig, core_fig, rna_v_dna_fig, rna_v_core_fig,
        rna_table_df.to_dict('records'), dna_table_df.to_dict('records'),
        core_table_df.to_dict('records'), rna_v_dna_table.to_dict('records'),
        rna_v_core_table.to_dict('records'), df_filtered.to_dict('records'))

@app.callback(
    Output('organism', 'options'),
    [Input('org_pool', 'value')]
)
def update_orgs(org_pool):
    if org_pool == 'top8':
        options = top_8_options
    elif org_pool == 'uti':
        options = uti_options
    elif org_pool == 'all':
        options = all_options
    return options

@app.callback(
    Output('download-link', 'href'),
    [Input('store', 'data')])
def update_download_link(df_dict_in):
    df_download = pd.DataFrame(df_dict_in)
    # dff = filter_data(filter_value)
    csv_string = df_download.to_csv(index=False, encoding='utf-8')
    csv_string = "data:text/csv;charset=utf-8," + urllib.parse.quote(csv_string)
    return csv_string

# @app.callback(
#     Output('core_slider', 'value'),
#     [Input('core_input', 'value')]
# )
# def update_core_slider(value):
#     return value

# @app.callback(
#     Output('dna_slider', 'value'),
#     [Input('dna_input', 'value')]
# )
# def update_dna_slider(value):
#     return value

# @app.callback(
#     Output('rna_slider', 'value'),
#     [Input('rna_input', 'value')]
# )
# def update_rna_slider(value):
#     return value

if __name__ == "__main__":
    app.run_server(port=9010)