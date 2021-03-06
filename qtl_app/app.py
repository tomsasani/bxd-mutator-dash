# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

import dash
import dash_core_components as dcc
import dash_html_components as html
import plotly.express as px
import plotly.graph_objects as go

import pandas as pd
import numpy as np
from dash.dependencies import Input, Output

app = dash.Dash(__name__)

gmap = pd.read_csv("data/bxd.gmap", skiprows=3)
pmap = pd.read_csv("data/bxd.pmap", skiprows=3)
geno = pd.read_csv("data/bxd.geno.new")

sbs_activity = pd.read_csv("data/COSMIC_SBS96_Activities_refit.txt", sep='\t')

muts = ("C>A", "C>T", "C>G", "A>T", "A>G", "A>C")

dfs = []

for mut in muts:
    for scan in ('rate', 'fraction'):
        formatted_mut = mut.replace('>', '.')
        lod_scores = pd.read_csv(f"data/{formatted_mut}.{scan}.lod.csv", skiprows=1, names=["marker", "lod"])
        lod_scores['mutation_type'] = mut
        lod_scores = lod_scores.merge(pmap, on="marker")
        lod_scores = lod_scores.query('chr != "X" and chr != "Y"')
        lod_scores['color'] = lod_scores['chr'].apply(lambda c: 'green' if int(c) % 2 == 0 else 'red')
        lod_scores['marker_number'] = np.arange(lod_scores.shape[0])
        lod_scores['scan_type'] = scan

        dfs.append(lod_scores)

df = pd.concat(dfs)

markers = list(gmap['marker'].values)

tidy_spectra = pd.read_csv("data/tidy_mutation_spectra.csv")

samples = tidy_spectra.bxd_strain_conv.unique()

# read in singletons
singleton = pd.read_csv("data/annotated_singleton_vars.csv")

app.layout = html.Div([
    dcc.Markdown("""
    ### A wild-derived antimutator drives germline mutation spectrum differences in a genetically diverse murine family

    Thomas A. Sasani, David G. Ashbrook, Annabel C. Beichman, Lu Lu, Abraham C. Palmer, Jonathan K. Pritchard, Robert W. Williams, Kelley Harris\
    """),
    dcc.Markdown("""###### **About this app**"""),
    dcc.Markdown(
        """We identified *de novo* germline mutations in a panel of recombinant inbred mouse lines (RILs) called the BXD, and mapped quantitative trait loci (QTL) for various phenotypes related to the mutation rate. You can find out more about the analysis in [our paper](https://www.biorxiv.org/content/10.1101/2021.03.12.435196v1.full)."""
    ),
    dcc.Markdown("""###### **You can use this app to:**"""),
    dcc.Markdown(""" 
    * Explore QTL mapping results
    * Visualize mutation fractions in mice with either haplotype at a particular locus
    * Identify the COSMIC mutation signatures active in a subset of BXD RILs


    """),
    html.Div([
        html.Div(
            className="row",
            children=[
                dcc.Markdown("""**Select a mutation type:**""",
                             style={'width': '33%'}),
                dcc.Markdown(
                    """**Select a phenotype to use for the QTL scan:**""",
                    style={'width': '33%'}),
                dcc.Markdown(
                    """**Select a marker below, or click on any point in the QTL plot:**""",
                    style={'width': '33%'}),
            ],
            style={'display': 'flex'},
        ),
        html.Div(className="row",
                 children=[
                     dcc.Dropdown(
                         id='mutation-type-dropdown',
                         options=[{
                             'label': i,
                             'value': i
                         } for i in muts],
                         value='C>A',
                         placeholder="Select a mutation type",
                         style={'width': '33%'},
                     ),
                     dcc.RadioItems(
                         id="scan-selection",
                         options=[
                             {
                                 'label': 'rate',
                                 'value': 'rate'
                             },
                             {
                                 'label': 'fraction',
                                 'value': 'fraction'
                             },
                         ],
                         value='fraction',
                         style={'width': '33%'},
                     ),
                     dcc.Dropdown(
                         id='marker-dropdown',
                         options=[{
                             'label': i,
                             'value': i
                         } for i in markers],
                         style={'width': '33%'},
                         placeholder="Select a marker",
                     ),
                 ],
                 style={'display': 'flex'}),
        html.Div(
            className="row",
            children=[
                dcc.Graph(
                    id='graph-with-dropdown',
                    hoverData={'points': [{
                        'customdata': 'rs31879829'
                    }]},
                    style={'width': '67%'},
                ),
                dcc.Graph(
                    id='mutation-stripplot',
                    style={'width': '33%'},
                ),
            ],
            style={'display': 'flex'},
        ),
    ]),
    html.Div(
        className="row",
        children=[
            dcc.Markdown(
                """**Select samples in panel B in which to visualize SBS activity.**"""
            ),
            dcc.Graph(
                id='sbs-activity-in-sample',
                hoverData={
                    'points': [
                        {
                            'customdata': [
                                'BXD001_TyJ_0361',
                                'BXD002_RwwJ_0430',
                                'BXD5_TyJ_0389',
                                'BXD006_TyJ_0474',
                                'BXD008_TyJ_0372',
                                'BXD009_TyJ_0383',
                            ],
                        },
                    ]
                },
            ),
        ],
    )
])

@app.callback(
    Output('graph-with-dropdown', 'figure'),
    [
        Input('mutation-type-dropdown', 'value'),
        Input('scan-selection', 'value'),
    ],
)
def update_qtl_plot(mutation_type, scan_type):
    filtered_df = df.query(f'mutation_type == "{mutation_type}" and scan_type == "{scan_type}"')
    chrom_to_max_idx = filtered_df.groupby('chr').max().to_dict()['marker_number']
    x_idxs, x_labels = [], []
    prev_idx = 0
    for chrom in np.arange(1, 20):
        max_idx = chrom_to_max_idx[str(chrom)]
        idx_to_use = (max_idx + prev_idx)  / 2
        x_idxs.append(idx_to_use)
        x_labels.append(str(chrom))
        prev_idx = max_idx

    fig = px.scatter(
        filtered_df,
        x="marker_number",
        y="lod",
        color="color",
        hover_name="marker",
        hover_data=["chr", "pos"],
        title="R/qtl2 mapping results (mm10/GRCm38)",
        labels={
            'marker_number': "Genotype marker on mm10",
            'lod': "LOD score",
        },
    )

    fig.update_layout(transition_duration=500,
                      showlegend=False,
                      xaxis=dict(
                          tickmode='array',
                          tickvals=x_idxs,
                          ticktext=x_labels,
                      ))

    return fig

@app.callback(
    Output('mutation-stripplot', 'figure'),
    [
        Input('graph-with-dropdown', 'clickData'),
        Input('mutation-type-dropdown', 'value'),
        Input('scan-selection', 'value'),
        Input('marker-dropdown', 'value'),
    ],
)
def update_stripplot(
    clickData,
    mutation_type,
    scan_type,
    marker_input,
):
    if clickData is None and marker_input is not None:
        marker_query = [marker_input]
    elif clickData is None and marker_input is None:
        marker_query = ["rs31879829"]
    else:
        marker_query = [clickData['points'][0]['hovertext']]
    genotypes_at_marker = geno.query(f'marker == {marker_query}')
    strain2hap = dict(zip(genotypes_at_marker, genotypes_at_marker.values[0, :]))

    spectra_mutation_type = tidy_spectra.query(f'base_mut == "{mutation_type}"')
    spectra_fraction = spectra_mutation_type.query(f'estimate_type == "{scan_type}"')
    spectra_with_hap = spectra_fraction.copy()
    spectra_with_hap['haplotype'] = spectra_with_hap['bxd_strain_conv'].apply(lambda b: strain2hap[b])

    fig = px.strip(
        spectra_with_hap,
        x="haplotype",
        y="estimate",
        color="haplotype",
        hover_name="bxd_strain_conv",
        title=f"Singleton {scan_type}s at {list(marker_query)[0]}",
        labels={
            "haplotype": "Haplotype at marker",
            "estimate": f"{mutation_type} {scan_type}"
        },
    )

    return fig

@app.callback(
    Output('sbs-activity-in-sample', 'figure'),
    [
        Input('mutation-stripplot', 'selectedData'),
    ],
)
def update_sbs_plot(selectedData):
    if selectedData is None:
        sample_names = list(np.random.choice(samples, 10))
    else:
        sample_names = [sd['hovertext'] for sd in selectedData['points']]
    sbs_activities = sbs_activity.query(f"Samples == {sample_names}")
    activity_tidy = sbs_activities.melt(
        id_vars=["Samples"],
        var_name=["signature"],
        value_name="mutnum",
    )
    fig = px.bar(
        activity_tidy,
        x="Samples",
        y="mutnum",
        color="signature",
        title="COSMIC mutation signature activity in selected samples",
        labels={"mutnum": "Number of singletons"},
    )
    return fig

if __name__ == '__main__':
    app.run_server(debug=True)
