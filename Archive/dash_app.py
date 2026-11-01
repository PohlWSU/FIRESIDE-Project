# Import packages
from dash import Dash, html, dash_table, dcc, callback, Output, Input
import dash
import numpy as np
import pandas as pd
import plotly.express as px
#from galaxy_functions import ElementInput, redshifted_wavelengthFunc, line_flux_at_z, sensScaled, Detected_Galaxies

# Incorporate data
line_table_data = np.load('/data/projects/FIRESIDE/data/line_table_cat_v9.npy', allow_pickle=True) #not sure if needed

# Initialize the app
app = Dash()

Assigned_Column_Element = {
'[CII] 158':line_table_data[:11],
'[OI] 145':line_table_data[:12], 
'[NII] 122':line_table_data[:13],
'[OIII] 88':line_table_data[:14],
'[OI] 63':line_table_data[:15] ,
'[NIII] 57':line_table_data[:16], 
'[OIII] 52':line_table_data[:17], 
'[SIII] 19':line_table_data[:18] ,
'[NeIII] 16':line_table_data[:18] ,
'[NeII] 13':line_table_data[:19],
'[OIV] 26':line_table_data[:20] ,
'[SIII] 33':line_table_data[:21],
'[SII] 35':line_table_data[:22],
'[SIV] 11':line_table_data[:23] ,
'PAH 7.7':line_table_data[:24],
'PAH 11.3':line_table_data[:25], 
'Bralpha':line_table_data[:26],
'Hmalpha':line_table_data[:27] ,
'Pfalpa':line_table_data[:28] ,
'Halpha':line_table_data[:29],
'[NeV] 14':line_table_data[:30], 
'[NeV] 24':line_table_data[:31] 
}
# App layout
app.layout = [
    html.Div(children='N of Galaxies Detected for Selected Element'),
    html.Hr(),
    dcc.RadioItems(options=['[CII] 158', '[OI] 145', '[NII] 122', '[OIII] 88', '[OI] 63', 
    '[NIII] 57', '[OIII] 52', '[SIII] 19', '[NeIII] 16', '[NeII] 13',
   '[OIV] 26', '[SIII] 33', '[SII] 35', '[SIV] 11', 'PAH 7.7',
   'PAH 11.3', 'Bralpha', 'Hmalpha', 'Pfalpa', 'Halpha',
   '[NeV] 14', '[NeV] 24'], value='[CII] 158', id='controls-and-radio-item', labelStyle={'display': 'inline-block', 'margin-Right': '10px'}),
 
    dcc.Graph(figure={}, id='controls-and-graph')
]

# Add controls to build the interaction
@callback(
    Output(component_id='controls-and-graph', component_property='figure'),
    Input(component_id='controls-and-radio-item', component_property='value')
)
def update_graph(selected_Element):
    col_index=Assigned_Column_Element[selected_Elelemt]
    #Galaxy_Detections_Element=line_table_data[:, col_index] not actual working function, placeholder
    x_valHist=line_table_data[:, col_index]
    fig = px.histogram(x=x_valHist, nbins=50, title=f"Luminosities for {selected_Element}")
    return fig

# Run the app

if __name__ == '__main__':
    app.run(debug=True, port=8070)
