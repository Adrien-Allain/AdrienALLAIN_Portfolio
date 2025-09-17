#======= Importation des bibliothèques =======#
import pandas as pd
import plotly.express as px
import dash
from dash import dcc, html, Input, Output
import dash_bootstrap_components as dbc
import datetime

#======= Chargement des données et modifications de la base =======#
df = pd.read_csv("supermarket_sales.csv")

# Convertir la colonne Date en format datetime
df['Date'] = pd.to_datetime(df['Date'])

# Créer une colonne pour la semaine
df['Week'] = df['Date'].dt.isocalendar().week
df['Year_Week'] = df['Date'].dt.strftime('%Y-%V')

#======= Initilisation de l'application Dash avec bootstrap =======#

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], title="Tableau de Bord des Ventes de Supermarché")

server = app.server


# Options pour les filtres
options_de_genre = [{'label': 'Tous', 'value': 'All'}] + [{'label': 'Femme', 'value': 'Female'}, {'label': 'Homme', 'value': 'Male'}]
options_de_villes = [{'label': 'Toutes', 'value': 'All'}] + [{'label': 'Yangon', 'value': 'Yangon'}, {'label': 'Naypyitaw', 'value': 'Naypyitaw'}, {'label': 'Mandalay', 'value': 'Mandalay'}]


    
#======= Maquette de l'application =======#

app.layout = dbc.Container([
    dbc.Row([
        dbc.Col([
            html.H1("Tableau de Bord des Ventes de Supermarché", className="text-center my-4")
        ])
    ], style={"background-color": "#4f75aa", "color": "white", "padding": "10px", "margin-bottom": "30px", "borderRadius": "15px"}),

    # Filtres
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H4("Filtres", className="text-center"),
                               style={"backgroundColor": "#cce5ff", "color": "#083c66"}),
                dbc.CardBody([
                    dbc.Row([
                        dbc.Col([
                            dbc.Label("Genre"),
                            dcc.Dropdown(
                                id='filtre-genre',
                                options=options_de_genre,
                                value='All',
                                clearable=False
                            ),
                        ], md=6),

                        dbc.Col([
                            dbc.Label("Ville"),
                            dcc.Dropdown(
                                id='filtre-ville',
                                options=options_de_villes,
                                value='All',
                                clearable=False
                            ),
                        ], md=6),
                    ])
                ])
            ])
        ])
    ], className="mb-4"),

    # Indicateurs
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H4("Montant total des achats", className="text-center"),
                               style={"backgroundColor": "#cce5ff", "color": "#083c66"}),
                dbc.CardBody([
                    html.H3(id='montant-total', className="text-center text-primary")
                ])
            ], className="h-100")
        ], md=6),

        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H4("Nombre total d'achats", className="text-center"),
                               style={"backgroundColor": "#cce5ff", "color": "#083c66"}),
                dbc.CardBody([
                    html.H3(id='total-achats', className="text-center text-primary")
                ])
            ], className="h-100")
        ], md=6),
    ], className="mb-4"),

    # histogramme
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H4("Répartition des montants totaux des achats", className="text-center"),
                               style={"backgroundColor": "#cce5ff", "color": "#083c66"}),
                dbc.CardBody([
                    dcc.Graph(id='histogramme')
                ])
            ])
        ], md=12),
    ], className="mb-4"),

    # Camembert et ligne sur la même ligne 
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H4("Répartition des catégories de produits", className="text-center"),
                               style={"backgroundColor": "#cce5ff", "color": "#083c66"}),
                dbc.CardBody([
                    dcc.Graph(id='camembert')
                ])
            ], className="h-100")
        ], md=6),

        dbc.Col([
            dbc.Card([
                dbc.CardHeader(html.H4("Évolution des achats par semaine", className="text-center"),
                               style={"backgroundColor": "#cce5ff", "color": "#083c66"}),
                dbc.CardBody([
                    dcc.Graph(id='plot-evolution')
                ])
            ], className="h-100")
        ], md=6),
    ]),

], fluid=True, className="p-4", style={"backgroundColor": "#4f75aa"}) # Couleur de fond 


def filter_data(gender, city):
    filtered_df = df.copy()
    
    if gender != 'All': #Si gender n'est pas all on prend le genre sélectionné dans le filtre
        filtered_df = filtered_df[filtered_df['Gender'] == gender]

    if city != 'All': #Même chose que pour gender
        filtered_df = filtered_df[filtered_df['City'] == city]

    return filtered_df

#======= Callbacks =======#
@app.callback(
    [Output('montant-total', 'children'),
     Output('total-achats', 'children'),
     Output('histogramme', 'figure'),
     Output('camembert', 'figure'),
     Output('plot-evolution', 'figure')],
    [Input('filtre-genre', 'value'), #Input pour les deux filtres
     Input('filtre-ville', 'value')]
)

def update_dashboard(gender, city):
    filtered_df = filter_data(gender, city)

    # Calcul des indicateurs
    montant_total = f"{filtered_df['Total'].sum():.2f} €" #somme de la colonne "Total" je mets .2f pour afficher 2 chiffres après la virgule
    achats_totaux = len(filtered_df)

    #Dictionnaire pour attribuer couleur aux villes
    couleurs_villes = {
    'Yangon': 'darkblue',
    'Mandalay': 'blue',
    'Naypyitaw': 'lightblue',
    }
    # Histogramme
    if gender == 'All' and city == 'All':
        hist_fig = px.histogram(
            filtered_df,
            x='Total',
            color='City',
            nbins=20,
            title="Répartition des montants totaux des achats par ville",
            labels={'Total': 'Montant total', 'City': 'Ville'},
            barmode='stack',
            color_discrete_map=couleurs_villes,
            template='plotly_white'
        )
    else:
        hist_fig = px.histogram(
            filtered_df,
            x='Total',
            nbins=20,
            title="Répartition des montants totaux des achats",
            labels={'Total': 'Montant total'},
            color= "City",
            color_discrete_map=couleurs_villes,
            template='plotly_white'
        )

        filter_info = [] #Création de la liste infos_filtres pour rendre titre graphique dynamique 
        if gender != 'All':
            filter_info.append(f"Sexe: {gender}") #Ajoute le genre sélectionné à la liste si pas "all"
        if city != 'All':
            filter_info.append(f"Ville: {city}")  #Ajoute la ville sélectionnée à la liste si pas "all"

        if filter_info:
            hist_fig.update_layout(
                title=f"Répartition des montants totaux des achats ({', '.join(filter_info)})" # La légende change en fonction des filtres 
            )

    hist_fig.update_layout(
        margin=dict(l=20, r=20, t=50, b=20),
        bargap=0.1,
        xaxis_title="Montant total",
        yaxis_title="Nombre d'achats",
        legend_title="Ville"
    )

    # Diagramme circulaire
    camembert_fig = px.pie(
        filtered_df,
        names='Product line',
        title="Répartition des catégories de produits",
        color_discrete_sequence=px.colors.sequential.Blues
    )
    filter_info = [] #Création de la liste infos_filtres pour rendre titre graphique dynamique 
    if gender != 'All':
        filter_info.append(f"Sexe: {gender}") #Ajoute le genre sélectionné à la liste si pas "all"
    if city != 'All':
        filter_info.append(f"Ville: {city}")  #Ajoute la ville sélectionnée à la liste si pas "all"

    if filter_info:
        camembert_fig.update_layout(
            title=f"Évolution des achats par semaine ({', '.join(filter_info)})" # La légende change en fonction des filtres 
            )
    camembert_fig.update_legends(title="Catégorie de produit")
    camembert_fig.update_layout(margin=dict(l=20, r=20, t=50, b=20))
    
    # courbe évolution temporelle
    if city == 'All':
        weekly_data = filtered_df.groupby(['Year_Week', 'City'])['Total'].sum().reset_index()
        courbe_temp_fig = px.line(
            weekly_data,
            x='Year_Week',
            y='Total',
            color='City',
            title="Évolution des achats par semaine",
            color_discrete_map=couleurs_villes,
            template='plotly_white'  
        )
    else:
        weekly_data = filtered_df.groupby(['Year_Week'])['Total'].sum().reset_index()
        weekly_data['City'] = city
        courbe_temp_fig = px.line(
            weekly_data,
            x='Year_Week',
            y='Total',
            color='City',
            title=f"Évolution des achats par semaine",
            color_discrete_map=couleurs_villes,
            template='plotly_white'  
        )
        
    filter_info = [] #Création de la liste infos_filtres pour rendre titre graphique dynamique 
    if gender != 'All':
        filter_info.append(f"Sexe: {gender}") #Ajoute le genre sélectionné à la liste si pas "all"
    if city != 'All':
        filter_info.append(f"Ville: {city}")  #Ajoute la ville sélectionnée à la liste si pas "all"

    if filter_info:
        courbe_temp_fig.update_layout(
            title=f"Évolution des achats par semaine ({', '.join(filter_info)})" # La légende change en fonction des filtres 
            )
        
    courbe_temp_fig.update_xaxes(title="Semaine")
    courbe_temp_fig.update_yaxes(title="Montant total des achats")
    courbe_temp_fig.update_layout(margin=dict(l=20, r=20, t=50, b=20))
    courbe_temp_fig.update_legends(title="Ville")

    return montant_total, achats_totaux, hist_fig, camembert_fig, courbe_temp_fig

# Lancer l'application
if __name__ == '__main__':
    app.run_server(debug=True, port=8067)