

def create_network_dashboard(prioritized_markers, threshold=0.5, top_n=75):
    """
    Creates an advanced interactive network dashboard visualization showing relationships
    between genera, proteins, and mechanisms with multiple viewing options and filters.

    Args:  prioritized_markers (DataFrame): DataFrame containing marker data
           threshold (float): Initial threshold for including edges (0-1)
           top_n (int): Number of top markers to include

    Returns:   dash.Dash application that can be run with app.run()
    """
    import networkx as nx
    import pandas as pd
    import numpy as np
    from dash import Dash, dcc, html, Input, Output, callback
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import dash_bootstrap_components as dbc
    from community import community_louvain
    import plotly.express as px

    # Style definitions
    colors = {
        'genus': '#2196F3',      # Blue
        'protein': '#4CAF50',    # Green
        'mechanism': '#FF5722',  # Orange-red
        'background': '#F7F7F7', # Light gray
        'text': '#333333'        # Dark gray
    }

    # Check required libraries
    try:
        from community import community_louvain
    except ImportError:
        print("Please install python-louvain: pip install python-louvain")

    # Process data and create network
    def process_network_data(markers_df, thresh, max_markers):
        # Check if the required columns exist
        required_cols = ['Genus', 'protein_name', 'combined_score', 'corrosion_mechanisms']
        missing_cols = [col for col in required_cols if col not in markers_df.columns]

        if missing_cols:
            print(f"Warning: Missing required columns: {missing_cols}")
            print("Available columns:", markers_df.columns.tolist())
            # Try to find alternative column names
            column_mapping = {}

            # Map for possible alternative column names
            alternatives = {
                'Genus': ['genus', 'GENUS', 'Genera', 'genera'],
                'protein_name': ['protein', 'Protein', 'protein_id', 'ProteinName', 'name'],
                'combined_score': ['score', 'total_score', 'Score', 'TotalScore'],
                'corrosion_mechanisms': ['mechanisms', 'Mechanisms', 'corrosion_mechanism']
            }

            # Try to find alternatives
            for missing in missing_cols:
                for alt in alternatives.get(missing, []):
                    if alt in markers_df.columns:
                        column_mapping[missing] = alt
                        print(f"Using '{alt}' instead of '{missing}'")
                        break

            # If we couldn't find alternatives for all missing columns, return empty data
            if len(column_mapping) < len(missing_cols):
                missing_after_mapping = [col for col in missing_cols if col not in column_mapping]
                print(f"Could not find alternatives for: {missing_after_mapping}")
                return None, None, None, None, None

            # Create a view with renamed columns
            markers_df = markers_df.rename(columns=column_mapping)

        # Create a networkx graph
        G = nx.Graph()

        # Node type tracking
        genera = set()
        proteins = set()
        mechanisms = set()

        # Node attributes for tracking info
        node_types = {}
        node_info = {}

        # Extract top markers
        try:
            top_markers = markers_df.sort_values('combined_score', ascending=False).head(max_markers)
        except Exception as e:
            print(f"Error sorting markers: {e}")
            # Try without sorting if there was an error
            top_markers = markers_df.head(max_markers)

        # Track edge weights for filtering
        edge_weights = {}

        # Process each row
        for _, row in top_markers.iterrows():
            try:
                genus = str(row['Genus'])
                protein = str(row['protein_name'])

                if pd.isna(genus) or pd.isna(protein):
                    continue

                # Trim long protein names for better visualization
                protein_display = protein
                if len(protein) > 30:
                    protein_display = protein[:27] + "..."

                # Store the node information
                genera.add(genus)
                proteins.add(protein_display)

                node_types[genus] = 'genus'
                node_info[genus] = {
                    'type': 'genus',
                    'full_name': genus,
                    'proteins': [] if genus not in node_info else node_info[genus]['proteins'],
                    'score': float(row.get('combined_score', 1.0)) if isinstance(row.get('combined_score', 1.0), (int, float)) else 1.0
                }
                node_info[genus]['proteins'].append(protein_display)

                node_types[protein_display] = 'protein'
                node_info[protein_display] = {
                    'type': 'protein',
                    'full_name': protein,
                    'genus': genus,
                    'mechanisms': [] if protein_display not in node_info else node_info[protein_display]['mechanisms'],
                    'score': float(row.get('combined_score', 1.0)) if isinstance(row.get('combined_score', 1.0), (int, float)) else 1.0
                }

                # Add genus-protein edge with normalized weight
                score = row.get('combined_score', 1.0)
                # Handle non-numeric scores
                if not isinstance(score, (int, float)) or pd.isna(score):
                    score = 1.0

                edge_key = (genus, protein_display)
                edge_weights[edge_key] = float(score)

                # Add mechanism nodes and edges
                mechanisms_col = row.get('corrosion_mechanisms', '')
                if isinstance(mechanisms_col, str) and mechanisms_col:
                    for mech in mechanisms_col.split(';'):
                        mech = mech.strip()
                        if mech:  # Only add non-empty mechanisms
                            mechanisms.add(mech)
                            node_types[mech] = 'mechanism'
                            if mech not in node_info:
                                node_info[mech] = {
                                    'type': 'mechanism',
                                    'full_name': mech,
                                    'proteins': [],
                                    'score': 1.0
                                }
                            node_info[mech]['proteins'].append(protein_display)
                            node_info[protein_display]['mechanisms'].append(mech)

                            edge_key = (protein_display, mech)
                            edge_weights[edge_key] = 1.0
            except Exception as e:
                print(f"Error processing row: {e}")
                continue

        # Normalize edge weights to 0-1 scale for filtering
        if edge_weights:
            max_weight = max(edge_weights.values())
            min_weight = min(edge_weights.values())
            weight_range = max_weight - min_weight

            if weight_range > 0:
                for edge, weight in edge_weights.items():
                    edge_weights[edge] = (weight - min_weight) / weight_range

            # Add edges that meet the threshold
            for edge, weight in edge_weights.items():
                if weight >= thresh:
                    G.add_edge(edge[0], edge[1], weight=weight)

        # Check if we have a valid graph
        if len(G.nodes()) == 0:
            print("No nodes in graph. Check your data or lower the threshold.")
            return None, None, None, None, None

        # Add node attributes to graph
        nx.set_node_attributes(G, node_types, 'type')

        # Add info attributes
        for node, info in node_info.items():
            if node in G.nodes():
                for key, value in info.items():
                    G.nodes[node][key] = value

        # Calculate community structure
        partition = community_louvain.best_partition(G)
        nx.set_node_attributes(G, partition, 'community')

        return G, node_types, node_info, genera, proteins, mechanisms

    # Create initial graph
    G, node_types, node_info, genera, proteins, mechanisms = process_network_data(
        prioritized_markers, threshold, top_n
    )

    if G is None:
        app = Dash(__name__)
        app.layout = html.Div([
            html.H4("Error: Could not create network visualization"),
            html.P("Check your data format and column names.")
        ])
        return app

    # Create main layout for the app
    app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

    app.layout = dbc.Container([
        dbc.Row([
            dbc.Col([
                html.H2("Corrosion Mechanism Network Explorer",
                        style={'textAlign': 'center', 'marginTop': '10px', 'color': colors['text']}),
                html.Hr()
            ], width=12)
        ]),

        dbc.Row([
            # Controls panel
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Visualization Controls"),
                    dbc.CardBody([
                        html.H6("Network Filtering"),
                        dbc.Label("Edge Weight Threshold:"),
                        dcc.Slider(
                            id='threshold-slider',
                            min=0.1,
                            max=0.9,
                            step=0.1,
                            value=threshold,
                            marks={i/10: str(i/10) for i in range(1, 10)},
                        ),
                        html.Br(),

                        dbc.Label("Node Types:"),
                        dbc.Checklist(
                            id='node-type-checklist',
                            options=[
                                {'label': 'Genera', 'value': 'genus'},
                                {'label': 'Proteins', 'value': 'protein'},
                                {'label': 'Mechanisms', 'value': 'mechanism'}
                            ],
                            value=['genus', 'protein', 'mechanism'],
                            inline=True
                        ),
                        html.Br(),

                        html.H6("Layout Options"),
                        dbc.RadioItems(
                            id='layout-type',
                            options=[
                                {'label': '2D Force-directed', 'value': '2d_force'},
                                {'label': '3D Force-directed', 'value': '3d_force'},
                                {'label': 'Circular', 'value': 'circular'},
                                {'label': 'Community-based', 'value': 'community'}
                            ],
                            value='2d_force'
                        ),
                        html.Br(),

                        html.H6("Node Sizing"),
                        dbc.RadioItems(
                            id='node-size-option',
                            options=[
                                {'label': 'By Connectivity', 'value': 'degree'},
                                {'label': 'By Importance (PageRank)', 'value': 'pagerank'},
                                {'label': 'By Score', 'value': 'score'},
                                {'label': 'Uniform', 'value': 'uniform'}
                            ],
                            value='degree'
                        ),
                        html.Br(),

                        html.H6("Node Search"),
                        dcc.Dropdown(
                            id='node-search',
                            options=[{'label': node, 'value': node} for node in sorted(list(G.nodes()))],
                            placeholder="Search for a node...",
                        ),
                        html.Br(),

                        dbc.Button("Reset View", id="reset-button", color="primary", className="mt-2"),
                    ])
                ], style={'height': '100%'}),
            ], width=3),

            # Main visualization panel
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Network Visualization"),
                    dbc.CardBody([
                        dcc.Graph(
                            id='network-graph',
                            figure={},
                            style={'height': '70vh'}
                        )
                    ])
                ]),
            ], width=9),
        ]),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Selected Node Details"),
                    dbc.CardBody([
                        html.Div(id='node-details', children="Click on a node to see details")
                    ])
                ]),
            ], width=12),
        ], className="mt-3"),

        dbc.Row([
            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Network Statistics"),
                    dbc.CardBody([
                        html.Div(id='network-stats')
                    ])
                ]),
            ], width=6),

            dbc.Col([
                dbc.Card([
                    dbc.CardHeader("Community Distribution"),
                    dbc.CardBody([
                        dcc.Graph(
                            id='community-chart',
                            figure={},
                            style={'height': '30vh'}
                        )
                    ])
                ]),
            ], width=6),
        ], className="mt-3 mb-5"),

    ], fluid=True)

    # Callbacks for interactivity
    @app.callback(
        [Output('network-graph', 'figure'),
         Output('community-chart', 'figure'),
         Output('network-stats', 'children')],
        [Input('threshold-slider', 'value'),
         Input('node-type-checklist', 'value'),
         Input('layout-type', 'value'),
         Input('node-size-option', 'value'),
         Input('node-search', 'value'),
         Input('reset-button', 'n_clicks')]
    )
    def update_graph(threshold_value, selected_node_types, layout_type, node_size_option, search_node, reset_clicks):
        # Recreate the graph with new threshold
        G, node_types, node_info, genera, proteins, mechanisms = process_network_data(
            prioritized_markers, threshold_value, top_n
        )

        if G is None:
            empty_fig = go.Figure()
            empty_fig.add_annotation(text="No data available with current settings",
                                    showarrow=False, font=dict(size=16))
            empty_bar = go.Figure()
            empty_bar.add_annotation(text="No communities to display",
                                    showarrow=False, font=dict(size=16))
            return empty_fig, empty_bar, "No network data available"

        # Filter by node types
        nodes_to_include = [node for node in G.nodes() if G.nodes[node]['type'] in selected_node_types]
        G_filtered = G.subgraph(nodes_to_include)

        if len(G_filtered.nodes()) == 0:
            empty_fig = go.Figure()
            empty_fig.add_annotation(text="No nodes match the current filter settings",
                                    showarrow=False, font=dict(size=16))
            empty_bar = go.Figure()
            empty_bar.add_annotation(text="No communities to display",
                                    showarrow=False, font=dict(size=16))
            return empty_fig, empty_bar, "No nodes match the current filter settings"

        # Calculate layout
        if layout_type == '2d_force':
            pos = nx.spring_layout(G_filtered, k=0.3, iterations=50, seed=42)
            is_3d = False
        elif layout_type == '3d_force':
            pos = nx.spring_layout(G_filtered, k=0.3, iterations=50, seed=42, dim=3)
            is_3d = True
        elif layout_type == 'circular':
            pos = nx.circular_layout(G_filtered)
            is_3d = False
        elif layout_type == 'community':
            # Group by community
            pos = {}
            communities = {}
            for node in G_filtered.nodes():
                comm = G_filtered.nodes[node].get('community', 0)
                if comm not in communities:
                    communities[comm] = []
                communities[comm].append(node)

            # Position communities in a grid
            grid_size = int(np.ceil(np.sqrt(len(communities))))
            for i, (comm, nodes) in enumerate(communities.items()):
                # Position each community in a grid cell
                x_base = (i % grid_size) * 2
                y_base = (i // grid_size) * 2

                # Position nodes within the community in a circular layout
                comm_pos = nx.circular_layout(G_filtered.subgraph(nodes))

                # Adjust the positions to place in the grid
                for node, coords in comm_pos.items():
                    pos[node] = (coords[0] * 0.8 + x_base, coords[1] * 0.8 + y_base)

            is_3d = False

        # Calculate node metrics for sizing
        if node_size_option == 'degree':
            node_sizes = dict(G_filtered.degree())
            max_size = max(node_sizes.values()) if node_sizes else 1
            node_sizes = {node: (1 + 3 * size / max_size) for node, size in node_sizes.items()}
        elif node_size_option == 'pagerank':
            node_sizes = nx.pagerank(G_filtered)
            max_size = max(node_sizes.values()) if node_sizes else 0.01
            node_sizes = {node: (1 + 5 * size / max_size) for node, size in node_sizes.items()}
        elif node_size_option == 'score':
            node_sizes = {node: G_filtered.nodes[node].get('score', 1.0) for node in G_filtered.nodes()}
            max_size = max(node_sizes.values()) if node_sizes else 1
            node_sizes = {node: (1 + 3 * size / max_size) for node, size in node_sizes.items()}
        else:  # uniform
            node_sizes = {node: 1.0 for node in G_filtered.nodes()}

        # Define size multipliers
        size_multipliers = {
            'genus': 30,
            'protein': 20,
            'mechanism': 25
        }

        # Create node traces by type
        node_traces = {}

        for node_type, color in colors.items():
            if node_type in ['genus', 'protein', 'mechanism']:
                mode = 'markers+text' if node_type != 'protein' else 'markers'

                if is_3d:
                    node_traces[node_type] = go.Scatter3d(
                        x=[], y=[], z=[],
                        text=[],
                        mode=mode,
                        hoverinfo='text',
                        marker=dict(
                            color=color,
                            size=[],
                            line=dict(width=1, color='rgba(50, 50, 50, 0.8)')
                        ),
                        textfont=dict(
                            color='black'
                        ),
                        name=node_type.capitalize()
                    )
                else:
                    node_traces[node_type] = go.Scatter(
                        x=[], y=[],
                        text=[],
                        mode=mode,
                        hoverinfo='text',
                        marker=dict(
                            color=color,
                            size=[],
                            line=dict(width=1, color='rgba(50, 50, 50, 0.8)')
                        ),
                        textfont=dict(
                            color='black'
                        ),
                        name=node_type.capitalize()
                    )

        # Add nodes to traces with possible highlighting
        hover_texts = {}
        for node in G_filtered.nodes():
            node_type = G_filtered.nodes[node]['type']
            info = G_filtered.nodes[node]

            # Highlight search node if specified
            highlight = (search_node is not None and node == search_node)
            node_color = 'yellow' if highlight else colors[node_type]

            # Create hover text
            if node_type == 'genus':
                hover_texts[node] = f"<b>Genus:</b> {info.get('full_name', node)}<br>" +                                    f"<b>Connected proteins:</b> {len(info.get('proteins', []))}<br>" +                                    f"<b>Score:</b> {info.get('score', 'N/A'):.2f}<br>" +                                    f"<b>Community:</b> {info.get('community', 'N/A')}"
            elif node_type == 'protein':
                hover_texts[node] = f"<b>Protein:</b> {info.get('full_name', node)}<br>" +                                    f"<b>Genus:</b> {info.get('genus', 'Unknown')}<br>" +                                    f"<b>Mechanisms:</b> {', '.join(info.get('mechanisms', ['None']))}<br>" +                                    f"<b>Score:</b> {info.get('score', 'N/A'):.2f}<br>" +                                    f"<b>Community:</b> {info.get('community', 'N/A')}"
            else:  # mechanism
                hover_texts[node] = f"<b>Mechanism:</b> {info.get('full_name', node)}<br>" +                                    f"<b>Connected proteins:</b> {len(info.get('proteins', []))}<br>" +                                    f"<b>Community:</b> {info.get('community', 'N/A')}"

            # Calculate size
            size = size_multipliers[node_type] * node_sizes[node]

            # Add coordinates based on dimensionality
            if is_3d:
                if node in pos:
                    node_traces[node_type]['x'] += [pos[node][0]]
                    node_traces[node_type]['y'] += [pos[node][1]]
                    node_traces[node_type]['z'] += [pos[node][2]]
                    node_traces[node_type]['text'] += [node]
                    node_traces[node_type]['marker']['size'] += [size]
                    node_traces[node_type]['marker']['color'] = node_color if highlight else colors[node_type]
            else:
                if node in pos:
                    node_traces[node_type]['x'] += [pos[node][0]]
                    node_traces[node_type]['y'] += [pos[node][1]]
                    node_traces[node_type]['text'] += [node]
                    node_traces[node_type]['marker']['size'] += [size]
                    node_traces[node_type]['marker']['color'] = node_color if highlight else colors[node_type]

        # Create edge trace
        if is_3d:
            edge_trace = go.Scatter3d(
                x=[], y=[], z=[],
                line=dict(width=1, color='rgba(150, 150, 150, 0.5)'),
                hoverinfo='none',
                mode='lines',
                showlegend=False
            )
        else:
            edge_trace = go.Scatter(
                x=[], y=[],
                line=dict(width=1, color='rgba(150, 150, 150, 0.5)'),
                hoverinfo='none',
                mode='lines',
                showlegend=False
            )

        # Add edges
        for edge in G_filtered.edges():
            if edge[0] in pos and edge[1] in pos:
                x0, y0 = pos[edge[0]][0], pos[edge[0]][1]
                x1, y1 = pos[edge[1]][0], pos[edge[1]][1]

                if is_3d:
                    z0, z1 = pos[edge[0]][2], pos[edge[1]][2]
                    edge_trace['x'] += [x0, x1, None]
                    edge_trace['y'] += [y0, y1, None]
                    edge_trace['z'] += [z0, z1, None]
                else:
                    edge_trace['x'] += [x0, x1, None]
                    edge_trace['y'] += [y0, y1, None]

        # Create the main figure
        figure_traces = [edge_trace] + [trace for trace in node_traces.values() if len(trace['x']) > 0]

        if is_3d:
            fig = go.Figure(data=figure_traces)
        else:
            fig = go.Figure(data=figure_traces)

        # Update hover information
        for node_type in node_traces:
            if len(node_traces[node_type]['x']) > 0:
                node_hover_texts = []
                for node_text in node_traces[node_type]['text']:
                    node_hover_texts.append(hover_texts.get(node_text, node_text))
                node_traces[node_type].update(hovertext=node_hover_texts)

        # Layout
        if is_3d:
            fig.update_layout(
                title=f"Network of Genera, Proteins, and Corrosion Mechanisms<br><sup>Threshold: {threshold_value}</sup>",
                showlegend=True,
                hovermode='closest',
                margin=dict(b=20, l=5, r=5, t=40),
                scene=dict(
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    zaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                ),
                height=700,
                legend=dict(
                    yanchor="top",
                    y=0.99,
                    xanchor="right",
                    x=0.99
                )
            )
        else:
            fig.update_layout(
                title=f"Network of Genera, Proteins, and Corrosion Mechanisms<br><sup>Threshold: {threshold_value}</sup>",
                showlegend=True,
                hovermode='closest',
                margin=dict(b=20, l=5, r=5, t=40),
                xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                height=700,
                legend=dict(
                    yanchor="top",
                    y=0.99,
                    xanchor="right",
                    x=0.99
                )
            )

        # Create community distribution chart
        community_counts = {}
        for node in G_filtered.nodes():
            comm = G_filtered.nodes[node].get('community', 0)
            node_type = G_filtered.nodes[node]['type']

            if comm not in community_counts:
                community_counts[comm] = {'genus': 0, 'protein': 0, 'mechanism': 0}

            community_counts[comm][node_type] += 1

        # Convert to dataframe for plotting
        comm_df = []
        for comm, counts in community_counts.items():
            for node_type, count in counts.items():
                comm_df.append({
                    'Community': f"Community {comm}",
                    'Node Type': node_type.capitalize(),
                    'Count': count
                })

        comm_df = pd.DataFrame(comm_df)

        if len(comm_df) > 0:
            community_fig = px.bar(comm_df, x='Community', y='Count', color='Node Type',
                                  barmode='stack',
                                  color_discrete_map={
                                      'Genus': colors['genus'],
                                      'Protein': colors['protein'],
                                      'Mechanism': colors['mechanism']
                                  })
            community_fig.update_layout(
                title="Node Distribution by Community",
                xaxis_title="Community",
                yaxis_title="Number of Nodes",
                legend_title="Node Type"
            )
        else:
            community_fig = go.Figure()
            community_fig.add_annotation(text="No communities to display",
                                       showarrow=False, font=dict(size=16))

        # Generate network statistics
        n_nodes = len(G_filtered.nodes())
        n_edges = len(G_filtered.edges())
        n_communities = len(set(nx.get_node_attributes(G_filtered, 'community').values()))
        avg_degree = sum(dict(G_filtered.degree()).values()) / n_nodes if n_nodes > 0 else 0

        node_type_counts = {
            'genus': sum(1 for node in G_filtered.nodes() if G_filtered.nodes[node]['type'] == 'genus'),
            'protein': sum(1 for node in G_filtered.nodes() if G_filtered.nodes[node]['type'] == 'protein'),
            'mechanism': sum(1 for node in G_filtered.nodes() if G_filtered.nodes[node]['type'] == 'mechanism')
        }

        stats_html = [
            html.H5("Summary Statistics"),
            html.Ul([
                html.Li(f"Total Nodes: {n_nodes}"),
                html.Li(f"Total Edges: {n_edges}"),
                html.Li(f"Communities Detected: {n_communities}"),
                html.Li(f"Average Connections: {avg_degree:.2f}"),
                html.Li([
                    "Node Types: ",
                    html.Span(f"{node_type_counts['genus']} Genera, ", style={'color': colors['genus']}),
                    html.Span(f"{node_type_counts['protein']} Proteins, ", style={'color': colors['protein']}),
                    html.Span(f"{node_type_counts['mechanism']} Mechanisms", style={'color': colors['mechanism']})
                ])
            ])
        ]

        return fig, community_fig, stats_html
