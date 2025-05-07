!pip install git+https://github.com/MagicAlex238/2_Micro.git

import plotly.graph_objects as go
import sys
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
from pathlib import Path
import ipywidgets as widgets
from IPython.display import display, HTML
#=====================================================================0000
# Own Scoring system
import corrosion_scoring as cs

# Determine the environment
'''if "google.colab" :
print("Running in Google Colab environment")
# for colab
base_dir = Path("/content/drive/MyDrive/MIC")
abundance_excel = base_dir / "data_picrust/merged_to_sequence.xlsx"
output_large = base_dir / "output_large"
output_base = base_dir
market_dir = base_dir / "output_large" 
#Directory to keep some Results
large_dir = base_dir / "2_Micro/data_visual"
large_dir.mkdir(parents=True, exist_ok=True)'''

if Path("/kaggle").exists():
    print("Running in Kaggle environment")
    # For Kaggle work# Input datasets (read-only in Kaggle) # Files in small input directory
    base_dir = Path("/kaggle/input/new-picrust/")  
    abundance_excel = base_dir / "merged_to_sequence.xlsx" # inside input small sizes input
    #Input market groups
    market_dir = Path("/kaggle/input/a-corrosion-markers")
    output_base = Path("/kaggle/working/")
    output_large = output_base
    #Directory to keep some Results
    large_dir = output_base
else:
    print("Running in local (VSCode) environment")
    base_dir = Path("/home/beatriz/MIC/")
    # Base Paths for local environment
    abundance_excel = base_dir / "2_Micro/data_Ref/merged_to_sequence.xlsx"
    #Input market groups
    market_dir = base_dir / "output_large/A_corrosion_markers.parquet"  # Directory, not file
    output_base = base_dir
    output_large = base_dir / "output_large"
    #Directory to keep some Results
    large_dir = base_dir / "data_visual"

parquet_files = list(market_dir.glob("*.parquet"))
print(f"Found {len(parquet_files)} parquet files")

# Load all files into a dictionary of DataFrames
corrosion_report = {file.stem: pd.read_parquet(file) for file in parquet_files}
#==============================================================================================================

# Integrated taxa from origin genus as headers with levels 6 for the genera, 7 for the GID, muss be cleaned
Integrated_T = pd.read_excel(abundance_excel, sheet_name='core_check_usual_taxa', header=[0,1,2,3,4,5,6,7], engine ='openpyxl')
# Drop first row (index 0) and first column in one chain
Integrated_T = Integrated_T.drop(index=0).drop(Integrated_T.columns[0], axis=1)
Integrated_T= Integrated_T.astype({'Sites': str})
Integrated_T['Sites'] = Integrated_T['Sites'].fillna('Source')
# Remove 'Unnamed' level names
Integrated_T.columns = Integrated_T.columns.map(lambda x: tuple('' if 'Unnamed' in str(level) else level for level in x))
# Changing dtypes to category whiles respecting structure
Integrated_T["Category"] = Integrated_T["Category"].astype("Int64")
Integrated_T= Integrated_T.set_index("Sites")
pre_Integrated = Integrated_T.T
# Define category dict outside
category_dict = Integrated_T.T.iloc[0, 0:-1].to_dict()
# Main analysis dataframes (core pipeline)
pattern_data = corrosion_report["pattern_data"] 
integrated_results = corrosion_report["integrated_results"]
classified_results = corrosion_report["classified_results"]
increasing_markers = corrosion_report["increasing_markers"]
prioritized_markers = corrosion_report["prioritized_markers"]
balanced_markers = corrosion_report["balanced_markers"]
inverse_markers = corrosion_report["inverse_markers"]

# Top marker groups (overall rankings)
top_markers = corrosion_report["group_top_markers"]
significant_markers = corrosion_report["group_significant_markers"]
high_prevalence = corrosion_report["group_high_prevalence"]
very_high_specificity = corrosion_report["group_very_high_specificity"]
# high_frequency = corrosion_report["group_high_frequency"]  # This key doesn't exist - commented out

# Relevance-based groups
high_metals_relevance = corrosion_report["group_high_metals_relevance"]
high_mechanism_relevance = corrosion_report["group_high_mechanism_relevance"]
high_pathway_relevance = corrosion_report["group_high_pathway_relevance"]
high_tier_relevance = corrosion_report["group_high_tier_relevance"]
high_corrosion_relevance = corrosion_report["group_high_corrosion_relevance"]
high_biological_relevance = corrosion_report["group_high_biological_relevance"]

# Consolidated groups
mechanism_all = corrosion_report["group_mechanism_all"]
metals_consolidated = corrosion_report["group_metals_consolidated"]
pathways_all = corrosion_report["group_pathways_all"]
functional_categories = corrosion_report["group_functional_categories"]

# Special groups
metal_iron_sulfur = corrosion_report["group_metal_iron_sulfur"]
high_synergy_markers = corrosion_report["group_high_synergy_markers"]
organic_metal_synergy = corrosion_report["group_organic_metal_synergy"]
corrosion_critical = corrosion_report["group_corrosion_critical"]
#===========================================================================================

# Define the category colors as in original code
category_colors = {
    1: '#008800',  # Dark green - Normal Operation
    2: '#FF8C00',  # Dark orange - Early Warning
    3: '#FF0000'   # Red - System Failure
}

categories_labels = {
    1: 'Normal Operation',
    2: 'Early Warning',
    3: 'System Failure'
}
relevant_categories = { 'hydrogen_metabolism',    'oxygen_metabolism',   'nitrogen_metabolism',    'manganese_processes',  'electron_transfer', 
                        'iron_sulfur_redox', 'ocre_formation',    'sulfur_metabolism',       'organic_acid_metabolism',  'metal_organic_interaction',
                          'biofilm_formation',  'carbon_metabolism', 'ph_modulation',     'temp_response',    'halogen_related',    'methanogenesis'}

relevant_categories = list(relevant_categories)

# 1 ----------------- VISUALIZATION FUNCTIONS -----------------
# 1 HEADMAP =============================================================================
def create_protein_distribution_heatmap(df, top_n=30, min_categories=2):
    """
    Creates a heatmap showing the distribution of top proteins across pathway categories
    
    Args:
        df: DataFrame with protein and pathway data
        top_n: Number of top proteins to display
        min_categories: Minimum number of categories a protein must appear in
    """
    import plotly.graph_objects as go
    import pandas as pd
    import numpy as np
    import corrosion_scoring as cs
    
    # Create unique protein-genus combinations
    df['protein_genus'] = df['protein_name'] + '___' + df['Genus']
    
    # Extract pathway categories from the corrosion scoring system
    unique_categories = sorted(cs.pathway_categories.keys())
    
    # Create a matrix of protein vs category
    protein_category_matrix = {}
    
    # Process each protein-genus pair
    for _, row in df.iterrows():
        protein_genus = row['protein_name'] + '___' + row['Genus']
        
        if protein_genus not in protein_category_matrix:
            protein_category_matrix[protein_genus] = {cat: 0 for cat in unique_categories}
        
        # Check for matches in pathway_categories, functional_categories, and corrosion_mechanisms
        for col in ['functional_categories', 'corrosion_mechanisms', 'pathway_matches']:
            if col in df.columns and isinstance(row.get(col), str):
                col_value = row[col]
                # For each pathway category in cs
                for category, terms in cs.pathway_categories.items():
                    # Check if any term from this category appears in the column
                    if any(term.lower() in col_value.lower() for term in terms):
                        protein_category_matrix[protein_genus][category] = 1
    
    # Convert to DataFrame
    heatmap_df = pd.DataFrame.from_dict(protein_category_matrix, orient='index')
    
    # Calculate number of categories for each protein
    heatmap_df['num_categories'] = heatmap_df.sum(axis=1)
    
    # Filter proteins that appear in at least min_categories
    heatmap_df = heatmap_df[heatmap_df['num_categories'] >= min_categories]
    
    # Get top proteins by number of categories
    top_proteins = heatmap_df.sort_values('num_categories', ascending=False).head(top_n)
    
    # Prepare data for heatmap
    heatmap_data = top_proteins.drop('num_categories', axis=1)
    
    # Create figure
    fig = go.Figure()
    
    # Add heatmap
    fig.add_trace(go.Heatmap(
        z=heatmap_data.values,
        x=heatmap_data.columns,
        y=heatmap_data.index,
        colorscale='YlGnBu',
        showscale=True,
        hovertemplate='%{y} in %{x}: %{z}<extra></extra>'
    ))
    
    # Add scatter plot for number of categories
    if len(heatmap_data) > 0:
        sizes = top_proteins['num_categories'] / top_proteins['num_categories'].max() * 20
    else:
        sizes = 10  # Default size if there's no data
    
    fig.add_trace(go.Scatter(
        x=['# Categories'],
        y=top_proteins.index,  # Use index which contains the protein_genus values
        mode='markers',
        marker=dict(
            color='black',
            size=sizes,
            symbol='circle',
            line=dict(width=0)
        ),
        text=top_proteins['num_categories'],
        hovertemplate='%{y} is in %{text} categories<extra></extra>',
        showlegend=False
    ))
    
    # Update layout
    fig.update_layout(
        title=f'Top {len(top_proteins)} Proteins by Number of Pathway Categories',
        xaxis=dict(
            title='Pathway Categories',
            tickangle=45
        ),
        yaxis=dict(
            title='Protein-Genus Combinations',
            tickfont=dict(size=10)
        ),
        height=800,
        width=1000,
        margin=dict(l=250, r=50, t=100, b=100)
    )
    
    return fig
#==========================
# 2 Protein network
 
def create_protein_network_plot(df, selected_pathway_category=None, min_count=2):
    """
    Creates a network visualization showing proteins shared between genera
    
    Args:
        df: DataFrame with protein-genus data
        selected_pathway_category: Optional filter for pathway category
        min_count: Minimum number of genera a protein must appear in
    """
    import corrosion_scoring as cs
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    
    # Filter by pathway category if specified and it exists in cs.pathway_categories
    if selected_pathway_category and selected_pathway_category in cs.pathway_categories:
        # Get the terms for this pathway category
        pathway_terms = cs.pathway_categories[selected_pathway_category]
        
        # Filter rows where any of these terms appear in functional_categories or corrosion_mechanisms
        filtered_indices = []
        for idx, row in df.iterrows():
            for col in ['functional_categories', 'corrosion_mechanisms', 'pathway_matches']:
                if col in df.columns and isinstance(row.get(col), str):
                    if any(term.lower() in row[col].lower() for term in pathway_terms):
                        filtered_indices.append(idx)
                        break
        
        if filtered_indices:
            df = df.loc[filtered_indices]
    
    # Create protein to genus mapping - keeping only the unique combinations
    protein_to_genus = {}
    genus_to_category = {}
    
    for idx, row in df.iterrows():
        protein = row['protein_name']
        genus = row['Genus']
        category = row['Category'] 
        
        if protein not in protein_to_genus:
            protein_to_genus[protein] = []
        
        # Only add genus if it's not already in the list for this protein
        if genus not in protein_to_genus[protein]:
            protein_to_genus[protein].append(genus)
        
        # Store the highest category for each genus (assuming higher number means higher risk)
        if genus not in genus_to_category or category > genus_to_category[genus]:
            genus_to_category[genus] = category
    
    # Filter proteins that appear in at least min_count genera
    shared_proteins = {p: g for p, g in protein_to_genus.items() if len(g) >= min_count}
    
    # If no proteins meet the criteria, return a message
    if not shared_proteins:
        fig = go.Figure()
        fig.add_annotation(
            text="No proteins found that match the criteria",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False
        )
        return fig
    
    # Create nodes and edges for the network
    nodes = []
    edges = []
    
    # Add protein nodes
    for protein in shared_proteins:
        # Get representative row for this protein
        protein_rows = df[df['protein_name'] == protein]
        protein_info = {}
        
        if len(protein_rows) > 0:
            # Combine multiple rows of the same protein
            if 'corrosion_mechanisms' in protein_rows:
                mechanisms = set()
                for mech_str in protein_rows['corrosion_mechanisms'].dropna():
                    if isinstance(mech_str, str):
                        mechanisms.update([m.strip() for m in mech_str.split(';')])
                protein_info['mechanisms'] = '; '.join(mechanisms)
                
            if 'functional_categories' in protein_rows:
                functions = set()
                for func_str in protein_rows['functional_categories'].dropna():
                    if isinstance(func_str, str):
                        functions.update([f.strip() for f in func_str.split(';')])
                protein_info['functions'] = '; '.join(functions)
                
            if 'combined_score' in protein_rows:
                protein_info['score'] = protein_rows['combined_score'].mean()
        
        nodes.append({
            'id': protein,
            'label': protein,
            'type': 'protein',
            'size': len(shared_proteins[protein]) * 5,  # Size based on number of genera
            'genera_count': len(shared_proteins[protein]),
            'info': protein_info
        })
    
    # Add genus nodes and edges
    all_genera = set()
    for protein, genera in shared_proteins.items():
        for genus in genera:
            if genus not in all_genera:
                all_genera.add(genus)
                category = genus_to_category.get(genus, 1)
                nodes.append({
                    'id': genus,
                    'label': genus,
                    'type': 'genus',
                    'category': category
                })
            
            edges.append({
                'from': protein,
                'to': genus
            })
    
    # Use a force-directed layout to position nodes
    pos = {}
    
    # Position proteins in a circle
    protein_nodes = [n for n in nodes if n['type'] == 'protein']
    genus_nodes = [n for n in nodes if n['type'] == 'genus']
    
    # Position protein nodes in an inner circle
    for i, node in enumerate(protein_nodes):
        angle = 2 * np.pi * i / len(protein_nodes) if len(protein_nodes) > 0 else 0
        radius = 0.4
        pos[node['id']] = (radius * np.cos(angle), radius * np.sin(angle))
    
    # Position genus nodes in an outer circle
    for i, node in enumerate(genus_nodes):
        angle = 2 * np.pi * i / len(genus_nodes) if len(genus_nodes) > 0 else 0
        radius = 0.8
        pos[node['id']] = (radius * np.cos(angle), radius * np.sin(angle))
    
    # Create the figure
    fig = go.Figure()
    
    # Add edges as lines
    for edge in edges:
        x0, y0 = pos[edge['from']]
        x1, y1 = pos[edge['to']]
        
        fig.add_trace(go.Scatter(
            x=[x0, x1],
            y=[y0, y1],
            mode='lines',
            line=dict(width=0.5, color='gray'),
            hoverinfo='none',
            showlegend=False
        ))
    
    # Add protein nodes
    for node in protein_nodes:
        x, y = pos[node['id']]
        
        # Create hover template
        hover_text = f"{node['label']}<br>Found in {node['genera_count']} genera<br>"
        
        if 'info' in node:
            info = node['info']
            if 'mechanisms' in info:
                hover_text += f"Mechanisms: {info['mechanisms']}<br>"
            if 'functions' in info:
                hover_text += f"Functions: {info['functions']}<br>"
            if 'score' in info:
                hover_text += f"Score: {info['score']:.2f}"
        
        hover_text += "<extra></extra>"
        
        fig.add_trace(go.Scatter(
            x=[x],
            y=[y],
            mode='markers+text',
            marker=dict(
                size=node['size'],
                color='skyblue',
                line=dict(width=1, color='black')
            ),
            text=node['label'],
            textposition="top center",
            name=f"{node['label']} (in {node['genera_count']} genera)",
            hovertemplate=hover_text
        ))
    
    # Add genus nodes, colored by category
    category_colors = {1: 'green', 2: 'orange', 3: 'red'}  # Adjust colors as needed
    categories_labels = {1: 'Category 1', 2: 'Category 2', 3: 'Category 3'}
    
    for node in genus_nodes:
        x, y = pos[node['id']]
        category = node['category']
        
        fig.add_trace(go.Scatter(
            x=[x],
            y=[y],
            mode='markers+text',
            marker=dict(
                size=15,
                color=category_colors.get(category, 'gray'),
                symbol='diamond',
                line=dict(width=1, color='black')
            ),
            text=node['label'],
            textposition="bottom center",
            name=f"{node['label']} (Category {category})",
            hovertemplate=(
                f"{node['label']}<br>" +
                f"Category: {category}" +
                f"<extra></extra>"
            )
        ))
    
    # Update layout
    title = "Protein-Genus Network"
    if selected_pathway_category:
        title += f" for {selected_pathway_category.replace('_', ' ').title()}"
    
    fig.update_layout(
        title=title,
        showlegend=True,
        hovermode='closest',
        xaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False
        ),
        yaxis=dict(
            showgrid=False,
            zeroline=False,
            showticklabels=False
        ),
        height=700,
        width=900,
        legend=dict(
            itemsizing='constant',
            orientation='v',
            x=1.05,
            y=0.5
        )
    )
    
    return fig
#===========================================================================================

# 3. Catetory overlap matrix
def create_category_overlap_matrix(df):
    """
    Creates a matrix showing the overlap between different pathway categories
    """
    # Extract all unique pathway categories
    all_categories = []
    
    for idx, row in df.iterrows():
        categories = str(row.get('pathway_categories', '')).split(',')
        categories = [c.strip() for c in categories if c.strip()]
        all_categories.extend(categories)
    
    all_categories = sorted(list(set(all_categories)))
    
    # Create a matrix to store overlaps
    overlap_matrix = np.zeros((len(all_categories), len(all_categories)))
    
    # For each protein-genus pair, update the matrix
    for idx, row in df.iterrows():
        protein_genus = f"{row['Genus']} - {row['protein_name']}"
        categories = str(row.get('pathway_categories', '')).split(',')
        categories = [c.strip() for c in categories if c.strip()]
        
        # Update the matrix for each pair of categories
        for i, cat1 in enumerate(categories):
            if cat1 not in all_categories:
                continue
                
            idx1 = all_categories.index(cat1)
            
            # Self-overlap
            overlap_matrix[idx1, idx1] += 1
            
            # Cross-overlap
            for j, cat2 in enumerate(categories[i+1:], i+1):
                if cat2 not in all_categories:
                    continue
                    
                idx2 = all_categories.index(cat2)
                overlap_matrix[idx1, idx2] += 1
                overlap_matrix[idx2, idx1] += 1
    
    # Create the heatmap
    fig = go.Figure(data=go.Heatmap(
        z=overlap_matrix,
        x=all_categories,
        y=all_categories,
        colorscale='Viridis',
        text=overlap_matrix.astype(int),
        texttemplate="%{text}",
        hovertemplate='%{y} and %{x} share %{z:.0f} protein-genus pairs<extra></extra>'
    ))
    
    # Update layout
    fig.update_layout(
        title='Pathway Category Overlap Matrix',
        xaxis=dict(
            title='Pathway Categories',
            tickangle=45
        ),
        yaxis=dict(
            title='Pathway Categories'
        ),
        height=700,
        width=900
    )
    
    return fig
# 4. category comparison
def plot_category_comparison_advanced(df, pathway_cat1, pathway_cat2, abundance_metric='log2fc_3vs1',
                                    max_items=20, highlight_shared=True, plot_height=800, plot_width=1200):
    """
    Create an improved mirror plot comparing pathway categories with better highlighting of shared items.
    """
    # Create a single figure with two subplots side by side
    
    # Build lookup for categories
    fig = make_subplots(rows=1, cols=2, 
                    subplot_titles=[pathway_cat1.replace('_', ' ').title(), 
                                    pathway_cat2.replace('_', ' ').title()],
                    horizontal_spacing=0.1)
    
    category_terms = {k: v for k, v in cs.pathway_categories.items()}

    # Define matching function
    def match_pathways(row, terms):
        return any(term.lower() in str(row.get('pathways', '')).lower() for term in terms)

    # Identify items in each pathway category
    items_cat1 = {}
    items_cat2 = {}
    all_proteins = set()

    for idx, row in df.iterrows():
        protein = row['protein_name']
        genus = row['Genus']
        item_id = f"{genus} - {protein}"
        all_proteins.add(protein)
        
        # Check pathway categories using the notebook logic
        if match_pathways(row, category_terms.get(pathway_cat1, [])):
            items_cat1[item_id] = row
        
        if match_pathways(row, category_terms.get(pathway_cat2, [])):
            items_cat2[item_id] = row

    #====================================0000         
    # Find items that are shared
    shared_items = set(items_cat1.keys()).intersection(set(items_cat2.keys()))
    
    # Find proteins that appear in both but in different genera
    proteins_in_both = {}
    for item_id in list(items_cat1.keys()) + list(items_cat2.keys()):
        if ' - ' in item_id:
            genus, protein = item_id.split(' - ', 1)
            if protein not in proteins_in_both:
                proteins_in_both[protein] = {'cat1': set(), 'cat2': set()}
            
            if item_id in items_cat1:
                proteins_in_both[protein]['cat1'].add(genus)
            
            if item_id in items_cat2:
                proteins_in_both[protein]['cat2'].add(genus)
    
    # Identify proteins in both categories but in different genera
    proteins_across_categories = {p: data for p, data in proteins_in_both.items() 
                                if data['cat1'] and data['cat2'] and len(data['cat1'].union(data['cat2'])) > 1}
    
    # Sort and limit items
    def sort_items(items_dict):
        return sorted(
            items_dict.items(),
            key=lambda x: (-x[1].get('Category', 1), -abs(x[1].get(abundance_metric, 0)))
        )[:max_items]
    
    sorted_items_cat1 = sort_items(items_cat1)
    sorted_items_cat2 = sort_items(items_cat2)
    
    # Plot left side
    for item_id, row in sorted_items_cat1:
        category = row.get('Category', 1)
        color = category_colors.get(category, 'gray')
        
        is_shared = item_id in shared_items
        
        # Check if protein is in both categories but different genera
        is_cross_category = False
        protein = item_id.split(' - ', 1)[1] if ' - ' in item_id else item_id
        if protein in proteins_across_categories and not is_shared:
            is_cross_category = True
        
        marker_size = 15 if category == 3 else 12 if is_shared or is_cross_category else 10
        marker_symbol = 'star' if category == 3 else 'circle'
        marker_line = dict(width=2, color='black') if is_shared else \
                     dict(width=1.5, color='purple') if is_cross_category else \
                     dict(width=1, color='white')
        
        x_val = row.get(abundance_metric, 0)
        
        fig.add_trace(
            go.Scatter(
                x=[x_val],
                y=[item_id],
                mode='markers',
                marker=dict(
                    color=color,
                    size=marker_size,
                    symbol=marker_symbol,
                    line=marker_line
                ),
                hovertemplate=(
                    f"{item_id}<br>"
                    f"Category: {categories_labels.get(category, 'Unknown')}<br>"
                    f"{abundance_metric}: {x_val:.2f}"
                    f"{'<br>Present in both categories' if is_shared else ''}"
                    f"{'<br>Protein in both categories across genera' if is_cross_category else ''}"
                    f"<extra></extra>"
                ),
                showlegend=False
            ),
            row=1, col=1
        )
    
    # Plot right side
    for item_id, row in sorted_items_cat2:
        category = row.get('Category', 1)
        color = category_colors.get(category, 'gray')
        
        is_shared = item_id in shared_items
        
        # Check if protein is in both categories but different genera
        is_cross_category = False
        protein = item_id.split(' - ', 1)[1] if ' - ' in item_id else item_id
        if protein in proteins_across_categories and not is_shared:
            is_cross_category = True
        
        marker_size = 15 if category == 3 else 12 if is_shared or is_cross_category else 10
        marker_symbol = 'star' if category == 3 else 'circle'
        marker_line = dict(width=2, color='black') if is_shared else \
                     dict(width=1.5, color='purple') if is_cross_category else \
                     dict(width=1, color='white')
        
        x_val = row.get(abundance_metric, 0)
        
        fig.add_trace(
            go.Scatter(
                x=[x_val],
                y=[item_id],
                mode='markers',
                marker=dict(
                    color=color,
                    size=marker_size,
                    symbol=marker_symbol,
                    line=marker_line
                ),
                hovertemplate=(
                    f"{item_id}<br>"
                    f"Category: {categories_labels.get(category, 'Unknown')}<br>"
                    f"{abundance_metric}: {x_val:.2f}"
                    f"{'<br>Present in both categories' if is_shared else ''}"
                    f"{'<br>Protein in both categories across genera' if is_cross_category else ''}"
                    f"<extra></extra>"
                ),
                showlegend=False
            ),
            row=1, col=2
        )
    
    # Add legend for the categories
    for cat, label in categories_labels.items():
        symbol = 'star' if cat == 3 else 'circle'
        size = 15 if cat == 3 else 10
        
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode='markers',
                marker=dict(
                    color=category_colors.get(cat, 'gray'),
                    size=size,
                    symbol=symbol
                ),
                name=label,
                showlegend=True
            ),
            row=1, col=1
        )
    
    # Add legend for shared items
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(
                color='gray',
                size=12,
                symbol='circle',
                line=dict(width=2, color='black')
            ),
            name="Same genus-protein in both categories",
            showlegend=True
        ),
        row=1, col=1
    )
    
    # Add legend for cross-category proteins
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(
                color='gray',
                size=12,
                symbol='circle',
                line=dict(width=1.5, color='purple')
            ),
            name="Same protein in both categories (different genera)",
            showlegend=True
        ),
        row=1, col=1
    )
    
    # Calculate appropriate x-axis range
    all_x_vals = []
    for item_dict in [items_cat1, items_cat2]:
        for row in item_dict.values():
            all_x_vals.append(row.get(abundance_metric, 0))
    
    if all_x_vals:
        x_min = min(all_x_vals)
        x_max = max(all_x_vals)
        x_padding = (x_max - x_min) * 0.1
        x_range = [x_min - x_padding, x_max + x_padding]
    else:
        x_range = [-1, 1]
    
    # Update layout
    fig.update_layout(
        title=f"Comparison of {pathway_cat1.replace('_', ' ').title()} vs {pathway_cat2.replace('_', ' ').title()}",
        height=plot_height,
        width=plot_width,
        margin=dict(l=200, r=50, t=100, b=50),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    # Update axes
    fig.update_xaxes(title_text=abundance_metric, range=x_range, row=1, col=1)
    fig.update_xaxes(title_text=abundance_metric, range=x_range, row=1, col=2)
    
    # Add a note about the number of items
    fig.add_annotation(
        text=f"{len(items_cat1)} items total ({len(shared_items)} shared)",
        xref="x domain",
        yref="paper",
        x=0.5,
        y=-0.15,
        showarrow=False,
        row=1, col=1
    )
    
    fig.add_annotation(
        text=f"{len(items_cat2)} items total ({len(shared_items)} shared)",
        xref="x2 domain",
        yref="paper",
        x=0.5,
        y=-0.15,
        showarrow=False,
        row=1, col=2
    )
    
    return fig

# 5 ----------------- MAIN DASHBOARD UI -----------------

def create_dashboard(df, relevant_categories):
    """
    Creates a dashboard with multiple visualization options for exploring protein-genus relationships
    """
    # Create tab structure
    tab_titles = ['Category Comparison', 'Protein Distribution', 'Network View', 'Category Overlap']
    tab_contents = [widgets.VBox([]) for _ in range(len(tab_titles))]
    tabs = widgets.Tab(children=tab_contents)
    for i, title in enumerate(tab_titles):
        tabs.set_title(i, title)
    
    # 6. --------------- TAB 1: CATEGORY COMPARISON ---------------
    
    # Create the widgets for category comparison
    cat1_dropdown = widgets.Dropdown(
        options=relevant_categories,
        value=relevant_categories[0],
        description='Left Category:',
        style={'description_width': 'initial'}
    )
    
    cat2_dropdown = widgets.Dropdown(
        options=relevant_categories,
        value=relevant_categories[1] if len(relevant_categories) > 1 else relevant_categories[0],
        description='Right Category:',
        style={'description_width': 'initial'}
    )
    
    abundance_dropdown = widgets.Dropdown(
        options=[
            ('Mean Cat 1', 'mean_cat1'), 
            ('Mean Cat 2', 'mean_cat2'), 
            ('Mean Cat 3', 'mean_cat3'),
            ('Fold Change 2vs1', 'fold_change_2vs1'), 
            ('Log2 FC 2vs1', 'log2fc_2vs1'),
            ('Fold Change 3vs2', 'fold_change_3vs2'), 
            ('Log2 FC 3vs2', 'log2fc_3vs2'),
            ('Fold Change 3vs1', 'fold_change_3vs1'), 
            ('Log2 FC 3vs1', 'log2fc_3vs1')
        ],
        value='log2fc_3vs1',
        description='Metric:',
        style={'description_width': 'initial'}
    )
    
    max_items_slider = widgets.IntSlider(
        value=20,
        min=5,
        max=50,
        step=5,
        description='Max items:',
        style={'description_width': 'initial'}
    )
    
    plot_height_slider = widgets.IntSlider(
        value=800,
        min=400,
        max=1500,
        step=100,
        description='Height:',
        style={'description_width': 'initial'}
    )
    
    plot_width_slider = widgets.IntSlider(
        value=1200,
        min=800,
        max=1800,
        step=100,
        description='Width:',
        style={'description_width': 'initial'}
    )
    
    highlight_shared_checkbox = widgets.Checkbox(
        value=True,
        description='Highlight shared items',
        style={'description_width': 'initial'}
    )
    
    # Output for category comparison tab
    output_comparison = widgets.Output()
    
    # Function to update the category comparison plot
    def update_comparison_plot(change):
        output_comparison.clear_output(wait=True)
        with output_comparison:
            fig = plot_category_comparison_advanced(
                df,
                cat1_dropdown.value,
                cat2_dropdown.value,
                abundance_dropdown.value,
                max_items_slider.value,
                highlight_shared_checkbox.value,
                plot_height_slider.value,
                plot_width_slider.value
            )
            fig.show()
    
    #  7. Register callbacks
    cat1_dropdown.observe(update_comparison_plot, names='value')
    cat2_dropdown.observe(update_comparison_plot, names='value')
    abundance_dropdown.observe(update_comparison_plot, names='value')
    max_items_slider.observe(update_comparison_plot, names='value')
    plot_height_slider.observe(update_comparison_plot, names='value')
    plot_width_slider.observe(update_comparison_plot, names='value')
    highlight_shared_checkbox.observe(update_comparison_plot, names='value')
    
    # Update button for explicit refresh
    update_comparison_button = widgets.Button(
        description='Update Plot',
        button_style='primary',
        tooltip='Click to refresh the plot',
        icon='refresh'
    )
    update_comparison_button.on_click(lambda b: update_comparison_plot(None))
    
    # Organize widgets for category comparison tab
    category_controls = widgets.HBox([cat1_dropdown, cat2_dropdown, abundance_dropdown])
    display_controls = widgets.HBox([max_items_slider, plot_height_slider, plot_width_slider])
    option_controls = widgets.HBox([highlight_shared_checkbox, update_comparison_button])
    
    # Assemble tab content
    tab_contents[0] = widgets.VBox([
        widgets.HTML("<h3>Pathway Category Comparison</h3>"),
        widgets.HTML("<p>Compare proteins across two pathway categories to identify shared elements.</p>"),
        category_controls,
        display_controls,
        option_controls,
        output_comparison
    ])
    tabs.children = tab_contents
    
    # Initial plot for comparison tab
    update_comparison_plot(None)
    
    # 8 --------------- TAB 2: PROTEIN DISTRIBUTION ---------------
    
    # Create widgets for protein distribution
    top_n_slider = widgets.IntSlider(
        value=30,
        min=10,
        max=100,
        step=5,
        description='Top proteins:',
        style={'description_width': 'initial'}
    )
    
    min_categories_slider = widgets.IntSlider(
        value=2,
        min=1,
        max=8,
        step=1,
        description='Min categories:',
        style={'description_width': 'initial'}
    )
    
    # Output for protein distribution tab
    output_distribution = widgets.Output()
    
    # Function to update the protein distribution plot
    def update_distribution_plot(change):
        output_distribution.clear_output(wait=True)
        with output_distribution:
            fig = create_protein_distribution_heatmap(
                df,
                top_n=top_n_slider.value,
                min_categories=min_categories_slider.value
            )
            fig.show()
    
    # Register callbacks
    top_n_slider.observe(update_distribution_plot, names='value')
    min_categories_slider.observe(update_distribution_plot, names='value')
    
    # Update button for explicit refresh
    update_distribution_button = widgets.Button(
        description='Update Plot',
        button_style='primary',
        tooltip='Click to refresh the plot',
        icon='refresh'
    )
    update_distribution_button.on_click(lambda b: update_distribution_plot(None))
    
    # Organize widgets for protein distribution tab
    distribution_controls = widgets.HBox([top_n_slider, min_categories_slider, update_distribution_button])
    
    # Assemble tab content
    tab_contents[0] = widgets.VBox([
        widgets.HTML("<h3>Protein Distribution Across Pathway Categories</h3>"),
        widgets.HTML("<p>Visualize which proteins appear in multiple pathway categories.</p>"),
        distribution_controls,
        output_distribution
    ])
    tabs.children = tab_contents

    # Initial plot for distribution tab
    update_distribution_plot(None)
    
    # 9 --------------- TAB 3: NETWORK VIEW ---------------
    
    # Create widgets for network view
    network_category_dropdown = widgets.Dropdown(
        options=[('All Categories', None)] + [(cat, cat) for cat in relevant_categories],
        value=None,
        description='Filter by category:',
        style={'description_width': 'initial'}
    )
    
    min_genera_slider = widgets.IntSlider(
        value=2,
        min=1,
        max=10,
        step=1,
        description='Min genera:',
        style={'description_width': 'initial'}
    )
    
    # Output for network view tab
    output_network = widgets.Output()
    
    # Function to update the network view
    def update_network_plot(change):
        output_network.clear_output(wait=True)
        with output_network:
            fig = create_protein_network_plot(
                df,
                selected_pathway_category=network_category_dropdown.value,
                min_count=min_genera_slider.value
            )
            fig.show()
    
    # Register callbacks
    network_category_dropdown.observe(update_network_plot, names='value')
    min_genera_slider.observe(update_network_plot, names='value')
    
    # Update button for explicit refresh
    update_network_button = widgets.Button(
        description='Update Plot',
        button_style='primary',
        tooltip='Click to refresh the plot',
        icon='refresh'
    )
    update_network_button.on_click(lambda b: update_network_plot(None))
    
    # Organize widgets for network view tab
    network_controls = widgets.HBox([network_category_dropdown, min_genera_slider, update_network_button])
    
    # Assemble tab content
    tab_contents[0] = widgets.VBox([
        widgets.HTML("<h3>Protein-Genus Network</h3>"),
        widgets.HTML("<p>Explore network connections between proteins and genera.</p>"),
        network_controls,
        output_network
    ])
    tabs.children = tab_contents

    # Initial plot for network tab
    update_network_plot(None)
    
    # 10 --------------- TAB 4: CATEGORY OVERLAP ---------------
    
    # Output for category overlap tab
    output_overlap = widgets.Output()
    
    # Function to update the category overlap matrix
    def update_overlap_matrix(change=None):
        output_overlap.clear_output(wait=True)
        with output_overlap:
            fig = create_category_overlap_matrix(df)
            fig.show()
    
    # Update button for explicit refresh
    update_overlap_button = widgets.Button(
        description='Update Matrix',
        button_style='primary',
        tooltip='Click to refresh the matrix',
        icon='refresh'
    )
    update_overlap_button.on_click(update_overlap_matrix)
    
    # Assemble tab content
    tab_contents[0] = widgets.VBox([
        widgets.HTML("<h3>Pathway Category Overlap Matrix</h3>"),
        widgets.HTML("<p>View how different pathway categories share protein-genus pairs.</p>"),
        update_overlap_button,
        output_overlap
    ])
    tabs.children = tab_contents

    
    # Initial plot for overlap tab
    update_overlap_matrix()
    
    # Display the dashboard
    display(tabs)
    
    return tabs

# 11 ----------------- UTILITY FUNCTIONS -----------------

def analyze_protein_genus_distribution(df):
    """
    Analyzes and summarizes the distribution of proteins across genera and pathway categories
    """
    # Count proteins per genus
    genus_counts = df['Genus'].value_counts()
    
    # Count proteins per pathway category
    category_counts = {}
    for idx, row in df.iterrows():
        categories = str(row.get('pathway_categories', '')).split(',')
        categories = [c.strip() for c in categories if c.strip()]
        
        for category in categories:
            if category not in category_counts:
                category_counts[category] = 0
            category_counts[category] += 1
    
    # Find proteins that appear in multiple genera
    protein_to_genera = {}
    for idx, row in df.iterrows():
        protein = row['protein_name']
        genus = row['Genus']
        
        if protein not in protein_to_genera:
            protein_to_genera[protein] = set()
        
        protein_to_genera[protein].add(genus)
    
    # Sort proteins by number of genera they appear in
    multi_genus_proteins = {
        protein: genera 
        for protein, genera in protein_to_genera.items() 
        if len(genera) > 1
    }
    
    sorted_multi_genus = sorted(
        multi_genus_proteins.items(),
        key=lambda x: len(x[1]),
        reverse=True
    )
    
    # Find proteins that appear in multiple pathway categories
    protein_to_categories = {}
    for idx, row in df.iterrows():
        protein = row['protein_name']
        categories = str(row.get('pathway_categories', '')).split(',')
        categories = [c.strip() for c in categories if c.strip()]
        
        if protein not in protein_to_categories:
            protein_to_categories[protein] = set()
        
        protein_to_categories[protein].update(categories)
    
    # Sort proteins by number of categories they appear in
    multi_category_proteins = {
        protein: categories 
        for protein, categories in protein_to_categories.items() 
        if len(categories) > 1
    }
    
    sorted_multi_category = sorted(
        multi_category_proteins.items(),
        key=lambda x: len(x[1]),
        reverse=True
    )
    
    # Create summary dataframes
    genus_summary = pd.DataFrame({
        'Genus': genus_counts.index,
        'Protein Count': genus_counts.values
    })
    
    category_summary = pd.DataFrame({
        'Category': list(category_counts.keys()),
        'Protein Count': list(category_counts.values())
    }).sort_values('Protein Count', ascending=False)
    
    multi_genus_summary = pd.DataFrame([
        {'Protein': protein, 'Genera Count': len(genera), 'Genera': ', '.join(sorted(genera))}
        for protein, genera in sorted_multi_genus[:20]  # Show top 20
    ])
    
    multi_category_summary = pd.DataFrame([
        {'Protein': protein, 'Category Count': len(categories), 'Categories': ', '.join(sorted(categories))}
        for protein, categories in sorted_multi_category[:20]  # Show top 20
    ])
    
    return {
        'genus_summary': genus_summary,
        'category_summary': category_summary,
        'multi_genus_summary': multi_genus_summary,
        'multi_category_summary': multi_category_summary
    }
#12
def export_visualizations(df, relevant_categories, output_dir ='protein_visualizations'):
    """
    Exports all visualizations to static HTML files for sharing
    """
    import os
    from IPython.display import HTML, display
    output_dir = output_share
    output_dir_str = str(output_dir)
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir_str):
        os.makedirs(output_dir_str)
    
    # Create index.html
    index_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Protein-Genus Visualization Dashboard</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1 {{ color: #333366; }}
            h2 {{ color: #666699; margin-top: 30px; }}
            ul {{ list-style-type: none; padding-left: 0; }}
            li {{ margin: 10px 0; }}
            a {{ color: #3366cc; text-decoration: none; }}
            a:hover {{ text-decoration: underline; }}
        </style>
    </head>
    <body>
        <h1>Protein-Genus Visualization Dashboard</h1>
        
        <h2>Pathway Category Comparisons</h2>
        <ul>
    """
    
    # 13. Generate category comparison plots
    metric = 'log2fc_3vs1'  # Default metric
    max_items = 20
    
    for i, cat1 in enumerate(relevant_categories):
        for j, cat2 in enumerate(relevant_categories[i+1:], i+1):
            filename = f"comparison_{cat1}_vs_{cat2}.html"
            filepath = os.path.join(output_dir_str, filename)
            
            fig = plot_category_comparison_advanced(
                df, cat1, cat2, metric, max_items, True, 800, 1200
            )
            fig.write_html(filepath)
            
            index_content += f'            <li><a href="{filename}" target="_blank">{cat1.replace("_", " ").title()} vs {cat2.replace("_", " ").title()}</a></li>\n'
    
    # Add protein distribution heatmap
    index_content += """
        </ul>
        
        <h2>Protein Distribution</h2>
        <ul>
    """
    
    filename = "protein_distribution.html"
    filepath = os.path.join(output_dir_str, filename) 
    fig = create_protein_distribution_heatmap(df, top_n=30, min_categories=2)
    fig.write_html(filepath)
    index_content += f'            <li><a href="{filename}" target="_blank">Protein Distribution Across Categories</a></li>\n'
    
    # Add network visualizations
    index_content += """
        </ul>
        
        <h2>Network Visualizations</h2>
        <ul>
    """

    # Overall network
    filename = "network_all.html"
    filepath = os.path.join(output_dir_str, filename)
    fig = create_protein_network_plot(df, selected_pathway_category=None, min_count=2)
    fig.write_html(filepath)
    index_content += f'            <li><a href="{filename}" target="_blank">All Categories Network</a></li>\n'
    
    # Individual category networks
    for cat in relevant_categories:
        filename = f"network_{cat}.html"
        filepath = os.path.join(output_dir_str, filename)
        fig = create_protein_network_plot(df, selected_pathway_category=cat, min_count=2)
        fig.write_html(filepath)
        index_content += f'            <li><a href="{filename}" target="_blank">{cat.replace("_", " ").title()} Network</a></li>\n'
    
    # Add category overlap matrix
    index_content += """
        </ul>
        
        <h2>Category Overlap Matrix</h2>
        <ul>
    """

    filename = "category_overlap.html"
    filepath = os.path.join(output_dir_str, filename)
    fig = create_category_overlap_matrix(df)
    fig.write_html(filepath)
    index_content += f'            <li><a href="{filename}" target="_blank">Pathway Category Overlap Matrix</a></li>\n'
    
    # Close HTML
    index_content += """
        </ul>
    </body>
    </html>
    """
    
    # Write index.html
    with open(os.path.join(output_dir_str, 'index.html'), 'w') as f:
        f.write(index_content)
    
    return HTML(f'<p>Visualizations exported to <code>{output_dir_str}</code> directory.</p>')

# 14 ----------------- CALLING THE DASHBOARD -----------------

output_share = Path("/home/beatriz/SharedFolder/protein_visualisations")
output_share.mkdir(parents=True, exist_ok=True)
# To use the dashboard with your data:
dashboard = create_dashboard(corrosion_critical, relevant_categories)

# To export static visualizations:
export_visualizations(corrosion_critical, relevant_categories, 'protein_visualizations')