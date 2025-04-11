#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
from pathlib import Path
import os
from pyarrow import parquet

output_large = Path(os.path.expanduser('~/MIC/output_large/'))
parquet_dir = output_large / "corrosion_markers_parquet"

# List all Parquet files in the directory
parquet_files = list(parquet_dir.glob("*.parquet"))

# Load all files into a dictionary of DataFrames
corrosion_report = {file.stem: pd.read_parquet(file) for file in parquet_files}

# Access a specific DataFrame
top_markers = corrosion_report["group_top_markers"]

# Create a copy of the DataFrame (deep=False avoids deep copying for performance)
df = top_markers.copy(deep=False)

# 1. Bar Chart comparing genera with abundance values
def create_genera_abundance_comparison(df, abundance_column='norm_abund_contri'):
    """
    Create a bar chart comparing genera with selected abundance values
    Color by category or mean categories
    """
    plt.figure(figsize=(12, 10))
    
    # First plot: Color by Category
    plt.subplot(2, 1, 1)
    ax1 = sns.barplot(x='Genus', y=abundance_column, data=df, hue='Category', palette='viridis')
    plt.title(f'Genera Abundance Comparison (Colored by Category)')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel(abundance_column)
    plt.legend(title='Category', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Second plot: Color by highest mean category
    plt.subplot(2, 1, 2)
    
    # Determine highest mean category for each row
    df['highest_mean'] = df[['mean_cat1', 'mean_cat2', 'mean_cat3']].idxmax(axis=1)
    df['highest_mean'] = df['highest_mean'].map({'mean_cat1': 'Category 1', 'mean_cat2': 'Category 2', 'mean_cat3': 'Category 3'})
    
    ax2 = sns.barplot(x='Genus', y=abundance_column, data=df, hue='highest_mean', palette='Set2')
    plt.title(f'Genera Abundance Comparison (Colored by Highest Mean Category)')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel(abundance_column)
    plt.legend(title='Highest Mean Category', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.savefig('genera_abundance_comparison.png', dpi=300, bbox_inches='tight')
    
    # Also create an interactive version with Plotly
    fig = make_subplots(rows=2, cols=1, 
                        subplot_titles=(f'Genera Abundance (by Category)', 
                                        f'Genera Abundance (by Highest Mean Category)'))
    
    # First plot with Plotly
    for category in df['Category'].unique():
        df_cat = df[df['Category'] == category]
        fig.add_trace(
            go.Bar(x=df_cat['Genus'], y=df_cat[abundance_column], name=category),
            row=1, col=1
        )
    
    # Second plot with Plotly
    for mean_cat in df['highest_mean'].unique():
        df_mean = df[df['highest_mean'] == mean_cat]
        fig.add_trace(
            go.Bar(x=df_mean['Genus'], y=df_mean[abundance_column], name=mean_cat),
            row=2, col=1
        )
    
    fig.update_layout(height=800, width=1000, title_text="Genera Abundance Comparison")
    fig.update_xaxes(tickangle=45)
    
    # Save as HTML file
    pio.write_html(fig, 'genera_abundance_comparison_interactive.html')
    
    return 'Generated static and interactive abundance comparison plots'

# 2. Sunburst Chart for hierarchical representation of pathways
def create_pathway_sunburst(df, color_by='Category'):
    """
    Create a sunburst chart showing the hierarchical representation of pathways
    """
    # Process the data for the sunburst chart
    # Extract hierarchy paths and prepare for visualization
    hierarchy_data = []
    
    for _, row in df.iterrows():
        hierarchy = row['hierarchy'].split(' > ')
        pathway = row['pathways'].split(' / ')
        
        # Create path for sunburst
        full_path = hierarchy.copy()
        
        # Add the base values for all levels
        for i in range(len(full_path)):
            sub_path = full_path[:i+1]
            hierarchy_data.append({
                'path': ' > '.join(sub_path),
                'value': 1,
                'Category': row['Category'],
                'corrosion_relevance': row['corrosion_relevance']
            })
    
    # Convert to DataFrame
    hierarchy_df = pd.DataFrame(hierarchy_data)
    
    # Create sunburst chart using Plotly
    if color_by == 'Category':
        fig = px.sunburst(
            hierarchy_df, 
            path=['path'],
            values='value',
            color='Category',
            title='Pathway Hierarchy Sunburst (Colored by Category)',
            width=800, 
            height=800
        )
    else:  # color_by == 'corrosion_relevance'
        fig = px.sunburst(
            hierarchy_df, 
            path=['path'],
            values='value',
            color='corrosion_relevance',
            title='Pathway Hierarchy Sunburst (Colored by Corrosion Relevance)',
            color_discrete_map={
                'Very High': '#ff0000',
                'High': '#ff8c00',
                'Medium': '#ffff00',
                'Low': '#008000'
            },
            width=800, 
            height=800
        )
    
    # Save as HTML file
    pio.write_html(fig, f'pathway_sunburst_{color_by}.html')
    
    # Also create a static version with matplotlib (simplified)
    plt.figure(figsize=(10, 10))
    
    # Simplified static plot - just showing the top level categories
    top_level_counts = hierarchy_df['path'].apply(lambda x: x.split(' > ')[0]).value_counts()
    plt.pie(top_level_counts, labels=top_level_counts.index, autopct='%1.1f%%')
    plt.title('Top-Level Pathway Distribution')
    plt.savefig('pathway_distribution_static.png', dpi=300, bbox_inches='tight')
    
    return f'Generated pathway sunburst chart colored by {color_by}'

# Run the functions to create the visualizations
print(create_genera_abundance_comparison(df, 'norm_abund_contri'))
print(create_pathway_sunburst(df, 'Category'))
print(create_pathway_sunburst(df, 'corrosion_relevance'))

# Create a custom sunburst chart for a more precise hierarchical representation
def create_detailed_sunburst(df):
    """
    Create a more detailed sunburst chart that properly represents the pathway hierarchy
    """
    # Process hierarchies to create proper parent-child relationships
    all_paths = []
    
    for _, row in df.iterrows():
        hierarchy_levels = row['hierarchy'].split(' > ')
        pathway_levels = row['pathways'].split(' / ')
        
        current_path = []
        # Add the hierarchy levels
        for level in hierarchy_levels:
            current_path.append(level)
            all_paths.append({
                'id': ' > '.join(current_path),
                'parent': ' > '.join(current_path[:-1]) if len(current_path) > 1 else '',
                'label': level,
                'value': 1,
                'Category': row['Category'],
                'corrosion_relevance': row['corrosion_relevance']
            })
        
        # Add the last pathway level connected to the last hierarchy level
        if pathway_levels:
            final_pathway = pathway_levels[-1]
            all_paths.append({
                'id': ' > '.join(current_path) + ' > ' + final_pathway,
                'parent': ' > '.join(current_path),
                'label': final_pathway,
                'value': 1,
                'Category': row['Category'],
                'corrosion_relevance': row['corrosion_relevance']
            })
    
    # Convert to DataFrame and aggregate duplicate paths
    sunburst_df = pd.DataFrame(all_paths)
    aggregated_df = sunburst_df.groupby(['id', 'parent', 'label', 'Category', 'corrosion_relevance']).sum().reset_index()
    
    # Create the sunburst chart
    fig = go.Figure()
    
    fig.add_trace(go.Sunburst(
        ids=aggregated_df['id'],
        labels=aggregated_df['label'],
        parents=aggregated_df['parent'],
        values=aggregated_df['value'],
        branchvalues='total',
        marker=dict(
            colors=[
                {'Very High': '#ff0000', 'High': '#ff8c00', 'Medium': '#ffff00', 'Low': '#008000'}[rel]
                for rel in aggregated_df['corrosion_relevance']
            ]
        ),
        hovertemplate='<b>%{label}</b><br>Path: %{id}<br>Value: %{value}<br>',
        name=''
    ))
    
    fig.update_layout(
        title_text='Detailed Pathway Hierarchy (Colored by Corrosion Relevance)',
        width=900,
        height=900
    )
    
    pio.write_html(fig, 'detailed_pathway_sunburst.html')
    
    return 'Generated detailed pathway sunburst chart'

print(create_detailed_sunburst(df))

# BONUS: Create a combined visualization that shows both approaches side by side
def create_combined_visualization(df):
    """
    Create a combined visualization showing multiple aspects of the data
    """
    fig = make_subplots(
        rows=2, cols=2,
        specs=[
            [{"type": "bar"}, {"type": "bar"}],
            [{"type": "pie", "colspan": 2}, None]
        ],
        subplot_titles=(
            "Abundance by Category", 
            "Abundance by Mean Categories",
            "Corrosion Relevance Distribution"
        )
    )
    
    # First subplot: Abundance by Category
    for cat in df['Category'].unique():
        df_cat = df[df['Category'] == cat]
        fig.add_trace(
            go.Bar(
                x=df_cat['Genus'],
                y=df_cat['norm_abund_contri'],
                name=cat
            ),
            row=1, col=1
        )
    
    # Second subplot: Abundance by highest mean category
    df['highest_mean'] = df[['mean_cat1', 'mean_cat2', 'mean_cat3']].idxmax(axis=1)
    df['highest_mean'] = df['highest_mean'].map({'mean_cat1': 'Category 1', 'mean_cat2': 'Category 2', 'mean_cat3': 'Category 3'})
    
    for mean_cat in df['highest_mean'].unique():
        df_mean = df[df['highest_mean'] == mean_cat]
        fig.add_trace(
            go.Bar(
                x=df_mean['Genus'],
                y=df_mean['norm_abund_contri'],
                name=mean_cat
            ),
            row=1, col=2
        )
    
    # Third subplot: Corrosion relevance distribution
    relevance_counts = df['corrosion_relevance'].value_counts()
    colors = {'Very High': '#ff0000', 'High': '#ff8c00', 'Medium': '#ffff00', 'Low': '#008000'}
    fig.add_trace(
        go.Pie(
            labels=relevance_counts.index,
            values=relevance_counts.values,
            marker_colors=[colors[rel] for rel in relevance_counts.index]
        ),
        row=2, col=1
    )
    
    fig.update_layout(
        height=1000,
        width=1200,
        title_text="Combined Data Visualization"
    )
    
    fig.update_xaxes(tickangle=45)
    
    pio.write_html(fig, 'combined_visualization.html')
    
    return "Generated combined visualization"

"""
INSTRUCTIONS:

1. To run this code:
   - Save it as a Python file (e.g., visualizations.py)
   - Make sure you have the required libraries installed:
     pip install pandas numpy matplotlib seaborn plotly

2. If you want to use your actual data:
   - Replace the sample DataFrame creation with:
     df = pd.read_csv('your_data_file.csv')
   
3. The script generates several output files:
   - Static images (.png files)
   - Interactive HTML files for the Plotly visualizations

4. For the sunburst chart:
   - The interactive HTML version provides the best experience
   - You can hover over sections to see details
   - You can click on sections to zoom in
"""