# Import necessary libraries
import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from pandas.errors import EmptyDataError
import os
import numpy as np

st.image('SlytheRINs-logo.svg', use_column_width=True)

# Title of the app
st.title('SlytheRINs')

# Description
st.write('SlytheRINs is a Streamlit app to analyze and compare Residue Interaction Networks (RINs) calculated from different protein structures and conformations.')
st.write('Developed by the [EvoMol-Lab](https://github.com/evomol-lab) team at the Bioinformatics Multidisciplinary Environment ([BioME](https://bioinfo.imd.ufrn.br)), Federal University of Rio Grande do Norte, Brazil')
st.write('Upload [RING](https://ring.biocomputingup.it/) generated `.edges.txt` files to analyze network parameters.')

# File uploader
uploaded_files = st.file_uploader("Choose edge files...", accept_multiple_files=True, type="txt")

# Option to use example files
use_example_files = st.checkbox("...or use example files to evaluate SlytheRINs functionality.")

st.title('First Part: Network Analysis and Comparison Results')

# Load the edge file into a dataframe
edge_dataframes = []

if use_example_files:
    example_dir = 'example_data'
    example_files = [os.path.join(example_dir, f) for f in os.listdir(example_dir) if f.endswith('.txt')]
    for example_file in example_files:
        try:
            df = pd.read_csv(example_file, sep='\t')
            edge_dataframes.append(df)
        except EmptyDataError:
            st.error(f"The file {example_file} is empty or has no columns to parse.")
        except Exception as e:
            st.error(f"An error occurred while reading the file {example_file}: {e}")
else:
    for uploaded_file in uploaded_files:
        try:
            df = pd.read_csv(uploaded_file, sep='\t')
            edge_dataframes.append(df)
        except EmptyDataError:
            st.error(f"The file {uploaded_file.name} is empty or has no columns to parse.")
        except Exception as e:
            st.error(f"An error occurred while reading the file {uploaded_file.name}: {e}")
# Ensure there are valid dataframes to process
if edge_dataframes:
    # Initialize lists to store results
    all_degrees = []
    all_betweenness = []
    all_clustering = []
    all_triangles = []
    assortativity_values = []

    for df in edge_dataframes:
        # Calculate parameters for each dataframe
        G = nx.from_pandas_edgelist(df, 'NodeId1', 'NodeId2')
        degree = dict(G.degree())
        betweenness = nx.betweenness_centrality(G)
        clustering = nx.clustering(G)
        triangles = nx.triangles(G)
        assortativity = nx.degree_assortativity_coefficient(G)
        
        all_degrees.append(degree)
        all_betweenness.append(betweenness)
        all_clustering.append(clustering)
        all_triangles.append(triangles)
        assortativity_values.append(assortativity)

    # Calculate variance of degree values
    all_degree_values = pd.DataFrame(all_degrees).fillna(0)
    degree_variance = all_degree_values.var(axis=0)

    # Calculate variance of betweenness centrality values
    all_betweenness_values = pd.DataFrame(all_betweenness).fillna(0)
    betweenness_variance = all_betweenness_values.var(axis=0)

    # Calculate variance of clustering coefficient values
    all_clustering_values = pd.DataFrame(all_clustering).fillna(0)
    clustering_variance = all_clustering_values.var(axis=0)
    
    # Calculate the actual degree values for each residue across all edge files
    all_degree_values = pd.DataFrame(all_degrees).fillna(0)

    # Calculate mean and standard deviation of degree values
    mean_degree = all_degree_values.mean(axis=0)
    std_degree = all_degree_values.std(axis=0)
    neg_std_degree = -std_degree

    # Sort the standard deviation values by residue positions
    sorted_std_degree = degree_variance.sort_index()
    residue_numbers = sorted_std_degree.index
    sorted_residue_numbers = sorted(residue_numbers, key=lambda x: int(x.split(':')[1]))

    # Modify the x-axis labels to the format 'Position:Amino acid'
    formatted_labels = [f'{x.split(":")[1]}:{x.split(":")[3]}' for x in sorted_residue_numbers]
    
    # Convert triangle counts to a DataFrame
    triangle_counts_df = pd.DataFrame(all_triangles).fillna(0)
    mean_triangles = triangle_counts_df.mean(axis=0)
    std_triangles = triangle_counts_df.std(axis=0)
    
    # Create interactive plots using Plotly
    # Plot for degree variance
    fig_degree = go.Figure()
    fig_degree.add_trace(go.Scatter(x=sorted_residue_numbers, y=degree_variance.loc[sorted_residue_numbers].values, mode='lines+markers'))
    fig_degree.update_layout(title='Variance of Degree Values Across Residues',
                             xaxis_title='Residue Position and Amino Acid',
                             yaxis_title='Variance of Degree Values',
                             xaxis_tickangle=-90)
    st.plotly_chart(fig_degree)

    # Plot for betweenness centrality variance
    fig_betweenness = go.Figure()
    fig_betweenness.add_trace(go.Scatter(x=sorted_residue_numbers, y=betweenness_variance.loc[sorted_residue_numbers].values, mode='lines+markers'))
    fig_betweenness.update_layout(title='Variance of Betweenness Centrality Across Residues',
                                  xaxis_title='Residue Position and Amino Acid',
                                  yaxis_title='Variance of Betweenness Centrality',
                                  xaxis_tickangle=-90)
    st.plotly_chart(fig_betweenness)

    # Plot for clustering coefficient variance
    fig_clustering = go.Figure()
    fig_clustering.add_trace(go.Scatter(x=sorted_residue_numbers, y=clustering_variance.loc[sorted_residue_numbers].values, mode='lines+markers'))
    fig_clustering.update_layout(title='Variance of Clustering Coefficient Across Residues',
                                 xaxis_title='Residue Position and Amino Acid',
                                 yaxis_title='Variance of Clustering Coefficient',
                                 xaxis_tickangle=-90)
    st.plotly_chart(fig_clustering)
    
    # Plot the number of triangles for each node
    fig_triangles = go.Figure()
    fig_triangles.add_trace(go.Scatter(x=sorted_residue_numbers, y=mean_triangles.values, mode='lines+markers', name='Mean Triangles'))
    fig_triangles.add_trace(go.Scatter(x=sorted_residue_numbers, y=(mean_triangles - std_triangles).values, fill=None, mode='lines', line_color='lightblue', name='Lower Bound'))
    fig_triangles.add_trace(go.Scatter(x=sorted_residue_numbers, y=(mean_triangles + std_triangles).values, fill='tonexty', mode='lines', line_color='lightblue', name='Upper Bound'))
    fig_triangles.update_layout(title='Mean and Standard Deviation of Triangle Counts per Residue',
                                xaxis_title='Residue',
                                yaxis_title='Number of Triangles')
    st.plotly_chart(fig_triangles)

    # Plot the Mean Degree Value per Residue Position
    fig_degree_sd = go.Figure()
    fig_degree_sd.add_trace(go.Scatter(x=sorted_residue_numbers, y=mean_degree.loc[sorted_residue_numbers].values, mode='lines+markers', name='Mean Degree'))
    fig_degree_sd.add_trace(go.Scatter(x=sorted_residue_numbers, y=(mean_degree - std_degree).loc[sorted_residue_numbers].values, fill=None, mode='lines', line_color='lightblue', name='Lower Bound'))
    fig_degree_sd.add_trace(go.Scatter(x=sorted_residue_numbers, y=(mean_degree + std_degree).loc[sorted_residue_numbers].values, fill='tonexty', mode='lines', line_color='lightblue', name='Upper Bound'))
    fig_degree_sd.update_layout(title='Mean Degree Value per Residue Position',
                         xaxis_title='Residue Position and Amino Acid',
                         yaxis_title='Degree',
                         xaxis_tickangle=-90)
    st.plotly_chart(fig_degree_sd)

    # Plot the assortativity values
    fig_assortativity = go.Figure()
    fig_assortativity.add_trace(go.Scatter(x=list(range(1, len(assortativity_values) + 1)), y=assortativity_values, mode='lines+markers', name='Assortativity Coefficient'))
    fig_assortativity.update_layout(title='Assortativity Coefficient for Each Graph',
                                xaxis_title='Graph Index',
                                yaxis_title='Assortativity Coefficient')
    st.plotly_chart(fig_assortativity)
    
    # Create an interactive plot using Plotly with shaded regions
    fig_dsd = go.Figure()
    fig_dsd.add_trace(go.Scatter(x=sorted_residue_numbers, y=std_degree.values, mode='lines+markers', name='Standard Deviation'))
    fig_dsd.add_trace(go.Scatter(x=sorted_residue_numbers, y=neg_std_degree.values, mode='lines+markers', name='Negative Standard Deviation'))

    # Add shaded regions
    fig_dsd.add_shape(type='rect', x0=sorted_residue_numbers[0], x1=sorted_residue_numbers[-1], y0=-0.5, y1=0.5,
                         fillcolor='lightgrey', opacity=0.5, layer='below', line_width=0)
    fig_dsd.add_shape(type='rect', x0=sorted_residue_numbers[0], x1=sorted_residue_numbers[-1], y0=-1.0, y1=1.0,
                         fillcolor='lightgrey', opacity=0.3, layer='below', line_width=0)
    fig_dsd.update_layout(title='Standard Deviation and Negative Standard Deviation of Degree Values',
                             xaxis_title='Residue Position',
                             yaxis_title='Relative Standard Deviation of Degree')

    st.plotly_chart(fig_dsd)
    
    # Function to format NodeId1 in the same way as formatted_labels
def format_node_id(node_id):
    parts = node_id.split(':')
    return f'{parts[1]}:{parts[3]}'

# Count the occurrences of each interaction type per residue position for each graph
interaction_counts = {}

for df in edge_dataframes:
    if 'Interaction' in df.columns and 'NodeId1' in df.columns:
        # Count occurrences of each interaction type per NodeId1
        counts = df.groupby(['Interaction', 'NodeId1']).size().reset_index(name='Count')
        # Iterate over each interaction type
        for interaction in counts['Interaction'].unique():
            if interaction not in interaction_counts:
                interaction_counts[interaction] = []
            # Append the counts for each residue position
            interaction_counts[interaction].append(counts[counts['Interaction'] == interaction])
    else:
        st.error("The file does not contain the required columns 'Interaction' and 'NodeId1'.")

# Calculate the mean and standard deviation of the counts for each interaction type per residue
interaction_stats_counts = {}

# Iterate over each interaction type
for interaction, data_list in interaction_counts.items():
    # Concatenate all dataframes in the list
    combined_data = pd.concat(data_list)

    # Group by NodeId1 and calculate mean and std
    stats = combined_data.groupby('NodeId1')['Count'].agg(['mean', 'std']).reset_index()
    stats['FormattedNodeId1'] = stats['NodeId1'].apply(format_node_id)
    interaction_stats_counts[interaction] = stats

# Function to create interaction count plot
def create_interaction_count_plot(category, stats):
    figInt = go.Figure()
    figInt.add_trace(go.Scatter(
        x=stats['FormattedNodeId1'],
        y=stats['mean'],
        error_y=dict(type='data', array=stats['std']),
        mode='markers',
        marker=dict(size=10),
        name=category
    ))
    figInt.update_layout(
        title=f'{category} Interaction Counts',
        xaxis_title='Residue Position',
        yaxis_title='Mean Count',
        hovermode='closest'
    )
    return figInt

st.title('Second Part: Analysis of the Chemical Interactions')
st.write('**Legend for the following plots:** HBOND - Hydrogen bonds; SSBOND - Disulphide bridges; IONIC - Ionic bond; VDW - van der Waals; PICATION - π-cation; PIPISTACK - π-π stacking; PIHBOND - π-hydrogen; METAL_ION - Metal ion coordination; HALOGEN - Halogen bond; IAC - Inter-Atomic Contact (fallback type); MC - Main chain; SC - Side chain.')

# Create and display plots for each interaction category
for category in interaction_stats_counts.keys():
    st.header(f'{category} Interactions')
    figInt = create_interaction_count_plot(category, interaction_stats_counts[category])
    st.plotly_chart(figInt)
    
st.title('Third Part: Integration with AlphaMissense predicted mutation effects')
st.write('This module only works with human proteins, with an Uniprot entry number.')
st.write('Upload a .tsv file with the AlphaMissense predicted mutation effects.')
st.write('This file can be obtained from this grep command: ```grep P03891 AlphaMissense_aa_substitutions.tsv > P03891-pred.tsv```.')

# File uploader for AM exported .tsv files
uploaded_files_am = st.file_uploader("Choose AM exported .tsv files...", accept_multiple_files=False, type="tsv")

# Option to use example files
use_example_files_am = st.checkbox("...or use example files.")

# Initialize the dataframe
am_df = pd.DataFrame()

if use_example_files_am:
    example_dir = 'example_data'
    example_files_am = [os.path.join(example_dir, f) for f in os.listdir(example_dir) if f.endswith('.tsv')]
    for example_file_am in example_files_am:
        try:
            temp_df = pd.read_csv(example_file_am, sep='\t')
            am_df = pd.concat([am_df, temp_df], ignore_index=True)
        except EmptyDataError:
            st.error(f"The file {example_file_am} is empty or has no columns to parse.")
        except Exception as e:
            st.error(f"An error occurred while reading the file {example_file_am}: {e}")
else:
    if uploaded_files_am is not None:
        try:
            am_df = pd.read_csv(uploaded_files_am, sep='\t')
        except EmptyDataError:
            st.error(f"The file {uploaded_files_am.name} is empty or has no columns to parse.")
        except Exception as e:
            st.error(f"An error occurred while reading the file {uploaded_files_am.name}: {e}")

# Ensure the dataframe is not empty before proceeding
if not am_df.empty:
    # Check if the expected columns are present
    if am_df.shape[1] > 1:
        # Extract the position from the amino acid position column (column 1)
        am_df['Position'] = am_df.iloc[:, 1].str.extract('(\d+)').astype(int)

        # Count the occurrences of each classification for each position
        classification_counts = am_df.groupby(['Position', am_df.columns[3]]).size().unstack(fill_value=0)

        # Define the color map for Plotly
        color_map = {'pathogenic': 'orange', 'benign': 'green', 'ambiguous': 'blue'}

        # Create a stacked bar plot using Plotly
        figAM1 = go.Figure()
        for classification in classification_counts.columns:
            figAM1.add_trace(go.Bar(
                x=classification_counts.index,
                y=classification_counts[classification],
                name=classification,
                marker_color=color_map.get(classification, 'grey')
            ))

        figAM1.update_layout(
            barmode='stack',
            title='Classification Counts per Position',
            xaxis_title='Position',
            yaxis_title='Count',
            legend_title='Classification'
        )

        st.plotly_chart(figAM1)

        # Extract the position number from the second column
        am_df['Position'] = am_df.iloc[:, 1].str.extract('(\d+)').astype(int)
        
        # Set the third column as am_score
        am_df['am_score'] = am_df.iloc[:, 2]
        
        # Get the last residue position
        last_residue_position = am_df['Position'].max()

        # Create a box plot using Plotly
        figAM2 = px.box(am_df, x='Position', y='am_score', title='Predicted AlphaMissense Scores of possible amino acid changes in each position')

        # Add colored regions to the y-axis
        figAM2.add_shape(type='rect', x0=0, x1=last_residue_position, y0=0, y1=0.33, fillcolor='blue', opacity=0.1, line_width=0)
        figAM2.add_shape(type='rect', x0=0, x1=last_residue_position, y0=0.34, y1=0.56, fillcolor='gray', opacity=0.1, line_width=0)
        figAM2.add_shape(type='rect', x0=0, x1=last_residue_position, y0=0.56, y1=1.0, fillcolor='red', opacity=0.1, line_width=0)

        st.plotly_chart(figAM2)
    else:
        st.error("The file does not contain the expected columns.")
else:
    st.warning("No valid AM files were uploaded or found in the example data.")