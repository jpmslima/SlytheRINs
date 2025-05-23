# Import necessary libraries
import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from pandas.errors import EmptyDataError
import os
import scipy.stats as st_sci  # Renamed to avoid conflict with streamlit 'st'
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
    if os.path.isdir(example_dir):
        example_files = [os.path.join(example_dir, f) for f in os.listdir(example_dir) if f.endswith('.txt')]
        if not example_files:
            st.warning(f"No '.txt' example files found in the '{example_dir}' directory.")
        for example_file in example_files:
            try:
                df = pd.read_csv(example_file, sep='\t')
                edge_dataframes.append(df)
            except EmptyDataError:
                st.error(f"The file {example_file} is empty or has no columns to parse.")
            except Exception as e:
                st.error(f"An error occurred while reading the file {example_file}: {e}")
    else:
        st.error(f"Example directory '{example_dir}' not found. Please create it and add example files or uncheck the box.")

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
    all_nodes = set()

    for df in edge_dataframes:
        # Check if essential columns exist
        if 'NodeId1' not in df.columns or 'NodeId2' not in df.columns:
            st.error(f"One of the files is missing 'NodeId1' or 'NodeId2' columns.")
            edge_dataframes = [] # Clear dataframes to stop processing
            break
        
        G = nx.from_pandas_edgelist(df, 'NodeId1', 'NodeId2')
        all_nodes.update(G.nodes()) # Collect all nodes across all graphs
        
        degree = dict(G.degree())
        betweenness = nx.betweenness_centrality(G)
        clustering = nx.clustering(G)
        triangles = nx.triangles(G)
        
        try:
           assortativity = nx.degree_assortativity_coefficient(G)
        except ZeroDivisionError:
            st.warning("Could not calculate assortativity for one graph (likely due to insufficient edges or structure). Assigning NaN.")
            assortativity = np.nan # Assign NaN if calculation fails

        all_degrees.append(degree)
        all_betweenness.append(betweenness)
        all_clustering.append(clustering)
        all_triangles.append(triangles)
        assortativity_values.append(assortativity)

if edge_dataframes: # Re-check if processing should continue
    
    # Ensure all nodes are present in each calculation using reindex
    all_nodes = sorted(list(all_nodes), key=lambda x: int(x.split(':')[1]))

    all_degree_values = pd.DataFrame(all_degrees).reindex(columns=all_nodes).fillna(0)
    all_betweenness_values = pd.DataFrame(all_betweenness).reindex(columns=all_nodes).fillna(0)
    all_clustering_values = pd.DataFrame(all_clustering).reindex(columns=all_nodes).fillna(0)
    triangle_counts_df = pd.DataFrame(all_triangles).reindex(columns=all_nodes).fillna(0)

    # Calculate variance
    degree_variance = all_degree_values.var(axis=0)
    betweenness_variance = all_betweenness_values.var(axis=0)
    clustering_variance = all_clustering_values.var(axis=0)
    
    # Calculate mean and standard deviation
    mean_degree = all_degree_values.mean(axis=0)
    std_degree = all_degree_values.std(axis=0)
    neg_std_degree = -std_degree
    mean_triangles = triangle_counts_df.mean(axis=0)
    std_triangles = triangle_counts_df.std(axis=0)

    # Sort the residue numbers using the collected all_nodes list
    sorted_residue_numbers = all_nodes

    # Modify the x-axis labels to the format 'Position:Amino acid'
    formatted_labels = [f'{x.split(":")[1]}:{x.split(":")[3]}' for x in sorted_residue_numbers]
    
    # --- PLOTS ---

    # Plot for degree variance
    fig_degree = go.Figure()
    fig_degree.add_trace(go.Scatter(x=formatted_labels, y=degree_variance.loc[sorted_residue_numbers].values, mode='lines+markers'))
    fig_degree.update_layout(title='Variance of Degree Values Across Residues',
                             xaxis_title='Residue Position and Amino Acid',
                             yaxis_title='Variance of Degree Values',
                             xaxis_tickangle=-90)
    st.plotly_chart(fig_degree)

    # Plot for betweenness centrality variance
    fig_betweenness = go.Figure()
    fig_betweenness.add_trace(go.Scatter(x=formatted_labels, y=betweenness_variance.loc[sorted_residue_numbers].values, mode='lines+markers'))
    fig_betweenness.update_layout(title='Variance of Betweenness Centrality Across Residues',
                                  xaxis_title='Residue Position and Amino Acid',
                                  yaxis_title='Variance of Betweenness Centrality',
                                  xaxis_tickangle=-90)
    st.plotly_chart(fig_betweenness)

    # Plot for clustering coefficient variance
    fig_clustering = go.Figure()
    fig_clustering.add_trace(go.Scatter(x=formatted_labels, y=clustering_variance.loc[sorted_residue_numbers].values, mode='lines+markers'))
    fig_clustering.update_layout(title='Variance of Clustering Coefficient Across Residues',
                                 xaxis_title='Residue Position and Amino Acid',
                                 yaxis_title='Variance of Clustering Coefficient',
                                 xaxis_tickangle=-90)
    st.plotly_chart(fig_clustering)
    
    # Plot the number of triangles for each node
    fig_triangles = go.Figure()
    fig_triangles.add_trace(go.Scatter(x=formatted_labels, y=mean_triangles.loc[sorted_residue_numbers].values, mode='lines+markers', name='Mean Triangles'))
    fig_triangles.add_trace(go.Scatter(x=formatted_labels, y=(mean_triangles - std_triangles).loc[sorted_residue_numbers].values, fill=None, mode='lines', line_color='lightblue', name='Lower Bound'))
    fig_triangles.add_trace(go.Scatter(x=formatted_labels, y=(mean_triangles + std_triangles).loc[sorted_residue_numbers].values, fill='tonexty', mode='lines', line_color='lightblue', name='Upper Bound'))
    fig_triangles.update_layout(title='Mean and Standard Deviation of Triangle Counts per Residue',
                                xaxis_title='Residue',
                                yaxis_title='Number of Triangles',
                                xaxis_tickangle=-90)
    st.plotly_chart(fig_triangles)

    # Plot the Mean Degree Value per Residue Position
    fig_degree_sd = go.Figure()
    fig_degree_sd.add_trace(go.Scatter(x=formatted_labels, y=mean_degree.loc[sorted_residue_numbers].values, mode='lines+markers', name='Mean Degree'))
    fig_degree_sd.add_trace(go.Scatter(x=formatted_labels, y=(mean_degree - std_degree).loc[sorted_residue_numbers].values, fill=None, mode='lines', line_color='lightblue', name='Lower Bound'))
    fig_degree_sd.add_trace(go.Scatter(x=formatted_labels, y=(mean_degree + std_degree).loc[sorted_residue_numbers].values, fill='tonexty', mode='lines', line_color='lightblue', name='Upper Bound'))
    fig_degree_sd.update_layout(title='Mean Degree Value per Residue Position',
                         xaxis_title='Residue Position and Amino Acid',
                         yaxis_title='Degree',
                         xaxis_tickangle=-90)
    st.plotly_chart(fig_degree_sd)

    # Plot the assortativity values (handle NaNs)
    fig_assortativity = go.Figure()
    fig_assortativity.add_trace(go.Scatter(x=list(range(1, len(assortativity_values) + 1)), y=assortativity_values, mode='lines+markers', name='Assortativity Coefficient'))
    fig_assortativity.update_layout(title='Assortativity Coefficient for Each Graph',
                                xaxis_title='Graph Index',
                                yaxis_title='Assortativity Coefficient')
    st.plotly_chart(fig_assortativity)
    
    # Plot the DSD plot
    fig_dsd = go.Figure()
    fig_dsd.add_trace(go.Scatter(x=formatted_labels, y=std_degree.loc[sorted_residue_numbers].values, mode='lines+markers', name='Standard Deviation'))
    fig_dsd.add_trace(go.Scatter(x=formatted_labels, y=neg_std_degree.loc[sorted_residue_numbers].values, mode='lines+markers', name='Negative Standard Deviation'))

    # Add shaded regions
    fig_dsd.add_shape(type='rect', x0=formatted_labels[0], x1=formatted_labels[-1], y0=-0.5, y1=0.5,
                         fillcolor='lightgrey', opacity=0.5, layer='below', line_width=0)
    fig_dsd.add_shape(type='rect', x0=formatted_labels[0], x1=formatted_labels[-1], y0=-1.0, y1=1.0,
                         fillcolor='lightgrey', opacity=0.3, layer='below', line_width=0)
    fig_dsd.update_layout(title='Standard Deviation and Negative Standard Deviation of Degree Values',
                             xaxis_title='Residue Position',
                             yaxis_title='Relative Standard Deviation of Degree',
                             xaxis_tickangle=-90)

    st.plotly_chart(fig_dsd)

    # --- START: New Code for Chi-Squared Test ---
    st.subheader("Significance of Degree Variation (vs. Poisson)")
    st.write("This test checks if the variation of degree values for each residue across all networks is significantly different from what would be expected by a random (Poisson) process. A low p-value (< 0.05) suggests the variation is structurally or dynamically relevant.")

    variation_p_values = []
    N = len(all_degree_values)  # Number of conformations/networks

    if N <= 1:
        st.warning("Need more than one network to calculate variation significance.")
    else:
        for residue in sorted_residue_numbers:
            degrees = all_degree_values[residue].values
            mean_deg = np.mean(degrees)
            var_deg = np.var(degrees, ddof=1)  # Use sample variance

            if var_deg == 0 or mean_deg == 0:
                p_value = 1.0  # No variation or zero mean -> not significantly variable
            else:
                chi2_stat = (N - 1) * var_deg / mean_deg
                cdf_val = st_sci.chi2.cdf(chi2_stat, N - 1)
                p_value = 2 * min(cdf_val, 1 - cdf_val) # Two-tailed p-value

            variation_p_values.append(p_value)

        # Create a DataFrame with the results (optional: add download button)
        df_variation_pvals = pd.DataFrame({
            "Residue": sorted_residue_numbers,
            "Formatted_Label": formatted_labels,
            "Degree_Variation_P_Value": variation_p_values,
        })
        # st.dataframe(df_variation_pvals) # Optionally display the table

        # Generate the plot
        fig_var_pvals = go.Figure()
        fig_var_pvals.add_trace(go.Scatter(
            x=formatted_labels,
            y=variation_p_values,
            mode='lines+markers',
            name="Degree Variation P-Value (vs. Poisson)"
        ))

        # Add 0.05 significance line
        fig_var_pvals.add_shape(
            dict(
                type="line",
                x0=formatted_labels[0],
                x1=formatted_labels[-1],
                y0=0.05,
                y1=0.05,
                line=dict(color="red", dash="dash")
            )
        )
        fig_var_pvals.update_layout(
            title="P-Values for Degree Variation per Residue (vs. Poisson)",
            xaxis_title="Residue Position and Amino Acid",
            yaxis_title="P-Value",
            xaxis_tickangle=-90
        )
        st.plotly_chart(fig_var_pvals)
    # --- END: New Code for Chi-Squared Test ---


    # --- Chemical Interactions Analysis ---
    st.title('Second Part: Analysis of the Chemical Interactions')
    st.write('**Legend for the following plots:** HBOND - Hydrogen bonds; SSBOND - Disulphide bridges; IONIC - Ionic bond; VDW - van der Waals; PICATION - π-cation; PIPISTACK - π-π stacking; PIHBOND - π-hydrogen; METAL_ION - Metal ion coordination; HALOGEN - Halogen bond; IAC - Inter-Atomic Contact (fallback type); MC - Main chain; SC - Side chain.')

    # Function to format NodeId1 in the same way as formatted_labels
    def format_node_id(node_id):
        parts = node_id.split(':')
        return f'{parts[1]}:{parts[3]}'

    interaction_counts = {}
    has_interaction_data = False

    for df in edge_dataframes:
        if 'Interaction' in df.columns and 'NodeId1' in df.columns:
            has_interaction_data = True
            counts = df.groupby(['Interaction', 'NodeId1']).size().reset_index(name='Count')
            for interaction in counts['Interaction'].unique():
                if interaction not in interaction_counts:
                    interaction_counts[interaction] = []
                interaction_counts[interaction].append(counts[counts['Interaction'] == interaction])
        # else: # No need to show error for every file, one warning is enough
            # st.error("The file does not contain the required columns 'Interaction' and 'NodeId1'.")

    if not has_interaction_data:
        st.warning("None of the uploaded files contain the 'Interaction' column. Skipping chemical interaction analysis.")
    else:
        interaction_stats_counts = {}
        all_interaction_nodes = set()

        # Collect all nodes from interaction data
        for interaction, data_list in interaction_counts.items():
            for df_int in data_list:
                all_interaction_nodes.update(df_int['NodeId1'].unique())
        
        sorted_interaction_nodes = sorted(list(all_interaction_nodes), key=lambda x: int(x.split(':')[1]))
        formatted_interaction_labels = [format_node_id(x) for x in sorted_interaction_nodes]


        for interaction, data_list in interaction_counts.items():
            combined_data = pd.concat(data_list)
            # Group by NodeId1 and calculate mean and std, reindex to ensure all nodes are present
            stats = combined_data.groupby('NodeId1')['Count'].agg(['mean', 'std']).reindex(sorted_interaction_nodes).fillna(0).reset_index()
            stats['FormattedNodeId1'] = stats['NodeId1'].apply(format_node_id)
            # Ensure stats are sorted for plotting
            stats = stats.set_index('NodeId1').reindex(sorted_interaction_nodes).reset_index()
            stats['FormattedNodeId1'] = stats['NodeId1'].apply(format_node_id) # Re-apply after reindex
            interaction_stats_counts[interaction] = stats.fillna(0) # Fill NaNs again after reindex


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
                hovermode='closest',
                xaxis_tickangle=-90
            )
            return figInt

        # Create and display plots for each interaction category
        for category, stats_df in interaction_stats_counts.items():
            st.header(f'{category} Interactions')
            figInt = create_interaction_count_plot(category, stats_df)
            st.plotly_chart(figInt)

# --- AlphaMissense Integration ---
st.title('Third Part: Integration with AlphaMissense predicted mutation effects')
st.write('This module only works with human proteins, with an Uniprot entry number.')
st.write('Upload a .tsv file with the AlphaMissense predicted mutation effects.')
st.write('This file can be obtained from the [AlphaMissense website](https://alphamissense.hegelab.org/) or using commands like: ```grep P03891 AlphaMissense_aa_substitutions.tsv > P03891-pred.tsv```.')

# File uploader for AM exported .tsv files
uploaded_files_am = st.file_uploader("Choose AM exported .tsv files...", accept_multiple_files=False, type="tsv")

# Option to use example files
use_example_files_am = st.checkbox("...or use an example AlphaMissense file.")

# Initialize the dataframe
am_df = pd.DataFrame()

if use_example_files_am:
    example_dir = 'example_data'
    am_example_path = os.path.join(example_dir, 'P03891-pred.tsv') # Assuming a specific name for example
    if os.path.exists(am_example_path):
        try:
            am_df = pd.read_csv(am_example_path, sep='\t', header=None) # Assume no header for grep output
            # Add headers manually if they are missing
            am_df.columns = ['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']
        except EmptyDataError:
            st.error(f"The example file {am_example_path} is empty or has no columns to parse.")
        except Exception as e:
            st.error(f"An error occurred while reading the file {am_example_path}: {e}")
    else:
        st.warning(f"Example AlphaMissense file '{am_example_path}' not found.")

else:
    if uploaded_files_am is not None:
        try:
            am_df = pd.read_csv(uploaded_files_am, sep='\t', header=None) # Assume no header
            am_df.columns = ['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']
        except EmptyDataError:
            st.error(f"The file {uploaded_files_am.name} is empty or has no columns to parse.")
        except Exception as e:
            st.error(f"An error occurred while reading the file {uploaded_files_am.name}: {e}")

# Ensure the dataframe is not empty before proceeding
if not am_df.empty:
    try:
        # Check if the expected columns are present (now we've added them)
        if 'Mutation' in am_df.columns and 'AM_Classification' in am_df.columns and 'AM_Score' in am_df.columns:
            # Extract the position from the mutation column (e.g., A123C -> 123)
            am_df['Position'] = am_df['Mutation'].str.extract('(\d+)').astype(int)

            # --- Plot 1: Classification Counts ---
            classification_counts = am_df.groupby(['Position', 'AM_Classification']).size().unstack(fill_value=0)
            color_map = {'likely_pathogenic': 'red', 'likely_benign': 'green', 'ambiguous': 'grey'}

            figAM1 = go.Figure()
            for classification in classification_counts.columns:
                figAM1.add_trace(go.Bar(
                    x=classification_counts.index,
                    y=classification_counts[classification],
                    name=classification,
                    marker_color=color_map.get(classification, 'blue') # Default color
                ))

            figAM1.update_layout(
                barmode='stack',
                title='AlphaMissense Classification Counts per Position',
                xaxis_title='Position',
                yaxis_title='Count',
                legend_title='Classification'
            )
            st.plotly_chart(figAM1)

            # --- Plot 2: Box Plot Scores ---
            am_df['am_score'] = am_df['AM_Score'].astype(float) # Ensure score is float
            last_residue_position = am_df['Position'].max()

            figAM2 = px.box(am_df, x='Position', y='am_score', title='Predicted AlphaMissense Scores of possible amino acid changes in each position')

            # Add colored regions (Updated thresholds based on AM paper)
            figAM2.add_shape(type='rect', x0=0, x1=last_residue_position, y0=0, y1=0.34, fillcolor='green', opacity=0.1, line_width=0, name='Likely Benign')
            figAM2.add_shape(type='rect', x0=0, x1=last_residue_position, y0=0.34, y1=0.564, fillcolor='gray', opacity=0.1, line_width=0, name='Ambiguous')
            figAM2.add_shape(type='rect', x0=0, x1=last_residue_position, y0=0.564, y1=1.0, fillcolor='red', opacity=0.1, line_width=0, name='Likely Pathogenic')

            st.plotly_chart(figAM2)
        else:
            st.error("The AlphaMissense file does not seem to have the expected structure (ProteinID, Mutation, Score, Classification). Please check the file format or ensure it was generated correctly.")
    except Exception as e:
        st.error(f"An error occurred while processing the AlphaMissense file: {e}. Please check its format.")

elif uploaded_files_am or use_example_files_am:
    st.warning("Could not load AlphaMissense data. Please check the file or the example file path.")
else:
    st.info("Upload an AlphaMissense file or select the example file to see the integration results.")