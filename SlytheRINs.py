# Import necessary libraries
import streamlit as st
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from pandas.errors import EmptyDataError
import os
import scipy.stats as st_sci # Renamed to avoid conflict with streamlit 'st'
import numpy as np

# --- Helper function for downloads ---
@st.cache_data # Use Streamlit's caching for efficiency
def convert_df_to_tsv(df):
    """Converts a DataFrame to a TSV string for download."""
    return df.to_csv(sep='\t', index=False).encode('utf-8')

# --- Streamlit App Layout ---
# Attempt to load the logo, handle if not found
try:
    st.image('SlytheRINs-logo.svg', use_column_width=True)
except FileNotFoundError:
    st.warning("SlytheRINs-logo.svg not found. The logo will not be displayed.")
except Exception as e:
    st.warning(f"Could not load logo: {e}")


# Title of the app
st.title('SlytheRINs')

# Description
st.write('SlytheRINs is a Streamlit app to analyze and compare Residue Interaction Networks (RINs) calculated from different protein structures and conformations.')
st.write('Developed by the [EvoMol-Lab](https://github.com/evomol-lab) team at the Bioinformatics Multidisciplinary Environment ([BioME](https://bioinfo.imd.ufrn.br)), Federal University of Rio Grande do Norte, Brazil')
st.write('Upload [RING](https://ring.biocomputingup.it/) generated `.edges.txt` files to analyze network parameters.')

# --- File Uploader ---
uploaded_files = st.file_uploader("Choose .edges.txt files...", accept_multiple_files=True, type="txt")
use_example_files = st.checkbox("...or use example files to evaluate SlytheRINs functionality.")

# --- Data Loading ---
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
            except EmptyDataError: st.error(f"Example file {example_file} is empty.")
            except Exception as e: st.error(f"Error reading example file {example_file}: {e}")
    else:
        st.error(f"Example directory '{example_dir}' not found. Please create it or uncheck the box.")
else:
    for uploaded_file in uploaded_files:
        try:
            df = pd.read_csv(uploaded_file, sep='\t')
            edge_dataframes.append(df)
        except EmptyDataError: st.error(f"Uploaded file {uploaded_file.name} is empty.")
        except Exception as e: st.error(f"Error reading uploaded file {uploaded_file.name}: {e}")

# --- Part 1: Network Analysis ---
st.title('Part 1: Network Analysis and Comparison Results')

if edge_dataframes:
    all_degrees, all_betweenness, all_clustering, all_triangles = [], [], [], []
    assortativity_values = []
    all_nodes_set = set() # Use a more descriptive name
    valid_dfs_for_processing = True # Flag to control processing flow

    for i, df_iter in enumerate(edge_dataframes): # Use df_iter to avoid conflict
        if 'NodeId1' not in df_iter.columns or 'NodeId2' not in df_iter.columns:
            st.error(f"File {i+1} (original index) is missing 'NodeId1' or 'NodeId2' columns. Halting analysis for this part.")
            valid_dfs_for_processing = False
            break
        try:
            G = nx.from_pandas_edgelist(df_iter, 'NodeId1', 'NodeId2')
            all_nodes_set.update(G.nodes())
            all_degrees.append(dict(G.degree()))
            all_betweenness.append(nx.betweenness_centrality(G))
            all_clustering.append(nx.clustering(G))
            all_triangles.append(nx.triangles(G))
            try:
                assortativity = nx.degree_assortativity_coefficient(G)
                assortativity_values.append(assortativity)
            except (ZeroDivisionError, nx.NetworkXError): # Catch specific NetworkX errors too
                st.warning(f"Could not calculate assortativity for graph {i+1}. Assigning NaN.")
                assortativity_values.append(np.nan)
        except Exception as e:
            st.error(f"Error processing graph from file {i+1}: {e}")
            valid_dfs_for_processing = False
            break


    if valid_dfs_for_processing and all_nodes_set:
        # Sort nodes by position number for consistent ordering
        sorted_residue_numbers = sorted(list(all_nodes_set), key=lambda x: int(x.split(':')[1]))
        formatted_labels = [f'{x.split(":")[1]}:{x.split(":")[3]}' for x in sorted_residue_numbers]

        # Ensure consistent DataFrames by reindexing with all found nodes
        all_degree_values = pd.DataFrame(all_degrees).reindex(columns=sorted_residue_numbers).fillna(0)
        all_betweenness_values = pd.DataFrame(all_betweenness).reindex(columns=sorted_residue_numbers).fillna(0)
        all_clustering_values = pd.DataFrame(all_clustering).reindex(columns=sorted_residue_numbers).fillna(0)
        triangle_counts_df = pd.DataFrame(all_triangles).reindex(columns=sorted_residue_numbers).fillna(0)

        # Calculations
        degree_variance = all_degree_values.var(axis=0)
        betweenness_variance = all_betweenness_values.var(axis=0)
        clustering_variance = all_clustering_values.var(axis=0)
        mean_degree = all_degree_values.mean(axis=0)
        std_degree = all_degree_values.std(axis=0)
        neg_std_degree = -std_degree
        mean_triangles = triangle_counts_df.mean(axis=0)
        std_triangles = triangle_counts_df.std(axis=0)
        N_graphs = len(all_degree_values) # Number of graphs/dataframes successfully processed

        # --- DataFrames for Download ---
        df_degree_variance_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Degree_Variance': degree_variance.loc[sorted_residue_numbers]})
        df_betweenness_variance_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Betweenness_Variance': betweenness_variance.loc[sorted_residue_numbers]})
        df_clustering_variance_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Clustering_Variance': clustering_variance.loc[sorted_residue_numbers]})
        df_mean_std_degree_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Mean_Degree': mean_degree.loc[sorted_residue_numbers], 'STD_Degree': std_degree.loc[sorted_residue_numbers]})
        df_triangle_counts_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Mean_Triangles': mean_triangles.loc[sorted_residue_numbers], 'STD_Triangles': std_triangles.loc[sorted_residue_numbers]})
        df_assortativity_dl = pd.DataFrame({'Graph_Index': list(range(1, len(assortativity_values) + 1)), 'Assortativity': assortativity_values})

        # --- Plotting and Downloads ---
        st.subheader("Variance Plots")
        fig_degree = go.Figure(go.Scatter(x=formatted_labels, y=degree_variance.loc[sorted_residue_numbers], mode='lines+markers'))
        fig_degree.update_layout(title='Variance of Degree Values', xaxis_title='Residue', yaxis_title='Variance', xaxis_tickangle=-90)
        st.plotly_chart(fig_degree)
        st.download_button("Download Degree Variance Data", convert_df_to_tsv(df_degree_variance_dl), "degree_variance.tsv", "text/tab-separated-values", key="dl_deg_var")

        fig_betweenness = go.Figure(go.Scatter(x=formatted_labels, y=betweenness_variance.loc[sorted_residue_numbers], mode='lines+markers'))
        fig_betweenness.update_layout(title='Variance of Betweenness Centrality', xaxis_title='Residue', yaxis_title='Variance', xaxis_tickangle=-90)
        st.plotly_chart(fig_betweenness)
        st.download_button("Download Betweenness Variance Data", convert_df_to_tsv(df_betweenness_variance_dl), "betweenness_variance.tsv", "text/tab-separated-values", key="dl_bet_var")

        fig_clustering = go.Figure(go.Scatter(x=formatted_labels, y=clustering_variance.loc[sorted_residue_numbers], mode='lines+markers'))
        fig_clustering.update_layout(title='Variance of Clustering Coefficient', xaxis_title='Residue', yaxis_title='Variance', xaxis_tickangle=-90)
        st.plotly_chart(fig_clustering)
        st.download_button("Download Clustering Variance Data", convert_df_to_tsv(df_clustering_variance_dl), "clustering_variance.tsv", "text/tab-separated-values", key="dl_clust_var")

        st.subheader("Mean & STD Plots")
        fig_triangles = go.Figure([
            go.Scatter(x=formatted_labels, y=mean_triangles.loc[sorted_residue_numbers], mode='lines+markers', name='Mean'),
            go.Scatter(x=formatted_labels, y=(mean_triangles - std_triangles).loc[sorted_residue_numbers], fill=None, mode='lines', line_color='lightblue', name='-STD', showlegend=False),
            go.Scatter(x=formatted_labels, y=(mean_triangles + std_triangles).loc[sorted_residue_numbers], fill='tonexty', mode='lines', line_color='lightblue', name='+STD', showlegend=False)
        ])
        fig_triangles.update_layout(title='Mean & STD of Triangle Counts', xaxis_title='Residue', yaxis_title='Count', xaxis_tickangle=-90)
        st.plotly_chart(fig_triangles)
        st.download_button("Download Triangle Counts Data", convert_df_to_tsv(df_triangle_counts_dl), "triangle_counts.tsv", "text/tab-separated-values", key="dl_tri_counts")

        fig_degree_sd = go.Figure([
            go.Scatter(x=formatted_labels, y=mean_degree.loc[sorted_residue_numbers], mode='lines+markers', name='Mean'),
            go.Scatter(x=formatted_labels, y=(mean_degree - std_degree).loc[sorted_residue_numbers], fill=None, mode='lines', line_color='lightblue', name='-STD', showlegend=False),
            go.Scatter(x=formatted_labels, y=(mean_degree + std_degree).loc[sorted_residue_numbers], fill='tonexty', mode='lines', line_color='lightblue', name='+STD', showlegend=False)
        ])
        fig_degree_sd.update_layout(title='Mean & STD of Degree Values', xaxis_title='Residue', yaxis_title='Degree', xaxis_tickangle=-90)
        st.plotly_chart(fig_degree_sd)
        st.download_button("Download Mean/STD Degree Data", convert_df_to_tsv(df_mean_std_degree_dl), "mean_std_degree.tsv", "text/tab-separated-values", key="dl_mean_std_deg")

        fig_dsd = go.Figure([
            go.Scatter(x=formatted_labels, y=std_degree.loc[sorted_residue_numbers], mode='lines+markers', name='STD'),
            go.Scatter(x=formatted_labels, y=neg_std_degree.loc[sorted_residue_numbers], mode='lines+markers', name='-STD')
        ])
        if formatted_labels: # Ensure labels exist before adding shapes
            fig_dsd.add_shape(type='rect', x0=formatted_labels[0], x1=formatted_labels[-1], y0=-0.5, y1=0.5, fillcolor='lightgrey', opacity=0.5, layer='below', line_width=0)
            fig_dsd.add_shape(type='rect', x0=formatted_labels[0], x1=formatted_labels[-1], y0=-1.0, y1=1.0, fillcolor='lightgrey', opacity=0.3, layer='below', line_width=0)
        fig_dsd.update_layout(title='Degree Standard Deviation (DSD)', xaxis_title='Residue', yaxis_title='STD', xaxis_tickangle=-90)
        st.plotly_chart(fig_dsd)

        st.subheader("Assortativity & Significance")
        assort_valid = [v for v in assortativity_values if not np.isnan(v)]
        if len(assort_valid) > 1:
            stat_assort, p_assort = st_sci.ttest_1samp(assort_valid, popmean=0)
            signif_assort = p_assort < 0.05
            texto_assort = f"{'* ' if signif_assort else ''}p = {p_assort:.3e}"
            st.write(f"**Assortativity T-Test:** Mean = {np.mean(assort_valid):.4f} | t = {stat_assort:.3f} | p = {p_assort:.3e} {'(*p < 0.05*)' if signif_assort else ''}")
        else:
            st.warning("Not enough valid assortativity values for t-test (requires at least 2).")
            texto_assort = "N/A (t-test requires >1 value)"

        fig_assortativity = go.Figure(go.Scatter(x=list(range(1, len(assortativity_values) + 1)), y=assortativity_values, mode='lines+markers'))
        fig_assortativity.add_hline(y=0, line_dash="dash", annotation_text="zero")
        fig_assortativity.add_annotation(xref="paper", yref="paper", x=0.98, y=0.98, text=texto_assort, showarrow=False, font=dict(size=10))
        fig_assortativity.update_layout(title="Assortativity Coefficient & T-Test", xaxis_title='Graph Index', yaxis_title='Assortativity')
        st.plotly_chart(fig_assortativity)
        st.download_button("Download Assortativity Data", convert_df_to_tsv(df_assortativity_dl), "assortativity.tsv", "text/tab-separated-values", key="dl_assort")

        # --- Paired T-test for Degree (vs Others) ---
        st.subheader("Significance of Degree (Paired T-test vs. Others)")
        st.write("Tests if a residue's degree is significantly different from the average degree of all *other* residues across the networks.")
        paired_t_p_values = []
        if N_graphs > 1:
            for residue in sorted_residue_numbers:
                x_degrees = all_degree_values[residue].values
                overall_means_excluding = []
                for idx, row in all_degree_values.iterrows():
                    # Ensure residue is in row.columns before trying to access/subtract
                    if residue in row.index:
                        total = row.sum() - row[residue]
                        count = len(row) -1 if residue in row.index else len(row) # Adjust count if residue was present
                    else: # Should not happen if reindexed correctly, but as a safeguard
                        total = row.sum()
                        count = len(row)
                    overall_means_excluding.append(total / count if count > 0 else 0)
                
                # Ensure x_degrees and overall_means_excluding have same length and are not constant
                if len(x_degrees) == len(overall_means_excluding) and len(x_degrees) > 1 and (np.std(x_degrees) > 0 or np.std(overall_means_excluding) > 0):
                    try:
                        t_stat, p_val_t = st_sci.ttest_rel(x_degrees, overall_means_excluding)
                        paired_t_p_values.append(p_val_t)
                    except ValueError: # Handles cases like constant arrays if std check fails
                        paired_t_p_values.append(np.nan)
                else:
                    paired_t_p_values.append(np.nan) # Not enough data or constant data

            df_significance_paired_t = pd.DataFrame({"Residue": sorted_residue_numbers, "Formatted_Label": formatted_labels, "Degree_Variance": degree_variance.loc[sorted_residue_numbers], "Paired_T_p_value": paired_t_p_values})
            fig_pvalues_paired = go.Figure(go.Scatter(x=formatted_labels, y=paired_t_p_values, mode='lines+markers'))
            if formatted_labels:
                fig_pvalues_paired.add_shape(type="line", x0=formatted_labels[0], x1=formatted_labels[-1], y0=0.05, y1=0.05, line=dict(color="red", dash="dash"))
            fig_pvalues_paired.update_layout(title="P-Values per Residue (Paired T-Test vs. Others)", xaxis_title="Residue", yaxis_title="P-Value", xaxis_tickangle=-90)
            st.plotly_chart(fig_pvalues_paired)
            st.download_button("Download Paired T-Test Data", convert_df_to_tsv(df_significance_paired_t), "residue_significance_paired_t.tsv", "text/tab-separated-values", key="dl_paired_t")
        else:
            st.warning("Need more than one network for Paired T-test.")

        # --- Chi-Squared Test for Degree Variation (vs Poisson) ---
        st.subheader("Significance of Degree Variation (vs. Poisson)")
        st.write("Tests if a residue's degree variation *within itself* across networks is significantly different from a random (Poisson) process.")
        variation_p_values = []
        if N_graphs > 1:
            for residue in sorted_residue_numbers:
                degrees_for_residue = all_degree_values[residue].values
                mean_deg_residue = np.mean(degrees_for_residue)
                var_deg_residue = np.var(degrees_for_residue, ddof=1) # sample variance
                if var_deg_residue == 0 or mean_deg_residue == 0:
                    p_value_chi2 = 1.0
                else:
                    chi2_stat = (N_graphs - 1) * var_deg_residue / mean_deg_residue
                    cdf_val = st_sci.chi2.cdf(chi2_stat, N_graphs - 1)
                    p_value_chi2 = 2 * min(cdf_val, 1 - cdf_val) # two-tailed
                variation_p_values.append(p_value_chi2)

            df_variation_pvals_dl = pd.DataFrame({"Residue": sorted_residue_numbers, "Formatted_Label": formatted_labels, "Degree_Variation_P_Value": variation_p_values})
            fig_var_pvals = go.Figure(go.Scatter(x=formatted_labels, y=variation_p_values, mode='lines+markers'))
            if formatted_labels:
                fig_var_pvals.add_shape(type="line", x0=formatted_labels[0], x1=formatted_labels[-1], y0=0.05, y1=0.05, line=dict(color="red", dash="dash"))
            fig_var_pvals.update_layout(title="P-Values for Degree Variation (vs. Poisson)", xaxis_title="Residue", yaxis_title="P-Value", xaxis_tickangle=-90)
            st.plotly_chart(fig_var_pvals)
            st.download_button("Download Degree Variation P-Value Data", convert_df_to_tsv(df_variation_pvals_dl), "degree_variation_significance.tsv", "text/tab-separated-values", key="dl_chi2_var")
        else:
            st.warning("Need more than one network for Degree Variation test.")
    elif valid_dfs_for_processing and not all_nodes_set:
         st.warning("Dataframes were loaded, but no nodes could be extracted. Check your .edges.txt files for valid NodeId formats (e.g., Chain:Position:ResID:AA).")


    # --- Part 2: Chemical Interactions Analysis ---
    st.title('Part 2: Analysis of the Chemical Interactions')
    st.write('**Legend examples:** HBOND - Hydrogen bonds; SSBOND - Disulphide bridges; IONIC - Ionic bond; VDW - van der Waals; etc.')

    def format_node_id_for_interactions(node_id): # Potentially different formatting if needed
        parts = node_id.split(':')
        if len(parts) >= 4:
            return f'{parts[1]}:{parts[3]}' # Position:AA
        return node_id # Fallback

    interaction_counts_dict = {}
    has_interaction_data = False
    all_interaction_nodes_set = set()

    for df_iter in edge_dataframes: # Iterate again over successfully loaded dataframes
        if 'Interaction' in df_iter.columns and 'NodeId1' in df_iter.columns:
            has_interaction_data = True
            all_interaction_nodes_set.update(df_iter['NodeId1'].unique())
            # Ensure NodeId1 is string for splitting
            counts = df_iter[df_iter['NodeId1'].astype(str).str.contains(':')].groupby(['Interaction', 'NodeId1']).size().reset_index(name='Count')
            for interaction_type in counts['Interaction'].unique():
                if interaction_type not in interaction_counts_dict:
                    interaction_counts_dict[interaction_type] = []
                interaction_counts_dict[interaction_type].append(counts[counts['Interaction'] == interaction_type])

    if not has_interaction_data:
        st.warning("None of the uploaded files contain the 'Interaction' and 'NodeId1' columns, or NodeId1 format is incorrect. Skipping chemical interaction analysis.")
    elif not all_interaction_nodes_set:
        st.warning("Interaction data column found, but no valid NodeIds could be extracted for interaction analysis. Skipping.")
    else:
        interaction_stats_counts_dict = {}
        # Sort interaction nodes by position for consistent ordering
        sorted_interaction_nodes = sorted(list(all_interaction_nodes_set), key=lambda x: int(x.split(':')[1]) if isinstance(x, str) and ':' in x and x.split(':')[1].isdigit() else 0)


        for interaction_type, data_list in interaction_counts_dict.items():
            if not data_list: continue # Skip if no data for this interaction type
            combined_data = pd.concat(data_list)
            if 'NodeId1' not in combined_data.columns: continue

            # Calculate stats, reindex with all interaction nodes, fill NaNs
            stats = combined_data.groupby('NodeId1')['Count'].agg(['mean', 'std']).reindex(sorted_interaction_nodes).fillna(0).reset_index()
            stats['FormattedNodeId1'] = stats['NodeId1'].apply(format_node_id_for_interactions)
            
            # Ensure sorting for plotting by re-setting index and reindexing
            stats = stats.set_index('NodeId1').reindex(sorted_interaction_nodes).reset_index()
            stats['FormattedNodeId1'] = stats['NodeId1'].apply(format_node_id_for_interactions) # Re-apply after reindex
            interaction_stats_counts_dict[interaction_type] = stats.fillna(0)


        def create_interaction_count_plot(category_name, stats_df):
            figInt = go.Figure(go.Scatter(
                x=stats_df['FormattedNodeId1'], y=stats_df['mean'],
                error_y=dict(type='data', array=stats_df['std'].fillna(0)), # Ensure error array has no NaNs
                mode='markers', marker=dict(size=10), name=category_name
            ))
            figInt.update_layout(title=f'{category_name} Interaction Counts', xaxis_title='Residue', yaxis_title='Mean Count', xaxis_tickangle=-90, hovermode='closest')
            return figInt

        if not interaction_stats_counts_dict:
            st.warning("No interaction statistics could be calculated.")
        else:
            for category_name, stats_df_for_plot in interaction_stats_counts_dict.items():
                st.subheader(f'{category_name} Interactions')
                if stats_df_for_plot.empty or 'FormattedNodeId1' not in stats_df_for_plot.columns or 'mean' not in stats_df_for_plot.columns:
                    st.write(f"No data or incomplete data to plot for {category_name}.")
                    continue
                figInt = create_interaction_count_plot(category_name, stats_df_for_plot)
                st.plotly_chart(figInt)
                st.download_button(
                    f"Download {category_name} Data",
                    convert_df_to_tsv(stats_df_for_plot),
                    f"interaction_{category_name}_counts.tsv",
                    "text/tab-separated-values",
                    key=f"dl_int_{category_name}"
                )

    # --- Part 3: AlphaMissense Integration ---
    st.title('Part 3: Integration with AlphaMissense predicted mutation effects')
    st.write('Upload a .tsv file with AlphaMissense predictions (e.g., from `grep PXXXXX AlphaMissense_aa_substitutions.tsv`). The file should have 4 columns: ProteinID, Mutation (e.g. A123G), AM_Score, AM_Classification.')
    uploaded_file_am = st.file_uploader("Choose AM .tsv file...", type="tsv", key="am_uploader") # Unique key
    use_example_file_am = st.checkbox("...or use an example AlphaMissense file.", key="am_example_checkbox") # Unique key

    am_df_final = pd.DataFrame() # Use a more descriptive name
    if use_example_file_am:
        am_example_path = 'example_data/P03891-pred.tsv' # Define your example file name
        if os.path.exists(am_example_path):
            try:
                am_df_final = pd.read_csv(am_example_path, sep='\t', header=None)
                # Assume grep output: ProteinID, Mutation, AM_Score, AM_Classification
                if am_df_final.shape[1] == 4:
                    am_df_final.columns = ['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']
                else:
                    st.error(f"Example AlphaMissense file {am_example_path} does not have 4 columns as expected.")
                    am_df_final = pd.DataFrame() # Reset if format is wrong
            except Exception as e: st.error(f"Error reading example AlphaMissense file {am_example_path}: {e}")
        else: st.warning(f"Example AlphaMissense file '{am_example_path}' not found in 'example_data' folder.")
    elif uploaded_file_am:
        try:
            am_df_final = pd.read_csv(uploaded_file_am, sep='\t', header=None)
            if am_df_final.shape[1] == 4:
                 am_df_final.columns = ['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']
            else:
                st.error(f"Uploaded AlphaMissense file {uploaded_file_am.name} does not have 4 columns. Please ensure it's ProteinID, Mutation, AM_Score, AM_Classification.")
                am_df_final = pd.DataFrame() # Reset
        except Exception as e: st.error(f"Error reading uploaded AlphaMissense file {uploaded_file_am.name}: {e}")

    if not am_df_final.empty:
        try:
            # Ensure correct types and extract position
            am_df_final['Position'] = am_df_final['Mutation'].astype(str).str.extract(r'([A-Z])(\d+)([A-Z])')[1].astype(int) # Extract digits for position
            am_df_final['AM_Score'] = pd.to_numeric(am_df_final['AM_Score'], errors='coerce')
            am_df_final.dropna(subset=['Position', 'AM_Score'], inplace=True) # Drop rows where conversion failed

            if am_df_final.empty:
                st.warning("No valid AlphaMissense data could be processed after type conversion and cleaning.")
            else:
                st.subheader("AlphaMissense Classification Counts")
                # Ensure AM_Classification is string for grouping
                classification_counts = am_df_final.groupby(['Position', am_df_final['AM_Classification'].astype(str)]).size().unstack(fill_value=0)
                color_map_am = {'likely_pathogenic': 'red', 'likely_benign': 'green', 'ambiguous': 'grey'}
                
                # Create a new DataFrame for plotting to handle potentially missing classification columns
                plot_df_am_class = pd.DataFrame(index=classification_counts.index)
                for col in color_map_am.keys(): # Ensure all desired columns are present
                    if col in classification_counts.columns:
                        plot_df_am_class[col] = classification_counts[col]
                    else:
                        plot_df_am_class[col] = 0 # Add missing columns with zeros
                
                figAM1 = px.bar(plot_df_am_class, title='AlphaMissense Classification Counts per Position', color_discrete_map=color_map_am)
                figAM1.update_layout(barmode='stack', xaxis_title='Position', yaxis_title='Count', legend_title_text='AM Classification')
                st.plotly_chart(figAM1)
                st.download_button("Download AM Classification Counts", convert_df_to_tsv(classification_counts.reset_index()), "AlphaMissense_classification_counts.tsv", "text/tab-separated-values", key="dl_am_class")

                st.subheader("AlphaMissense Score Distribution")
                last_residue_position_am = am_df_final['Position'].max()
                figAM2 = px.box(am_df_final, x='Position', y='AM_Score', title='Predicted AlphaMissense Scores per Position')
                # AlphaMissense score interpretation thresholds
                figAM2.add_shape(type='rect', x0=0, x1=last_residue_position_am, y0=0, y1=0.34, fillcolor='green', opacity=0.1, line_width=0, name='Likely Benign')
                figAM2.add_shape(type='rect', x0=0, x1=last_residue_position_am, y0=0.34, y1=0.564, fillcolor='grey', opacity=0.1, line_width=0, name='Ambiguous')
                figAM2.add_shape(type='rect', x0=0, x1=last_residue_position_am, y0=0.564, y1=1.0, fillcolor='red', opacity=0.1, line_width=0, name='Likely Pathogenic')
                st.plotly_chart(figAM2)
                st.download_button("Download Full AM Data", convert_df_to_tsv(am_df_final), "AlphaMissense_data.tsv", "text/tab-separated-values", key="dl_am_full")

        except Exception as e:
            st.error(f"Error processing AlphaMissense file: {e}. Please check its format and content.")
            st.error("Expected format: ProteinID, Mutation (e.g., A123G), AM_Score (numeric), AM_Classification (text).")

    elif uploaded_file_am or use_example_file_am : # Only show if an attempt was made
        st.warning("Could not load or process AlphaMissense data. Please check the file or the example file path/format.")
    else:
        st.info("Upload an AlphaMissense .tsv file or select the example file to see the integration results.")

elif uploaded_files or use_example_files: # If files were selected but processing failed early
     st.error("No valid data could be loaded or processed from the provided .edges.txt files. Please check file format and content.")
else:
    st.info("Upload `.edges.txt` files or check the 'use example files' box to begin analysis.")

