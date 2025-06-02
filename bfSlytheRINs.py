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
import io # Required for BytesIO for image download
from collections import Counter # For degree distribution
import requests # For downloading AlphaMissense data
import re # For regex operations (extracting UniProt ID from filename)

# Helper function for downloads
@st.cache_data # Use Streamlit's caching for efficiency
def convert_df_to_tsv(df):
    """Converts a DataFrame to a TSV string for download."""
    return df.to_csv(sep='\t', index=False).encode('utf-8')

@st.cache_data
def convert_fig_to_png(fig):
    """Converts a Plotly figure to PNG bytes for download."""
    try:
        img_bytes = fig.to_image(format="png", scale=2) # Increase scale for better resolution
        return img_bytes
    except Exception as e:
        # Fallback or error message if Kaleido is not installed or fails
        st.error(f"Could not generate PNG. Ensure 'kaleido' is installed (pip install kaleido). Error: {e}")
        # Return an empty BytesIO object or None to prevent app crash
        return io.BytesIO()

# Page Configuration
st.set_page_config(layout="wide")


# Sidebar Setup
with st.sidebar:
    try:
        st.image('SlytheRINs-logo.svg', use_container_width=True)
    except FileNotFoundError:
        st.warning("SlytheRINs-logo.svg not found.")
    except Exception as e:
        st.warning(f"Could not load logo: {e}")

    st.title('**SlytheRINs**')
    st.write('*SlytheRINs is a Streamlit app to analyze and compare Residue Interaction Networks (RINs) calculated from different protein structures and conformations.*')
    st.caption('Developed by the [EvoMol-Lab](https://github.com/evomol-lab) - BioME, UFRN, Brazil')
    
    st.header("Data Upload")
    # Section for RIN data
    st.subheader("1. RIN Data")
    uploaded_files = st.file_uploader("Upload [RING](https://ring.biocomputingup.it/) generated .edges.txt files", accept_multiple_files=True, type="txt")
    use_example_files = st.checkbox("...or use example RIN files")

    # Section for AlphaMissense data
    st.subheader("2. AlphaMissense Data")
    st.info("The download of AlphaMissense data from the [AlphaFold DB](https://alphafold.ebi.ac.uk/) only works for human proteins.")
    uniprot_id_input = st.text_input("Enter UniProt Accession (e.g., P03891)", key="uniprot_id_input")
    download_am_button = st.button("Download & Process from AlphaFold DB", key="download_am_button")
    
    uploaded_file_am = st.file_uploader("...or upload AlphaMissense .tsv/.csv file", type=["tsv", "csv"], key="am_uploader_sidebar")
    use_example_file_am = st.checkbox("...or use an example AlphaMissense file", key="am_example_checkbox_sidebar")


# Data Loading Logic
edge_dataframes = []
edge_df_names = [] 

if use_example_files:
    example_dir = 'example_data'
    if os.path.isdir(example_dir):
        example_files_paths = [os.path.join(example_dir, f) for f in os.listdir(example_dir) if f.endswith('.txt')]
        if not example_files_paths:
            st.sidebar.warning(f"No '.txt' example files found in '{example_dir}'.")
        for example_file_path in example_files_paths:
            try:
                df = pd.read_csv(example_file_path, sep='\t')
                edge_dataframes.append(df)
                edge_df_names.append(os.path.basename(example_file_path))
            except EmptyDataError: st.sidebar.error(f"Example file {example_file_path} is empty.")
            except Exception as e: st.sidebar.error(f"Error reading example file {example_file_path}: {e}")
    else:
        st.sidebar.error(f"Example directory '{example_dir}' not found.")
elif uploaded_files: 
    for uploaded_file_obj in uploaded_files:
        try:
            df = pd.read_csv(uploaded_file_obj, sep='\t')
            edge_dataframes.append(df)
            edge_df_names.append(uploaded_file_obj.name)
        except EmptyDataError: st.sidebar.error(f"Uploaded file {uploaded_file_obj.name} is empty.")
        except Exception as e: st.sidebar.error(f"Error reading uploaded file {uploaded_file_obj.name}: {e}")

# AlphaMissense Data Loading (prioritized)
am_df_final = pd.DataFrame()
am_data_source_message = ""

# Initialize session state for downloaded AM data if it doesn't exist
if 'downloaded_am_data' not in st.session_state:
    st.session_state.downloaded_am_data = None
if 'downloaded_uniprot_id' not in st.session_state:
    st.session_state.downloaded_uniprot_id = None

# Define mapping for AlphaFold DB class values
am_class_mapping = {'LPath': 'pathogenic', 'LBen': 'benign', 'Amb': 'ambiguous'}

if download_am_button and uniprot_id_input:
    uniprot_id = uniprot_id_input.strip().upper()
    am_url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-aa-substitutions.csv"
    try:
        with st.spinner(f"Downloading AlphaMissense data for {uniprot_id}..."):
            response = requests.get(am_url, timeout=30) 
            response.raise_for_status()  
            csv_data = io.StringIO(response.text)
            temp_df = pd.read_csv(csv_data) # AlphaFold DB CSV has headers
            
            # Check for expected AlphaFold DB columns
            if 'protein_variant' in temp_df.columns and \
               'am_pathogenicity' in temp_df.columns and 'am_class' in temp_df.columns:
                
                am_df_final = temp_df.rename(columns={
                    'protein_variant': 'Mutation', 
                    'am_pathogenicity': 'AM_Score', 
                    'am_class': 'AM_Classification'
                })
                am_df_final['ProteinID'] = uniprot_id # Use the downloaded UniProt ID
                am_df_final = am_df_final[['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']] # Select and order
                
                # Map classification values
                am_df_final['AM_Classification'] = am_df_final['AM_Classification'].map(am_class_mapping).fillna(am_df_final['AM_Classification'])
                
                st.session_state.downloaded_am_data = am_df_final.copy() 
                st.session_state.downloaded_uniprot_id = uniprot_id
                am_data_source_message = f"Using AlphaMissense data downloaded for {uniprot_id}."
                st.sidebar.success(f"Successfully downloaded and processed data for {uniprot_id}.")
            else:
                st.sidebar.error(f"Downloaded CSV for {uniprot_id} has unexpected column names. Expected 'protein_variant', 'am_pathogenicity', 'am_class'.")
                st.session_state.downloaded_am_data = None
                st.session_state.downloaded_uniprot_id = None
    except requests.exceptions.HTTPError as e:
        st.sidebar.error(f"Could not download data for {uniprot_id}. HTTP Error: {e.response.status_code}. Check UniProt ID or try uploading.")
        st.session_state.downloaded_am_data = None
        st.session_state.downloaded_uniprot_id = None
    except requests.exceptions.RequestException as e:
        st.sidebar.error(f"Network error downloading data for {uniprot_id}: {e}")
        st.session_state.downloaded_am_data = None
        st.session_state.downloaded_uniprot_id = None
    except Exception as e:
        st.sidebar.error(f"Error processing downloaded AlphaMissense data for {uniprot_id}: {e}")
        st.session_state.downloaded_am_data = None
        st.session_state.downloaded_uniprot_id = None

# Use downloaded data if available and no new download was triggered for a different ID
if st.session_state.downloaded_am_data is not None and not st.session_state.downloaded_am_data.empty and \
   (not (download_am_button and uniprot_id_input) or (uniprot_id_input.strip().upper() == st.session_state.downloaded_uniprot_id)):
    am_df_final = st.session_state.downloaded_am_data
    if not am_data_source_message: # Ensure message is set if reloaded from session state
         am_data_source_message = f"Using previously downloaded AlphaMissense data for {st.session_state.downloaded_uniprot_id}."

elif uploaded_file_am:
    try:
        # Try reading as CSV first (AlphaFold DB format)
        temp_df = pd.read_csv(uploaded_file_am)
        uploaded_file_am.seek(0) # Reset pointer for potential re-read

        if 'protein_variant' in temp_df.columns and \
           'am_pathogenicity' in temp_df.columns and 'am_class' in temp_df.columns:
            am_df_final = temp_df.rename(columns={
                'protein_variant': 'Mutation', 
                'am_pathogenicity': 'AM_Score', 
                'am_class': 'AM_Classification'
            })
            # Try to extract UniProt ID from filename or set to Unknown
            match = re.search(r'AF-([A-Z0-9]+)-F1', uploaded_file_am.name, re.IGNORECASE)
            protein_id_from_filename = match.group(1) if match else "Unknown"
            am_df_final['ProteinID'] = protein_id_from_filename
            
            am_df_final = am_df_final[['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']] # Select and order
            am_df_final['AM_Classification'] = am_df_final['AM_Classification'].map(am_class_mapping).fillna(am_df_final['AM_Classification'])
            am_data_source_message = f"Using uploaded file: {uploaded_file_am.name} (Detected AlphaFold DB CSV format)."
            st.session_state.downloaded_am_data = None # Clear any previously downloaded data
            st.session_state.downloaded_uniprot_id = None
        else: # Assume it's the 4-column TSV/CSV format (grep output)
            temp_df = pd.read_csv(uploaded_file_am, sep=None, header=None, engine='python', on_bad_lines='skip') # sep=None tries to infer
            uploaded_file_am.seek(0)
            if temp_df.shape[1] == 4:
                 am_df_final = temp_df.copy()
                 am_df_final.columns = ['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']
                 # No need to map classification here, assume it's already 'pathogenic', 'benign', etc.
                 am_data_source_message = f"Using uploaded file: {uploaded_file_am.name} (Assumed 4-column format)."
                 st.session_state.downloaded_am_data = None 
                 st.session_state.downloaded_uniprot_id = None
            else:
                st.sidebar.error(f"Uploaded file {uploaded_file_am.name} is not in a recognized AlphaMissense format (neither AlphaFold DB CSV nor 4-column format).")
                am_df_final = pd.DataFrame()
    except Exception as e: 
        st.sidebar.error(f"Error reading uploaded AlphaMissense file {uploaded_file_am.name}: {e}")
        am_df_final = pd.DataFrame()
elif use_example_file_am:
    am_example_path = 'example_data/P03891-pred.tsv' # This is likely a grep output
    if os.path.exists(am_example_path):
        try:
            am_df_final = pd.read_csv(am_example_path, sep='\t', header=None) 
            if am_df_final.shape[1] == 4:
                am_df_final.columns = ['ProteinID', 'Mutation', 'AM_Score', 'AM_Classification']
                am_data_source_message = "Using example AlphaMissense file (P03891-pred.tsv)."
                st.session_state.downloaded_am_data = None
                st.session_state.downloaded_uniprot_id = None
            else:
                st.sidebar.error(f"Example AlphaMissense file {am_example_path} does not have 4 columns.")
                am_df_final = pd.DataFrame()
        except Exception as e: st.sidebar.error(f"Error reading example AlphaMissense file {am_example_path}: {e}")
    else: st.sidebar.warning(f"Example AlphaMissense file '{am_example_path}' not found in 'example_data' folder.")


# Main Dashboard Area
st.header('**SlytheRINs Dashboard**')

# Create tabs for different analysis parts
tab1, tab2, tab3 = st.tabs(["RIN Comparison & Network Complexity", "Chemical Interactions", "AlphaMissense Integration"])

# Define these flags outside the tab1 specific block if tab2 needs them
# However, it's better if tab2 relies on edge_dataframes directly
# and then does its own specific checks.
# For now, we'll define them here for tab1.
all_nodes_set_tab1 = set() 
valid_dfs_for_processing_tab1 = True 


with tab1:
    st.subheader('RIN Comparison & Network Complexity')
    if edge_dataframes:
        all_degrees_list = [] 
        all_graphs_nx = []    
        all_degrees_tab1, all_betweenness_tab1, all_clustering_tab1, all_triangles_tab1 = [], [], [], [] 
        assortativity_values_tab1 = []
        # all_nodes_set_tab1 = set() # Already defined above, or should be local to this block
        # valid_dfs_for_processing_tab1 = True # Already defined above

        for i, df_iter in enumerate(edge_dataframes): 
            if 'NodeId1' not in df_iter.columns or 'NodeId2' not in df_iter.columns:
                st.error(f"File '{edge_df_names[i]}' is missing 'NodeId1' or 'NodeId2' columns. Halting analysis for this part.")
                valid_dfs_for_processing_tab1 = False
                break
            try:
                G = nx.from_pandas_edgelist(df_iter, 'NodeId1', 'NodeId2')
                all_graphs_nx.append(G) 
                all_nodes_set_tab1.update(G.nodes())
                current_degrees = dict(G.degree())
                all_degrees_list.append(current_degrees) 
                all_degrees_tab1.append(current_degrees)      
                all_betweenness_tab1.append(nx.betweenness_centrality(G))
                all_clustering_tab1.append(nx.clustering(G)) 
                all_triangles_tab1.append(nx.triangles(G))
                try:
                    assortativity = nx.degree_assortativity_coefficient(G)
                    assortativity_values_tab1.append(assortativity)
                except (ZeroDivisionError, nx.NetworkXError): 
                    st.warning(f"Could not calculate assortativity for graph '{edge_df_names[i]}'. Assigning NaN.")
                    assortativity_values_tab1.append(np.nan)
            except Exception as e:
                st.error(f"Error processing graph from file '{edge_df_names[i]}': {e}")
                valid_dfs_for_processing_tab1 = False
                break

        if valid_dfs_for_processing_tab1 and all_nodes_set_tab1:
            sorted_residue_numbers = sorted(list(all_nodes_set_tab1), key=lambda x: int(x.split(':')[1]))
            formatted_labels = [f'{x.split(":")[1]}:{x.split(":")[3]}' for x in sorted_residue_numbers]
            node_to_formatted_label = dict(zip(sorted_residue_numbers, formatted_labels))
            
            all_degree_values = pd.DataFrame(all_degrees_tab1).reindex(columns=sorted_residue_numbers).fillna(0)
            all_betweenness_values = pd.DataFrame(all_betweenness_tab1).reindex(columns=sorted_residue_numbers).fillna(0)
            all_clustering_values = pd.DataFrame(all_clustering_tab1).reindex(columns=sorted_residue_numbers).fillna(0)
            triangle_counts_df = pd.DataFrame(all_triangles_tab1).reindex(columns=sorted_residue_numbers).fillna(0)
            
            degree_variance = all_degree_values.var(axis=0)
            betweenness_variance = all_betweenness_values.var(axis=0)
            clustering_variance = all_clustering_values.var(axis=0)
            mean_degree = all_degree_values.mean(axis=0)
            std_degree = all_degree_values.std(axis=0)
            neg_std_degree = -std_degree
            mean_triangles = triangle_counts_df.mean(axis=0)
            std_triangles = triangle_counts_df.std(axis=0)
            N_graphs = len(all_degree_values) 

            st.metric(label="Number of Networks Loaded", value=N_graphs)
            st.metric(label="Total Unique Residues Identified", value=len(all_nodes_set_tab1))

            df_degree_variance_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Degree_Variance': degree_variance.loc[sorted_residue_numbers]})
            df_betweenness_variance_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Betweenness_Variance': betweenness_variance.loc[sorted_residue_numbers]})
            df_clustering_variance_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Clustering_Variance': clustering_variance.loc[sorted_residue_numbers]})
            df_mean_std_degree_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Mean_Degree': mean_degree.loc[sorted_residue_numbers], 'STD_Degree': std_degree.loc[sorted_residue_numbers]})
            df_triangle_counts_dl = pd.DataFrame({'Residue': sorted_residue_numbers, 'Formatted_Label': formatted_labels, 'Mean_Triangles': mean_triangles.loc[sorted_residue_numbers], 'STD_Triangles': std_triangles.loc[sorted_residue_numbers]})
            df_assortativity_dl = pd.DataFrame({'Graph_Index': list(range(1, len(assortativity_values_tab1) + 1)), 
                                                'Graph_Name': edge_df_names[:len(assortativity_values_tab1)], 
                                                'Assortativity': assortativity_values_tab1})
            
            st.markdown("#### Variance Plots (Across All Networks)")
            col_var1, col_var2, col_var3 = st.columns(3)
            with col_var1:
                st.write("##### Degree Variance")
                st.write("*Variation in connectivity (degree) of residues across conformations. High variance highlights residues with constant change in the number of interactions.*")
                fig_degree = go.Figure(go.Scatter(x=formatted_labels, y=degree_variance.loc[sorted_residue_numbers], mode='lines+markers'))
                fig_degree.update_layout(title='Degree Variance', xaxis_title='Residue', yaxis_title='Variance', xaxis_tickangle=-90, height=400)
                st.plotly_chart(fig_degree, use_container_width=True)
                dl_c1, dl_c2 = st.columns(2)
                with dl_c1:
                    st.download_button("Data (TSV)", convert_df_to_tsv(df_degree_variance_dl), "degree_variance.tsv", "text/tab-separated-values", key="dl_deg_var_tsv_tab1")
                with dl_c2:
                    st.download_button("Plot (PNG)", convert_fig_to_png(fig_degree), "degree_variance_plot.png", "image/png", key="dl_deg_var_png_tab1")

            with col_var2:
                st.write("##### Betweenness-centrality Variance")
                st.write("*Variation in betweenness centrality of residues across conformations. High variance indicates residues that frequently change their role as bridges or bottlenecks between other residues in protein networks.*")
                fig_betweenness = go.Figure(go.Scatter(x=formatted_labels, y=betweenness_variance.loc[sorted_residue_numbers], mode='lines+markers'))
                fig_betweenness.update_layout(title='Betweenness Variance', xaxis_title='Residue', yaxis_title='Variance', xaxis_tickangle=-90, height=400)
                st.plotly_chart(fig_betweenness, use_container_width=True)
                dl_c1, dl_c2 = st.columns(2)
                with dl_c1:
                    st.download_button("Data (TSV)", convert_df_to_tsv(df_betweenness_variance_dl), "betweenness_variance.tsv", "text/tab-separated-values", key="dl_bet_var_tsv_tab1")
                with dl_c2:
                    st.download_button("Plot (PNG)", convert_fig_to_png(fig_betweenness), "betweenness_variance_plot.png", "image/png", key="dl_bet_var_png_tab1")

            with col_var3:
                st.write("##### Clustering Coefficient Variance")
                st.write("*Variation in clustering coefficient of residues across conformations. High variance indicates residues that frequently change their local connectivity patterns.*")
                fig_clustering = go.Figure(go.Scatter(x=formatted_labels, y=clustering_variance.loc[sorted_residue_numbers], mode='lines+markers'))
                fig_clustering.update_layout(title='Clustering Coeff. Variance', xaxis_title='Residue', yaxis_title='Variance', xaxis_tickangle=-90, height=400)
                st.plotly_chart(fig_clustering, use_container_width=True)
                dl_c1, dl_c2 = st.columns(2)
                with dl_c1:
                    st.download_button("Data (TSV)", convert_df_to_tsv(df_clustering_variance_dl), "clustering_variance.tsv", "text/tab-separated-values", key="dl_clust_var_tsv_tab1")
                with dl_c2:
                    st.download_button("Plot (PNG)", convert_fig_to_png(fig_clustering), "clustering_variance_plot.png", "image/png", key="dl_clust_var_png_tab1")


            st.markdown("#### Mean & STD Plots (Across All Networks)")
            col_mean_std1, col_mean_std2 = st.columns(2)
            with col_mean_std1:
                st.write("##### Triangle Counts")
                st.write("*Display each residue’s average (mean) involvement in triangles (3 mutually connected residues) and its variability across conformations. Shaded areas display the standard deviation, indicating residue participation in local clusters.*")
                fig_triangles = go.Figure([
                    go.Scatter(x=formatted_labels, y=mean_triangles.loc[sorted_residue_numbers], mode='lines+markers', name='Mean'),
                    go.Scatter(x=formatted_labels, y=(mean_triangles - std_triangles).loc[sorted_residue_numbers], fill=None, mode='lines', line_color='lightblue', name='-STD', showlegend=False),
                    go.Scatter(x=formatted_labels, y=(mean_triangles + std_triangles).loc[sorted_residue_numbers], fill='tonexty', mode='lines', line_color='lightblue', name='+STD', showlegend=False)
                ])
                fig_triangles.update_layout(title='Mean & STD of Triangle Counts', xaxis_title='Residue', yaxis_title='Count', xaxis_tickangle=-90, height=400)
                st.plotly_chart(fig_triangles, use_container_width=True)
                dl_c1, dl_c2 = st.columns(2)
                with dl_c1:
                    st.download_button("Data (TSV)", convert_df_to_tsv(df_triangle_counts_dl), "triangle_counts.tsv", "text/tab-separated-values", key="dl_tri_counts_tsv_tab1")
                with dl_c2:
                    st.download_button("Plot (PNG)", convert_fig_to_png(fig_triangles), "triangle_counts_plot.png", "image/png", key="dl_tri_counts_png_tab1")

            with col_mean_std2:
                st.write("##### Degree Values")
                st.write("*Display each residue’s average (mean) degree value and variability across distinct conformations. Shaded areas display the standard deviation, indicating residue connectivity stability.*")
                fig_degree_sd = go.Figure([
                    go.Scatter(x=formatted_labels, y=mean_degree.loc[sorted_residue_numbers], mode='lines+markers', name='Mean'),
                    go.Scatter(x=formatted_labels, y=(mean_degree - std_degree).loc[sorted_residue_numbers], fill=None, mode='lines', line_color='lightblue', name='-STD', showlegend=False),
                    go.Scatter(x=formatted_labels, y=(mean_degree + std_degree).loc[sorted_residue_numbers], fill='tonexty', mode='lines', line_color='lightblue', name='+STD', showlegend=False)
                ])
                fig_degree_sd.update_layout(title='Mean & STD of Degree Values', xaxis_title='Residue', yaxis_title='Degree', xaxis_tickangle=-90, height=400)
                st.plotly_chart(fig_degree_sd, use_container_width=True)
                dl_c1, dl_c2 = st.columns(2)
                with dl_c1:
                    st.download_button("Data (TSV)", convert_df_to_tsv(df_mean_std_degree_dl), "mean_std_degree.tsv", "text/tab-separated-values", key="dl_mean_std_deg_tsv_tab1")
                with dl_c2:
                    st.download_button("Plot (PNG)", convert_fig_to_png(fig_degree_sd), "mean_std_degree_plot.png", "image/png", key="dl_mean_std_deg_png_tab1")

            st.markdown("#### Degree Distribution per Residue (Box Plot - Across All Networks)")
            df_degree_boxplot_data = all_degree_values.melt(var_name='Residue_Raw', value_name='Degree')
            df_degree_boxplot_data['Residue'] = df_degree_boxplot_data['Residue_Raw'].map(node_to_formatted_label)
            df_degree_boxplot_data['Residue_Raw_Cat'] = pd.Categorical(df_degree_boxplot_data['Residue_Raw'], categories=sorted_residue_numbers, ordered=True)
            df_degree_boxplot_data = df_degree_boxplot_data.sort_values('Residue_Raw_Cat')

            if not df_degree_boxplot_data.empty:
                fig_degree_boxplot = px.box(df_degree_boxplot_data, x='Residue', y='Degree', 
                                            title='Distribution of Degree Values per Residue')
                fig_degree_boxplot.update_layout(xaxis_title='Residue', yaxis_title='Degree', xaxis_tickangle=-90, height=500)
                st.plotly_chart(fig_degree_boxplot, use_container_width=True)
                col1_bp, col2_bp = st.columns(2)
                with col1_bp:
                    st.download_button("Data (TSV)", convert_df_to_tsv(df_degree_boxplot_data[['Residue', 'Degree', 'Residue_Raw']]), "degree_boxplot_data.tsv", "text/tab-separated-values", key="dl_deg_boxplot_tsv_tab1")
                with col2_bp:
                    st.download_button("Plot (PNG)", convert_fig_to_png(fig_degree_boxplot), "degree_boxplot.png", "image/png", key="dl_deg_boxplot_png_tab1")
            else:
                st.warning("Could not generate data for the degree distribution box plot.")
            
            st.markdown("#### Degree Standard Deviation (DSD) Plot")
            st.write("*A SlytheRINs' specific metric. It displays each residue's Degree Standard Deviation (DSD), highlighting how variable its number of interactions is across conformations. Zero or low DSD values (shaded regions of the plot) indicate stable connections, while high DSD values indicates dynamic behavior.* ")
            fig_dsd = go.Figure([
                go.Scatter(x=formatted_labels, y=std_degree.loc[sorted_residue_numbers], mode='lines+markers', name='STD'),
                go.Scatter(x=formatted_labels, y=neg_std_degree.loc[sorted_residue_numbers], mode='lines+markers', name='-STD')
            ])
            if formatted_labels: 
                fig_dsd.add_shape(type='rect', x0=formatted_labels[0], x1=formatted_labels[-1], y0=-0.5, y1=0.5, fillcolor='lightgrey', opacity=0.5, layer='below', line_width=0)
                fig_dsd.add_shape(type='rect', x0=formatted_labels[0], x1=formatted_labels[-1], y0=-1.0, y1=1.0, fillcolor='lightgrey', opacity=0.3, layer='below', line_width=0)
            fig_dsd.update_layout(title='Degree Standard Deviation (DSD)', xaxis_title='Residue', yaxis_title='STD', xaxis_tickangle=-90, height=500)
            st.plotly_chart(fig_dsd, use_container_width=True)
            st.download_button("Download DSD Plot (PNG)", convert_fig_to_png(fig_dsd), "dsd_plot.png", "image/png", key="dl_dsd_png_tab1")


            st.markdown("#### Assortativity & Significance (Per Network)")
            assort_valid = [v for v in assortativity_values_tab1 if not np.isnan(v)] 
            if len(assort_valid) > 1:
                stat_assort, p_assort = st_sci.ttest_1samp(assort_valid, popmean=0)
                signif_assort = p_assort < 0.05
                texto_assort = f"{'* ' if signif_assort else ''}p = {p_assort:.3e}"
                st.write(f"**Overall Assortativity T-Test (vs 0):** Mean = {np.mean(assort_valid):.4f} | t = {stat_assort:.3f} | p = {p_assort:.3e} {'(*p < 0.05*)' if signif_assort else ''}")
            else:
                st.warning("Not enough valid assortativity values for overall t-test (requires at least 2).")
                texto_assort = "N/A (t-test requires >1 value)"

            fig_assortativity = go.Figure(go.Scatter(x=list(range(1, len(assortativity_values_tab1) + 1)), y=assortativity_values_tab1, mode='lines+markers'))
            fig_assortativity.add_hline(y=0, line_dash="dash", annotation_text="zero")
            fig_assortativity.add_annotation(xref="paper", yref="paper", x=0.98, y=0.98, text=texto_assort, showarrow=False, font=dict(size=10))
            fig_assortativity.update_layout(title="Assortativity Coefficient per Network", xaxis_title='Network Index', yaxis_title='Assortativity', height=400) 
            st.plotly_chart(fig_assortativity, use_container_width=True)
            col1, col2 = st.columns(2)
            with col1:
                st.download_button("Data (TSV)", convert_df_to_tsv(df_assortativity_dl), "assortativity.tsv", "text/tab-separated-values", key="dl_assort_tsv_tab1")
            with col2:
                st.download_button("Plot (PNG)", convert_fig_to_png(fig_assortativity), "assortativity_plot.png", "image/png", key="dl_assort_png_tab1")
            
            st.markdown("#### Significance Plots (Across All Networks)")
            col_sig1, col_sig2 = st.columns(2)
            with col_sig1:
                st.markdown("###### Degree Significance (Paired T-test vs. Others)")
                st.caption("Tests if a residue's degree is significantly different from the average degree of all *other* residues across the networks.")
                paired_t_p_values = []
                if N_graphs > 1:
                    for residue in sorted_residue_numbers:
                        x_degrees = all_degree_values[residue].values
                        overall_means_excluding = []
                        for idx, row in all_degree_values.iterrows():
                            if residue in row.index:
                                total = row.sum() - row[residue]
                                count = len(row) -1 if residue in row.index else len(row) 
                            else: 
                                total = row.sum()
                                count = len(row)
                            overall_means_excluding.append(total / count if count > 0 else 0)
                        if len(x_degrees) == len(overall_means_excluding) and len(x_degrees) > 1 and (np.std(x_degrees) > 0 or np.std(overall_means_excluding) > 0):
                            try:
                                t_stat, p_val_t = st_sci.ttest_rel(x_degrees, overall_means_excluding)
                                paired_t_p_values.append(p_val_t)
                            except ValueError: 
                                paired_t_p_values.append(np.nan)
                        else:
                            paired_t_p_values.append(np.nan) 

                    df_significance_paired_t = pd.DataFrame({"Residue": sorted_residue_numbers, "Formatted_Label": formatted_labels, "Degree_Variance": degree_variance.loc[sorted_residue_numbers], "Paired_T_p_value": paired_t_p_values})
                    fig_pvalues_paired = go.Figure(go.Scatter(x=formatted_labels, y=paired_t_p_values, mode='lines+markers'))
                    if formatted_labels:
                        fig_pvalues_paired.add_shape(type="line", x0=formatted_labels[0], x1=formatted_labels[-1], y0=0.05, y1=0.05, line=dict(color="red", dash="dash"))
                    fig_pvalues_paired.update_layout(title="P-Values (Paired T-Test)", xaxis_title="Residue", yaxis_title="P-Value", xaxis_tickangle=-90, height=400)
                    st.plotly_chart(fig_pvalues_paired, use_container_width=True)
                    dl_c1_sub, dl_c2_sub = st.columns(2)
                    with dl_c1_sub:
                        st.download_button("Data (TSV)", convert_df_to_tsv(df_significance_paired_t), "residue_significance_paired_t.tsv", "text/tab-separated-values", key="dl_paired_t_tsv_tab1")
                    with dl_c2_sub:
                        st.download_button("Plot (PNG)", convert_fig_to_png(fig_pvalues_paired), "paired_t_test_plot.png", "image/png", key="dl_paired_t_png_tab1")
                else:
                    st.warning("Need >1 network.")
            
            with col_sig2:
                st.markdown("###### Degree Variation Significance (vs. Poisson)")
                st.caption("Tests if a residue's degree variation *within itself* across networks is significantly different from a random (Poisson) process.")
                variation_p_values = []
                if N_graphs > 1:
                    for residue in sorted_residue_numbers:
                        degrees_for_residue = all_degree_values[residue].values
                        mean_deg_residue = np.mean(degrees_for_residue)
                        var_deg_residue = np.var(degrees_for_residue, ddof=1) 
                        if var_deg_residue == 0 or mean_deg_residue == 0:
                            p_value_chi2 = 1.0
                        else:
                            chi2_stat = (N_graphs - 1) * var_deg_residue / mean_deg_residue
                            cdf_val = st_sci.chi2.cdf(chi2_stat, N_graphs - 1)
                            p_value_chi2 = 2 * min(cdf_val, 1 - cdf_val) 
                        variation_p_values.append(p_value_chi2)

                    df_variation_pvals_dl = pd.DataFrame({"Residue": sorted_residue_numbers, "Formatted_Label": formatted_labels, "Degree_Variation_P_Value": variation_p_values})
                    fig_var_pvals = go.Figure(go.Scatter(x=formatted_labels, y=variation_p_values, mode='lines+markers'))
                    if formatted_labels:
                        fig_var_pvals.add_shape(type="line", x0=formatted_labels[0], x1=formatted_labels[-1], y0=0.05, y1=0.05, line=dict(color="red", dash="dash"))
                    fig_var_pvals.update_layout(title="P-Values (Degree Variation)", xaxis_title="Residue", yaxis_title="P-Value", xaxis_tickangle=-90, height=400)
                    st.plotly_chart(fig_var_pvals, use_container_width=True)
                    dl_c1_sub, dl_c2_sub = st.columns(2)
                    with dl_c1_sub:
                        st.download_button("Data (TSV)", convert_df_to_tsv(df_variation_pvals_dl), "degree_variation_significance.tsv", "text/tab-separated-values", key="dl_chi2_var_tsv_tab1")
                    with dl_c2_sub:
                        st.download_button("Plot (PNG)", convert_fig_to_png(fig_var_pvals), "degree_variation_plot.png", "image/png", key="dl_chi2_var_png_tab1")
                else:
                    st.warning("Need >1 network.")

            st.markdown("#### Network Complexity Analysis (Per Selected Network)")
            st.caption("This section provides metrics for a selected individual network to help assess its complexity (e.g., random vs. scale-free characteristics).")

            if not all_graphs_nx:
                st.warning("No graph objects were successfully created for complexity analysis.")
            else:
                graph_name_options = edge_df_names
                if graph_name_options:
                    selected_graph_name_complex = st.selectbox(
                        "Select a network to analyze for complexity:",
                        graph_name_options,
                        key="complexity_graph_selector"
                    )

                    if selected_graph_name_complex:
                        try:
                            selected_index = edge_df_names.index(selected_graph_name_complex)
                            G_complex = all_graphs_nx[selected_index] 
                            
                            st.markdown(f"##### Metrics for Network: `{selected_graph_name_complex}`")
                            
                            col_metric1, col_metric2, col_metric3 = st.columns(3)
                            num_nodes = G_complex.number_of_nodes()
                            num_edges = G_complex.number_of_edges()
                            density = nx.density(G_complex)
                            avg_clustering_graph = nx.average_clustering(G_complex)
                            
                            col_metric1.metric(label="Number of Nodes", value=num_nodes)
                            col_metric2.metric(label="Number of Edges", value=num_edges)
                            col_metric3.metric(label="Network Density", value=f"{density:.4f}")
                            col_metric1.metric(label="Avg. Clustering Coeff.", value=f"{avg_clustering_graph:.4f}")


                            if nx.is_connected(G_complex):
                                avg_shortest_path = nx.average_shortest_path_length(G_complex)
                                col_metric2.metric(label="Avg. Shortest Path", value=f"{avg_shortest_path:.4f}")
                            else:
                                if num_nodes > 0 and G_complex.number_of_edges() > 0: 
                                    largest_cc_nodes = max(nx.connected_components(G_complex), key=len)
                                    subgraph = G_complex.subgraph(largest_cc_nodes)
                                    if subgraph.number_of_nodes() > 1 : 
                                        avg_shortest_path = nx.average_shortest_path_length(subgraph)
                                        col_metric2.metric(label="Avg. Shortest Path (LCC)", value=f"{avg_shortest_path:.4f}")
                                        st.caption("*Avg. Shortest Path calculated on the Largest Connected Component (LCC) as graph is disconnected.*")

                                    else:
                                        col_metric2.metric(label="Avg. Shortest Path (LCC)", value="N/A")
                                else:
                                    col_metric2.metric(label="Avg. Shortest Path", value="N/A")
                            
                            st.markdown("###### Degree Distribution")
                            if num_nodes > 0:
                                degrees_current_graph = [d for n, d in G_complex.degree()]
                                if not degrees_current_graph:
                                    st.write("No degrees to plot for this graph (empty or no edges).")
                                else:
                                    degree_counts = Counter(degrees_current_graph)
                                    deg, cnt = zip(*degree_counts.items())
                                    
                                    df_degree_dist_dl = pd.DataFrame({'Degree': deg, 'Count': cnt}).sort_values(by='Degree')

                                    fig_deg_dist = go.Figure()
                                    fig_deg_dist.add_trace(go.Bar(x=deg, y=cnt))
                                    
                                    log_x_dd = st.checkbox(f"Log X-axis", key=f"logx_dd_{selected_graph_name_complex.replace('.', '_')}")
                                    log_y_dd = st.checkbox(f"Log Y-axis", key=f"logy_dd_{selected_graph_name_complex.replace('.', '_')}")

                                    xaxis_type_dd = "log" if log_x_dd else "linear"
                                    yaxis_type_dd = "log" if log_y_dd else "linear"
                                    
                                    fig_deg_dist.update_layout(
                                        title=f"Degree Distribution for '{selected_graph_name_complex}'",
                                        xaxis_title="Degree (k)",
                                        yaxis_title="Number of Nodes P(k)",
                                        bargap=0.1,
                                        xaxis_type=xaxis_type_dd,
                                        yaxis_type=yaxis_type_dd,
                                        height=400
                                    )
                                    
                                    st.plotly_chart(fig_deg_dist, use_container_width=True)
                                    dl_c1_sub, dl_c2_sub = st.columns(2)
                                    with dl_c1_sub:
                                        st.download_button(f"Data (TSV)", convert_df_to_tsv(df_degree_dist_dl), f"degree_distribution_{selected_graph_name_complex}.tsv", "text/tab-separated-values", key=f"dl_deg_dist_tsv_sel_{selected_graph_name_complex.replace('.', '_')}")
                                    with dl_c2_sub:
                                        st.download_button(f"Plot (PNG)", convert_fig_to_png(fig_deg_dist), f"degree_distribution_{selected_graph_name_complex}_plot.png", "image/png", key=f"dl_deg_dist_png_sel_{selected_graph_name_complex.replace('.', '_')}")
                            else:
                                st.write("Graph is empty, cannot compute degree distribution.")
                        except ValueError:
                            st.error(f"Could not find the selected graph '{selected_graph_name_complex}' in the processed list.")
                        except IndexError:
                             st.error(f"Index error trying to access graph '{selected_graph_name_complex}'.")
                else: 
                    st.info("No networks available to select for complexity analysis.")
            # --- END OF Network Complexity Analysis ---

        elif valid_dfs_for_processing_tab1 and not all_nodes_set_tab1: 
             st.warning("Dataframes were loaded, but no nodes could be extracted. Check your .edges.txt files for valid NodeId formats (e.g., Chain:Position:ResID:AA).")
        elif not edge_dataframes: 
            st.info("Upload .edges.txt files or select 'use example files' via the sidebar to see Part 1 results.")
        # else: (implicitly, if edge_dataframes exist but valid_dfs_for_processing is False, errors were already shown)
    else: # No edge_dataframes loaded
        st.info("Upload .edges.txt files or select 'use example files' via the sidebar to begin Part 1 analysis.")


with tab2:
    st.subheader('Chemical Interaction Analysis')
    # Check if edge_dataframes were loaded, and if the processing in tab1 was valid (valid_dfs_for_processing_tab1 and all_nodes_set_tab1)
    # This ensures that Part 2 only runs if Part 1 had meaningful data to process.
    if edge_dataframes and valid_dfs_for_processing_tab1 and all_nodes_set_tab1 : 
        st.write('**Legend examples:** HBOND - Hydrogen bonds; SSBOND - Disulphide bridges; IONIC - Ionic bond; VDW - van der Waals; etc.')

        def format_node_id_for_interactions(node_id):
            parts = node_id.split(':')
            if len(parts) >= 4:
                return f'{parts[1]}:{parts[3]}' 
            return node_id 

        interaction_counts_dict = {}
        has_interaction_data = False
        all_interaction_nodes_set = set()

        for df_iter in edge_dataframes: # Use the globally loaded edge_dataframes
            if 'Interaction' in df_iter.columns and 'NodeId1' in df_iter.columns:
                has_interaction_data = True
                valid_node_ids = df_iter[df_iter['NodeId1'].astype(str).str.count(':') == 3]['NodeId1']
                all_interaction_nodes_set.update(valid_node_ids.unique())
                counts = df_iter[df_iter['NodeId1'].astype(str).str.count(':') == 3].groupby(['Interaction', 'NodeId1']).size().reset_index(name='Count')
                for interaction_type in counts['Interaction'].unique():
                    if interaction_type not in interaction_counts_dict:
                        interaction_counts_dict[interaction_type] = []
                    interaction_counts_dict[interaction_type].append(counts[counts['Interaction'] == interaction_type])

        if not has_interaction_data:
            st.warning("None of the uploaded files suitable for Part 2 contain the 'Interaction' and 'NodeId1' columns, or NodeId1 format is incorrect. Skipping chemical interaction analysis.")
        elif not all_interaction_nodes_set:
            st.warning("Interaction data column found, but no valid NodeIds (e.g. A:123:ALA:N) could be extracted for interaction analysis. Skipping.")
        else:
            interaction_stats_counts_dict = {}
            if not all_interaction_nodes_set: # Should be redundant due to above check, but safe
                 st.warning("Interaction data found, but no nodes could be extracted for Part 2. Skipping analysis.")
            else:
                sorted_interaction_nodes = sorted(list(all_interaction_nodes_set), key=lambda x: int(x.split(':')[1]))
                for interaction_type, data_list in interaction_counts_dict.items():
                    if not data_list: continue
                    combined_data = pd.concat(data_list)
                    if 'NodeId1' not in combined_data.columns: continue
                    stats = combined_data.groupby('NodeId1')['Count'].agg(['mean', 'std']).reindex(sorted_interaction_nodes).fillna(0).reset_index()
                    stats['FormattedNodeId1'] = stats['NodeId1'].apply(format_node_id_for_interactions)
                    stats = stats.set_index('NodeId1').reindex(sorted_interaction_nodes).reset_index()
                    stats['FormattedNodeId1'] = stats['NodeId1'].apply(format_node_id_for_interactions) 
                    interaction_stats_counts_dict[interaction_type] = stats.fillna(0)

                def create_interaction_count_plot(category_name, stats_df):
                    figInt = go.Figure(go.Scatter(
                        x=stats_df['FormattedNodeId1'], y=stats_df['mean'],
                        error_y=dict(type='data', array=stats_df['std'].fillna(0)),
                        mode='markers', marker=dict(size=10), name=category_name
                    ))
                    figInt.update_layout(title=f'{category_name} Counts', xaxis_title='Residue', yaxis_title='Mean Count', xaxis_tickangle=-90, hovermode='closest', height=400)
                    return figInt

                if not interaction_stats_counts_dict:
                    st.warning("No interaction statistics could be calculated.")
                else:
                    num_interaction_cols = 2 
                    interaction_categories = list(interaction_stats_counts_dict.keys())
                    for i in range(0, len(interaction_categories), num_interaction_cols):
                        cols_interaction = st.columns(num_interaction_cols)
                        for j in range(num_interaction_cols):
                            if i + j < len(interaction_categories):
                                category_name = interaction_categories[i+j]
                                stats_df_for_plot = interaction_stats_counts_dict[category_name]
                                with cols_interaction[j]:
                                    # st.markdown(f"###### {category_name} Interactions") 
                                    if stats_df_for_plot.empty or 'FormattedNodeId1' not in stats_df_for_plot.columns or 'mean' not in stats_df_for_plot.columns:
                                        st.write(f"No data for {category_name}.")
                                        continue
                                    figInt = create_interaction_count_plot(category_name, stats_df_for_plot)
                                    st.plotly_chart(figInt, use_container_width=True)
                                    dl_c1_sub, dl_c2_sub = st.columns(2)
                                    with dl_c1_sub:
                                        st.download_button(
                                            f"Data (TSV)",
                                            convert_df_to_tsv(stats_df_for_plot),
                                            f"interaction_{category_name}_counts.tsv",
                                            "text/tab-separated-values",
                                            key=f"dl_int_tsv_{category_name}_tab2"
                                        )
                                    with dl_c2_sub:
                                        st.download_button(
                                            f"Plot (PNG)",
                                            convert_fig_to_png(figInt),
                                            f"interaction_{category_name}_plot.png",
                                            "image/png",
                                            key=f"dl_int_png_{category_name}_tab2"
                                        )
    else:
        st.info("Upload .edges.txt files via the sidebar and ensure they contain 'Interaction' and 'NodeId1' columns for Part 2 analysis.")


with tab3:
    st.subheader('AlphaMissense Integration')
    if am_data_source_message: 
        st.success(am_data_source_message)

    if not am_df_final.empty:
        try:
            am_df_final['Mutation'] = am_df_final['Mutation'].astype(str)
            am_df_final['Position_Extract'] = am_df_final['Mutation'].str.extract(r'[A-Za-z](\d+)[A-Za-z]')[0]
            am_df_final.dropna(subset=['Position_Extract'], inplace=True)
            am_df_final['Position'] = am_df_final['Position_Extract'].astype(int)


            am_df_final['AM_Score'] = pd.to_numeric(am_df_final['AM_Score'], errors='coerce')
            am_df_final.dropna(subset=['Position', 'AM_Score'], inplace=True)

            if am_df_final.empty:
                st.warning("No valid AlphaMissense data could be processed after type conversion and cleaning.")
            else:
                st.markdown("#### AlphaMissense Classification Counts")
                classification_counts = am_df_final.groupby(['Position', 'AM_Classification']).size().unstack(fill_value=0)
                color_map_am = {'pathogenic': 'orange', 'benign': 'green', 'ambiguous': 'blue'} 
                figAM1 = go.Figure()
                plot_order = [col for col in color_map_am.keys() if col in classification_counts.columns]
                for col in classification_counts.columns:
                    if col not in plot_order:
                        plot_order.append(col)
                for classification_type in plot_order:
                    if classification_type in classification_counts.columns: 
                        figAM1.add_trace(go.Bar(
                            name=classification_type,
                            x=classification_counts.index, 
                            y=classification_counts[classification_type],
                            marker_color=color_map_am.get(classification_type, 'grey') 
                        ))
                figAM1.update_layout(
                    barmode='stack', 
                    title='AlphaMissense Classification Counts per Position',
                    xaxis_title='Position', 
                    yaxis_title='Count', 
                    legend_title_text='AM Classification', height=500
                )
                st.plotly_chart(figAM1, use_container_width=True)
                col1, col2 = st.columns(2)
                with col1:
                    st.download_button("Data (TSV)", convert_df_to_tsv(classification_counts.reset_index()), "AlphaMissense_classification_counts.tsv", "text/tab-separated-values", key="dl_am_class_tsv_tab3")
                with col2:
                    st.download_button("Plot (PNG)", convert_fig_to_png(figAM1), "AlphaMissense_classification_plot.png", "image/png", key="dl_am_class_png_tab3")

                st.markdown("#### AlphaMissense Score Distribution")
                last_residue_position_am = am_df_final['Position'].max()
                figAM2 = px.box(am_df_final, x='Position', y='AM_Score', title='Predicted AlphaMissense Scores per Position')
                figAM2.add_shape(type='rect', x0=0, x1=last_residue_position_am, y0=0, y1=0.34, fillcolor='green', opacity=0.1, line_width=0, name='Likely Benign')
                figAM2.add_shape(type='rect', x0=0, x1=last_residue_position_am, y0=0.34, y1=0.564, fillcolor='grey', opacity=0.1, line_width=0, name='Ambiguous')
                figAM2.add_shape(type='rect', x0=0, x1=last_residue_position_am, y0=0.564, y1=1.0, fillcolor='red', opacity=0.1, line_width=0, name='Likely Pathogenic')
                figAM2.update_layout(height=500)
                st.plotly_chart(figAM2, use_container_width=True)
                col1, col2 = st.columns(2)
                with col1:
                    st.download_button("Data (TSV)", convert_df_to_tsv(am_df_final[['ProteinID', 'Mutation', 'Position', 'AM_Score', 'AM_Classification']]), "AlphaMissense_data.tsv", "text/tab-separated-values", key="dl_am_full_tsv_tab3")
                with col2:
                    st.download_button("Plot (PNG)", convert_fig_to_png(figAM2), "AlphaMissense_score_plot.png", "image/png", key="dl_am_full_png_tab3")
        except Exception as e:
            st.error(f"Error processing AlphaMissense file: {e}. Please check its format and content.")
            st.error("Expected format for downloaded/uploaded files: ProteinID, Mutation (e.g., A123G), AM_Score (numeric), AM_Classification (text).")

    elif uniprot_id_input or uploaded_file_am or use_example_file_am : 
        st.warning("Could not load or process AlphaMissense data. Please check the UniProt ID, the uploaded file, or the example file path/format via the sidebar.")
    else:
        st.info("Enter a UniProt ID and click 'Download & Process', upload an AlphaMissense .tsv/.csv file, or select the example file via the sidebar to see Part 3 integration results.")


# Final message if no data is loaded at all for RINs
if not edge_dataframes and not (uploaded_files or use_example_files): # This condition needs to be outside the tab structure
     st.info("Welcome to SlytheRINs! Please upload .edges.txt files or select 'use example files' via the sidebar to begin analysis.")
# End of the main Streamlit app code