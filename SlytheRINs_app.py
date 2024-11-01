# Create a new Python file for the Streamlit app
streamlit_app_code = '''
import streamlit as st
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns

# Title of the app
st.title('Network Analysis App')

# Description
st.write('This app calculates network parameters from edge files.')

# List all edge files in the current directory
edge_files = [file for file in os.listdir() if file.endswith('_edges.txt')]

# Display the list of edge files
st.write('Available edge files:', edge_files)

# Function to calculate network parameters
def calculate_network_parameters(edge_file):
    # Load the edge data
    df = pd.read_csv(edge_file, sep='\t')
    
    # Create a graph
    G = nx.from_pandas_edgelist(df, 'NodeId1', 'NodeId2')
    
    # Calculate network parameters (e.g., degree, centrality)
    degree = dict(G.degree())
    centrality = nx.betweenness_centrality(G)
    
    return degree, centrality

# Select an edge file to analyze
selected_file = st.selectbox('Select an edge file to analyze:', edge_files)

# Calculate and display network parameters
if selected_file:
    degree, centrality = calculate_network_parameters(selected_file)
    st.write('Degree:', degree)
    st.write('Centrality:', centrality)

'''

# Save the Streamlit app code to a file
with open('streamlit_app.py', 'w') as file:
    file.write(streamlit_app_code)

print('Streamlit app structure created and saved as streamlit_app.py')

