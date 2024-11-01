import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns

# List all edge files in the current directory
edge_files = [file for file in os.listdir() if file.endswith('_edges.txt')]

# Function to calculate network parameters
def calculate_network_parameters(edge_file):
    # Load the edge data
    df = pd.read_csv(edge_file, sep='\t')
    
    # Create a graph
    G = nx.from_pandas_edgelist(df, 'NodeId1', 'NodeId2')
    
    # Calculate parameters
    degree = dict(G.degree())
    betweenness = nx.betweenness_centrality(G)
    clustering = nx.clustering(G)
    
    return degree, betweenness, clustering

# Initialize lists to store results
all_degrees = []
all_betweenness = []
all_clustering = []

# Calculate parameters for each file
for file in edge_files:
    degree, betweenness, clustering = calculate_network_parameters(file)
    all_degrees.append(degree)
    all_betweenness.append(betweenness)
    all_clustering.append(clustering)

# Calculate variance of degree values
all_degree_values = pd.DataFrame(all_degrees).fillna(0)
degree_variance = all_degree_values.var(axis=0)

# Calculate the standard deviation of the degree values
std_degree = all_degree_values.std(axis=0)

# Sort the standard deviation values by residue positions
sorted_std_degree = std_degree.sort_index()
neg_sorted_std_degree = -sorted_std_degree

# Extract and sort the residue numbers from the index
residue_numbers = sorted_std_degree.index
sorted_residue_numbers = sorted(residue_numbers, key=lambda x: int(x.split(':')[1]))

# Display the sorted residue numbers to verify
print(sorted_residue_numbers[:10])

# Modify the x-axis labels to the format 'Position:Amino acid'
# Extract and format the labels from sorted_residue_numbers
formatted_labels = [f'{x.split(":")[1]}:{x.split(":")[3]}' for x in sorted_residue_numbers]

# Display the first few formatted labels to verify
print(formatted_labels[:10])

# Plot variance of degree values with updated x-axis labels
plt.figure(figsize=(36, 12))
plt.plot(formatted_labels, degree_variance.loc[sorted_residue_numbers].values, marker='o')
plt.xticks(rotation=90, fontsize=6)  # Set x-axis labels to vertical
plt.xlabel('Residue Position and Amino Acid')
plt.ylabel('Variance of Degree Values')
plt.title('Variance of Degree Values Across Residues')
plt.autoscale(enable=True, axis='x', tight=True)
plt.savefig('P03891-Fig1.svg', format='svg')
plt.show()

# Calculate variance of betweenness centrality values
all_betweenness_values = pd.DataFrame(all_betweenness).fillna(0)
betweenness_variance = all_betweenness_values.var(axis=0)

# Plot variance of betweenness centrality values
plt.figure(figsize=(36, 12))
plt.plot(formatted_labels, betweenness_variance.values, marker='o')
plt.xticks(rotation=90, fontsize=6)  # Set x-axis labels to vertical
plt.xlabel('Residue Position')
plt.ylabel('Variance of Betweenness Centrality')
plt.title('Variance of Betweenness Centrality Across Residues')
plt.autoscale(enable=True, axis='x', tight=True)
plt.savefig('P03891-Fig2.svg', format='svg')
plt.show()

# Calculate variance of clustering coefficient values
all_clustering_values = pd.DataFrame(all_clustering).fillna(0)
clustering_variance = all_clustering_values.var(axis=0)

# Plot variance of clustering coefficient values
plt.figure(figsize=(36, 12))
plt.plot(formatted_labels, clustering_variance.values, marker='o')
plt.xticks(rotation=90, fontsize=6)  # Set x-axis labels to vertical
plt.xlabel('Residue Position')
plt.ylabel('Variance of Clustering Coefficient')
plt.title('Variance of Clustering Coefficient Across Residues')
plt.autoscale(enable=True, axis='x', tight=True)
plt.savefig('P03891-Fig3.svg', format='svg')
plt.show()

# Calculate the number of triangles for each network
triangle_counts = []

for file in edge_files:
    # Load the edge data
    df = pd.read_csv(file, sep='\t')
    
    # Create a graph
    G = nx.from_pandas_edgelist(df, 'NodeId1', 'NodeId2')
    
    # Calculate the number of triangles
    triangles = nx.triangles(G)
    triangle_counts.append(triangles)

# Convert triangle counts to a DataFrame
triangle_counts_df = pd.DataFrame(triangle_counts).fillna(0)

# Plot the number of triangles for each residue position
plt.figure(figsize=(48, 12))
plt.plot(formatted_labels, triangle_counts_df.sum(axis=0), marker='o')
plt.xticks(rotation=90, fontsize=6)  # Set x-axis labels to vertical
plt.xlabel('Residue Position')
plt.ylabel('Number of Triangles')
plt.title('Number of Triangles Across Residues')
plt.autoscale(enable=True, axis='x', tight=True)
plt.savefig('P03891-Fig4.svg', format='svg')
plt.show()

# Calculate the actual degree values for each residue across all edge files
all_degree_values = pd.DataFrame(all_degrees).fillna(0)

# Display the first few rows to verify
print(all_degree_values.head())

# Calculate mean and standard deviation of degree values
mean_degree = all_degree_values.mean(axis=0)
std_degree = all_degree_values.std(axis=0)

# Plot the actual degree values with standard deviation shaded
plt.figure(figsize=(48, 12))
plt.plot(formatted_labels, mean_degree.values, label='Mean Degree', color='b')
plt.fill_between(formatted_labels, mean_degree - std_degree, mean_degree + std_degree, color='b', alpha=0.2, label='Standard Deviation')
plt.xticks(rotation=90, fontsize=6)  # Set x-axis labels to vertical
plt.xlabel('Residue Position')
plt.ylabel('Degree')
plt.title('Mean degree value per pesidue position')
plt.legend()
plt.autoscale(enable=True, axis='x', tight=True)
plt.savefig('P03891-Fig5.svg', format='svg')
plt.show()

# Calculate assortativity values for each graph
assortativity_values = []

for file in edge_files:
    # Load the edge data
    df = pd.read_csv(file, sep='\t')
    
    # Create a graph
    G = nx.from_pandas_edgelist(df, 'NodeId1', 'NodeId2')
    
    # Calculate assortativity
    assortativity = nx.degree_assortativity_coefficient(G)
    assortativity_values.append(assortativity)

# Plot the assortativity values
plt.figure(figsize=(20, 10))
plt.plot(range(1, len(assortativity_values) + 1), assortativity_values, marker='o')
plt.xlabel('Graph Index')
plt.ylabel('Assortativity Coefficient')
plt.title('Assortativity Coefficient for Each Graph')
plt.savefig('P03891-Fig6.svg', format='svg')
plt.show()

# List all edge files in the current directory
edge_files = [file for file in os.listdir() if file.endswith('_edges.txt')]

# Re-import necessary libraries
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import os

# List all edge files in the current directory
edge_files = [file for file in os.listdir() if file.endswith('_edges.txt')]

# Function to calculate degree values
all_degrees = []

for file in edge_files:
    # Load the edge data
    df = pd.read_csv(file, sep='\t')
    
    # Create a graph
    G = nx.from_pandas_edgelist(df, 'NodeId1', 'NodeId2')
    
    # Calculate degree
    degree = dict(G.degree())
    all_degrees.append(degree)

# Convert degree values to a DataFrame
all_degree_values = pd.DataFrame(all_degrees).fillna(0)

# Calculate the standard deviation of the degree values
std_degree = all_degree_values.std(axis=0)

# Calculate the negative standard deviation values
neg_std_degree = -std_degree

# Display the first few negative standard deviation values to verify
print(neg_std_degree.head())

# Extend the x-axis by increasing the figure size and set x-axis labels to vertical
plt.figure(figsize=(48, 12))
plt.plot(formatted_labels, std_degree.values, color='gray')
plt.plot(formatted_labels, neg_std_degree.values, color='gray', linestyle='--')
plt.fill_between(formatted_labels, -1.0, 1.0, color='lightgrey', alpha=0.3)
plt.fill_between(formatted_labels, -0.5, 0.5, color='lightgrey', alpha=0.5)
plt.axhline(0, color='black', linewidth=0.8)
plt.ylim(-max(std_degree.values), max(std_degree.values))  # Center zero on y-axis
plt.xticks(rotation=90, fontsize=6)  # Set x-axis labels to vertical
plt.xlabel('Residue Position')
plt.ylabel('Relative Standard Deviation of Degree')
plt.title('Dsd trend')
plt.autoscale(enable=True, axis='x', tight=True)
plt.savefig('P03891-Fig7.svg', format='svg')
plt.show()