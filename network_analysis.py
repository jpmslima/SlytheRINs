
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

## Re-plot the standard deviation values using sorted residue numbers
#plt.figure(figsize=(20, 10))
#plt.plot(sorted_residue_numbers, sorted_std_degree.values, color='gray')
#plt.plot(sorted_residue_numbers, neg_sorted_std_degree.values, color='gray', linestyle='--')
#plt.fill_between(sorted_residue_numbers, -1.0, 1.0, color='lightgrey', alpha=0.3)
#plt.fill_between(sorted_residue_numbers, -0.5, 0.5, color='lightgrey', alpha=0.5)
#plt.axhline(0, color='black', linewidth=0.8)
#plt.ylim(-max(sorted_std_degree.values), max(sorted_std_degree.values))  # Center zero on y-axis
#plt.xticks(rotation=90)  # Set x-axis labels to vertical
#plt.xlabel('Residue Position')
#plt.ylabel('Standard Deviation of Degree')
#plt.title('Standard Deviation of Degree with Centered Zero')
#plt.grid(True)
#plt.show()

## Load edge data from all relevant files for each interaction type
#interaction_types = ['HBOND', 'IONIC', 'PICATION', 'VDW', 'PIPISTACK', 'PIHBOND', 'SSBOND', 'METAL_ION', 'HALOGEN', 'IAC']
#
## Initialize a dictionary to store differences for each interaction type
#interaction_differences = {interaction: {} for interaction in interaction_types}
#
## Process each interaction type
#for interaction in interaction_types:
#    # Create a dictionary to store edges for each model
#    edge_dict = {}
#    
#    for file in edge_files:
#        model_number = file.split('_')[2]  # Extract model number from filename
#        df = pd.read_csv(file, sep='\t')
#        interaction_edges = df[df['Interaction'].str.contains(interaction)]
#        edge_dict[model_number] = set(zip(interaction_edges['NodeId1'], interaction_edges['NodeId2']))
#    
#    # Calculate differences in interactions for each node
#    for model1, edges1 in edge_dict.items():
#        for model2, edges2 in edge_dict.items():
#            if model1 != model2:
#                diff = edges1.symmetric_difference(edges2)
#                for edge in diff:
#                    node1, node2 = edge
#                    interaction_differences[interaction][node1] = interaction_differences[interaction].get(node1, 0) + 1
#                    interaction_differences[interaction][node2] = interaction_differences[interaction].get(node2, 0) + 1
#
## Convert differences to DataFrames for each interaction type
#interaction_diff_dfs = {interaction: pd.DataFrame(list(differences.items()), columns=['Node', 'DifferenceCount'])
#                        for interaction, differences in interaction_differences.items()}
#
## Display the first few rows of the HBOND differences DataFrame
#print(interaction_diff_dfs['HBOND'].head())
#
## Load interaction data for each type from the edge files
#interaction_data = {interaction: [] for interaction in interaction_types}
#
#for file in edge_files:
#    df = pd.read_csv(file, sep='\t')
#    for interaction in interaction_types:
#        interaction_edges = df[df['Interaction'].str.contains(interaction)]
#        interaction_data[interaction].append(interaction_edges)
#
## Concatenate data for each interaction type
#interaction_data = {interaction: pd.concat(data, ignore_index=True) for interaction, data in interaction_data.items()}
#
## Display the first few rows of the HBOND interaction data
#print(interaction_data['HBOND'].head())
#
## Load interaction data for each type from the edge files
#interaction_data = {interaction: [] for interaction in interaction_types}
#
#for file in edge_files:
#    df = pd.read_csv(file, sep='\t')
#    for interaction in interaction_types:
#        interaction_edges = df[df['Interaction'].str.contains(interaction)]
#        interaction_data[interaction].append(interaction_edges)
#
## Concatenate data for each interaction type
#interaction_data = {interaction: pd.concat(data, ignore_index=True) for interaction, data in interaction_data.items()}
#
## Display the first few rows of the HBOND interaction data
#print(interaction_data['HBOND'].head())
#
## Calculate mean, median, and standard deviation for each interaction type
#interaction_stats = {}
#
#for interaction, data in interaction_data.items():
#    # Calculate statistics
#    mean_values = data.groupby('NodeId1').size().mean()
#    median_values = data.groupby('NodeId1').size().median()
#    std_values = data.groupby('NodeId1').size().std()
#    
#    # Store the statistics
#    interaction_stats[interaction] = {
#        'mean': mean_values,
#        'median': median_values,
#        'std': std_values
#    }
#
## Display the calculated statistics for each interaction type
#print(interaction_stats)
#
## Create boxplots for each interaction type with residues ordered on the x-axis
#plt.figure(figsize=(80, 40))
#
## Iterate over each interaction type to create a boxplot
#for i, (interaction, data) in enumerate(interaction_data.items(), 1):
#    plt.subplot(2, 5, i)
#    sns.boxplot(x='NodeId1', y='Distance', data=data)
#    plt.title(f'Boxplot of {interaction} Interactions')
#    plt.xticks(rotation=90)
#    plt.xlabel('Residue')
#    plt.ylabel('Distance')
#    plt.grid(True)
#
#plt.tight_layout()
#plt.show()
#
## Redefine interaction types
#interaction_types = ['HBOND', 'IONIC', 'PICATION', 'VDW', 'PIPISTACK', 'PIHBOND', 'SSBOND', 'METAL_ION', 'HALOGEN', 'IAC']
#
## Reload interaction data for each type from the edge files
#interaction_data = {interaction: [] for interaction in interaction_types}
#
#for file in edge_files:
#    df = pd.read_csv(file, sep='\t')
#    for interaction in interaction_types:
#        interaction_edges = df[df['Interaction'].str.contains(interaction)]
#        interaction_data[interaction].append(interaction_edges)
#
## Concatenate data for each interaction type
#interaction_data = {interaction: pd.concat(data, ignore_index=True) for interaction, data in interaction_data.items()}
#
## Extract interaction subtypes from the interaction data
#for interaction, data in interaction_data.items():
#    # Add a new column for the subtype
#    data['Subtype'] = data['Interaction'].apply(lambda x: x.split(':')[1] if ':' in x else 'Unknown')
#    interaction_data[interaction] = data
#
## Display the first few rows of the HBOND interaction data with subtypes
#print(interaction_data['HBOND'].head())
#
## Create boxplots for each interaction subtype
#plt.figure(figsize=(40, 20))
#
## Iterate over each interaction type and its subtypes to create boxplots
#for i, (interaction, data) in enumerate(interaction_data.items(), 1):
#    plt.subplot(2, 5, i)
#    sns.boxplot(x='Subtype', y='Distance', data=data)
#    plt.title(f'Boxplot of {interaction} Interaction Subtypes')
#    plt.xticks(rotation=90)
#    plt.xlabel('Subtype')
#    plt.ylabel('Distance')
#    plt.grid(True)
#
#plt.tight_layout()
#plt.show()

## List of specified interactions
#specified_interactions = ['HBOND', 'IONIC', 'PICATION', 'VDW', 'PIPISTACK', 'PIHBOND', 'SSBOND', 'METAL_ION', 'HALOGEN', 'IAC']
#
## Initialize a dictionary to store data for specified interactions
#specified_interaction_data = {interaction: [] for interaction in specified_interactions}
#
## Load and filter interaction data
#for file in edge_files:
#    df = pd.read_csv(file, sep='\t')
#    for interaction in specified_interactions:
#        interaction_edges = df[df['Interaction'].str.contains(interaction, na=False)]
#        if not interaction_edges.empty:
#            specified_interaction_data[interaction].append(interaction_edges)
#
## Concatenate data for each specified interaction type
#specified_interaction_data = {interaction: pd.concat(data, ignore_index=True) for interaction, data in specified_interaction_data.items() if data}
#
## Display the first few rows of each interaction data if available
#if 'HBOND' in specified_interaction_data:
#    print(specified_interaction_data['HBOND'].head())
#else:
#    print('HBOND interaction not present in the data.')
#
#if 'IONIC' in specified_interaction_data:
#    print(specified_interaction_data['IONIC'].head())
#else:
#    print('IONIC interaction not present in the data.')
#
#if 'PICATION' in specified_interaction_data:
#    print(specified_interaction_data['PICATION'].head())
#else:
#    print('PICATION interaction not present in the data.')
#
#if 'VDW' in specified_interaction_data:
#    print(specified_interaction_data['VDW'].head())
#else:
#    print('VDW interaction not present in the data.')
#
#if 'PIPISTACK' in specified_interaction_data:
#    print(specified_interaction_data['PIPISTACK'].head())
#else:
#    print('PIPISTACK interaction not present in the data.')
#
#if 'PIHBOND' in specified_interaction_data:
#    print(specified_interaction_data['PIHBOND'].head())
#else:
#    print('PIHBOND interaction not present in the data.')
#
#if 'SSBOND' in specified_interaction_data:
#    print(specified_interaction_data['SSBOND'].head())
#else:
#    print('SSBOND interaction not present in the data.')
#
#if 'METAL_ION' in specified_interaction_data:
#    print(specified_interaction_data['METAL_ION'].head())
#else:
#    print('METAL_ION interaction not present in the data.')
#
#if 'HALOGEN' in specified_interaction_data:
#    print(specified_interaction_data['HALOGEN'].head())
#else:
#    print('HALOGEN interaction not present in the data.')
#
#if 'IAC' in specified_interaction_data:
#    print(specified_interaction_data['IAC'].head())
#else:
#    print('IAC interaction not present in the data.')
#
## Extract interaction data for each interaction type
#interaction_data_separated = {}
#
#for interaction in specified_interactions:
#    if interaction in specified_interaction_data:
#        interaction_data_separated[interaction] = specified_interaction_data[interaction]
#
## Display the first few rows of the HBOND interaction data if available
#if 'HBOND' in interaction_data_separated:
#    print(interaction_data_separated['HBOND'].head())
#else:
#    print('HBOND interaction not present in the data.')
#
## Calculate the count of each interaction type per residue and order residues
#interaction_counts_ordered = {}
#
#for interaction, data in interaction_data_separated.items():
#    # Group by NodeId1 and count the number of interactions, then sort by NodeId1
#    counts = data.groupby('NodeId1').size().reset_index(name='Count').sort_values(by='NodeId1')
#    interaction_counts_ordered[interaction] = counts
#
### Create individual boxplots for each interaction type
##for interaction, counts in interaction_counts_ordered.items():
##    plt.figure(figsize=(10, 6))
##    sns.boxplot(x='NodeId1', y='Count', data=counts)
##    plt.title(f'{interaction} Interaction Counts per Residue')
##    plt.xlabel('Residue')
##    plt.ylabel('Count')
##    plt.xticks(rotation=90)
##    plt.grid(axis='y')
##    plt.tight_layout()
##    plt.show()
#
### Create individual boxplots for each interaction type with outliers
##for interaction, counts in interaction_counts_ordered.items():
##    plt.figure(figsize=(10, 6))
##    sns.boxplot(x='NodeId1', y='Count', data=counts, showfliers=True)
##    plt.title(f'{interaction} Interaction Counts per Residue')
##    plt.xlabel('Residue')
##    plt.ylabel('Count')
##    plt.xticks(rotation=90)
##    plt.grid(axis='y')
##    plt.tight_layout()
##    plt.show()
##
## Calculate mean, median, and standard deviation for each interaction type and residue
#interaction_stats = {}
#
#for interaction, counts in interaction_counts_ordered.items():
#    stats = counts.groupby('NodeId1')['Count'].agg(['mean', 'median', 'std']).reset_index()
#    interaction_stats[interaction] = stats
#
### Create boxplots using the calculated statistics for each interaction type
##for interaction, stats in interaction_stats.items():
##    plt.figure(figsize=(10, 6))
##    sns.boxplot(x='NodeId1', y='mean', data=stats)
##    plt.title(f'{interaction} Interaction Mean Counts per Residue')
##    plt.xlabel('Residue')
##    plt.ylabel('Mean Count')
##    plt.xticks(rotation=90)
##    plt.grid(axis='y')
##    plt.tight_layout()
##    plt.show()
##
### Create boxplots using the calculated statistics for each interaction type
##for interaction, stats in interaction_stats.items():
##    plt.figure(figsize=(10, 6))
##    # Plot mean, median, and standard deviation as separate points
##    plt.errorbar(stats['NodeId1'], stats['mean'], yerr=stats['std'], fmt='o', label='Mean ± SD')
##    plt.scatter(stats['NodeId1'], stats['median'], color='red', label='Median')
##    plt.title(f'{interaction} Interaction Statistics per Residue')
##    plt.xlabel('Residue')
##    plt.ylabel('Count')
##    plt.xticks(rotation=90)
##    plt.grid(axis='y')
##    plt.legend()
##    plt.tight_layout()
##    plt.show()
#
## Prepare data for boxplot visualization using calculated statistics
#import numpy as np
#
## Create a DataFrame to hold the statistics for plotting
#boxplot_data = {}
#
#for interaction, stats in interaction_stats.items():
#    # Calculate the lower and upper bounds for the boxplot using mean and std
#    stats['lower'] = stats['mean'] - stats['std']
#    stats['upper'] = stats['mean'] + stats['std']
#    # Ensure no negative values
#    stats['lower'] = stats['lower'].clip(lower=0)
#    # Store the data for plotting
#    boxplot_data[interaction] = stats[['NodeId1', 'mean', 'median', 'lower', 'upper']]
#
## Display the prepared data for HBOND interactions
#if 'HBOND' in boxplot_data:
#    print(boxplot_data['HBOND'].head())
#else:
#    print('HBOND interaction not present in the data.')
#
## Generate boxplots for each interaction type using the calculated statistics
#for interaction, stats in boxplot_data.items():
#    plt.figure(figsize=(10, 6))
#    sns.boxplot(x='NodeId1', y='mean', data=stats)
#    plt.errorbar(stats['NodeId1'], stats['mean'], yerr=[stats['mean'] - stats['lower'], stats['upper'] - stats['mean']], fmt='o', label='Mean ± SD')
#    plt.scatter(stats['NodeId1'], stats['median'], color='red', label='Median')
#    plt.title(f'{interaction} Interaction Statistics per Residue')
#    plt.xlabel('Residue')
#    plt.ylabel('Mean Count')
#    plt.xticks(rotation=90)
#    plt.grid(axis='y')
#    plt.legend()
#    plt.tight_layout()
#    plt.show()
#
## Extract the number of interactions for each residue from the interaction data
#interaction_counts = {}
#
#for interaction, data in interaction_data_separated.items():
#    # Count the number of interactions per residue
#    counts = data['NodeId1'].value_counts().reset_index()
#    counts.columns = ['Residue', 'Count']
#    interaction_counts[interaction] = counts
#
## Order the residues in ascending order for each interaction type
#for interaction, counts in interaction_counts.items():
#    counts.sort_values(by='Residue', inplace=True)
#
## Create individual boxplots for each interaction type with ordered residues and outliers
#for interaction, counts in interaction_counts.items():
#    plt.figure(figsize=(10, 6))
#    sns.boxplot(x='Residue', y='Count', data=counts, showfliers=True)
#    plt.title(f'Boxplot of {interaction} Interaction Counts')
#    plt.xticks(rotation=90)
#    plt.xlabel('Residue')
#    plt.ylabel('Number of Interactions')
#    plt.grid(True)
#    plt.tight_layout()
#    plt.show()
