
import os
import pandas as pd

# Extract the pdbcode from the filenames in the directory
# If you have RINs from several PDB codes, the following lines will extract the PDB codes.
# Do not forget to verify if the PDB codes are correct.
pdbcode = [f[:4] for f in os.listdir() if f.endswith('_nodes.txt')]
print('RINs files from the following PDBs were found:', pdbcode)

# Extract unique values from  from both edges and nodes files
for code in pdbcode:
    data_nodes = pd.read_csv(code + "_nodes.txt", delimiter='\t')
    data_edges = pd.read_csv(code + "_edges.txt", delimiter='\t')
    unique_models_N = data_nodes['Model'].unique()
    unique_models_E = data_edges['Model'].unique()
    print(f"Unique models in node file for {code}: {unique_models_N}")
    print(f"Unique models in edges file for {code}: {unique_models_E}")

# Split the data based on the 'Model' column and save each subset into a separate .txt file
for model1 in unique_models_N:
    subset = data_nodes[data_nodes['Model'] == model1]
    filename = output_directory+code+'_model_' + str(model1) + '_nodes.txt'
    subset.to_csv(filename, sep='\t', index=False)
    print('Saved:', filename)
for model2 in unique_models_E:
    subset = data_edges[data_edges['Model'] == model2]
    filename = output_directory+code+'_model_' + str(model2) + '_edges.txt'
    subset.to_csv(filename, sep='\t', index=False)
    print('Saved:', filename)
# Renaming the original files
os.rename(code + "_edges.txt", input_directory + code + "_edges_original.txt.original")
