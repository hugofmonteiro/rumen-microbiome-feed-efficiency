# This method using only DAA methods is fully explained in https://microbiome.github.io/OMA/
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import networkx as nx
from fa2 import ForceAtlas2
from networkx.linalg.graphmatrix import adjacency_matrix

# First: Upload the main taxonomy dataframe naming the variables
ALL_variableNames_MLRidge = pd.read_csv('ALL_variableNames_MLRidge.csv')
ALL_variableNames_MLRidge.head

# Second: Upload the best variables from ML
Ridge_BestVariables = pd.read_csv('Ridge_BestVariables.csv')

# Get the unique values of the "Variable" column in df2
variable_values = Ridge_BestVariables["Variable"].unique()

# Subset df1 based on the values of the "Variable" column in df2
Ridge_BestVariablesNames = ALL_variableNames_MLRidge[ALL_variableNames_MLRidge["Variable"].isin(variable_values)]

# Third: Upload best variables from LR

ALL_variableNames_LR = pd.read_csv('ALL_variableNames_LR.csv')
LR_BestVariables = pd.read_csv('LR_BestVariables.csv')

# Get the unique values of the "Variable" column in df2
variable_values_LR = LR_BestVariables["Variable"].unique()

# Fourth: Merge dataframes

# Subset df1 based on the values of the "Variable" column in df2
LR_BestVariablesNames = ALL_variableNames_LR[ALL_variableNames_LR["Variable"].isin(variable_values_LR)]
merged_df = pd.concat([Ridge_BestVariablesNames, LR_BestVariablesNames], axis=0, ignore_index=True)
print(merged_df)

# Drop the 'Variable' column
merged_df = merged_df.drop(columns=['Variable'])

# Prepare merged dataframe for network analysis

# Split the "Name" column based on the "_" character
split_names = merged_df["Name"].str.split('_', expand=True)

# Assign the split values to new columns in the DataFrame
merged_df['Name'] = split_names[0]
merged_df['Method'] = split_names[1]

# Fifth: Upload the best variables from ALDEx2, ANCOM-BC, MaAslin2, and LinDA
DAA_BestVariables = pd.read_csv('DAA_BestVariables.csv')

# Sixth: Merge the last dataset with previous ones. Do not merge the Fifth earlier
merged_complete = pd.concat([merged_df, DAA_BestVariables], axis=0, ignore_index=True)

# Network analysis on major microbial taxa contributing to a trait
# Add your renaming dictionary here
name_replacements = {
    'CountML': 'Count ML',
    'CountLR': 'Count LR',
    'relativeAbundanceML': 'RA ML',
    'relativeAbundanceLR': 'RA LR',
    'CLRML': 'CLR ML',
    'CLRLR': 'CLR LR',
    'aldex2': 'ALDEx2',
    'maaslin2': 'MaAslin2',
    'ancombc': 'ANCOM-BC',
    'LinDA': 'LinDA',
    # Add more name replacements if you desire to have more DAA tests
}

# Replace the old names with the new names in the DataFrame
merged_complete['Method'] = merged_complete['Method'].replace(name_replacements)

# Split the DataFrame based on the "Name" categories
Count_ML = merged_complete[merged_complete["Method"] == "Count ML"]
RA_ML = merged_complete[merged_complete["Method"] == "RA ML"]
CLR_ML = merged_complete[merged_complete["Method"] == "CLR ML"]
Count_LR = merged_complete[merged_complete["Method"] == "Count LR"]
RA_LR = merged_complete[merged_complete["Method"] == "RA LR"]
CLR_LR = merged_complete[merged_complete["Method"] == "CLR LR"]
ALDEx2 = merged_complete[merged_complete["Method"] == "ALDEx2"]
MaAslin2 = merged_complete[merged_complete["Method"] == "MaAslin2"]
ANCOMBC = merged_complete[merged_complete["Method"] == "ANCOM-BC"]
LinDA = merged_complete[merged_complete["Method"] == "LinDA"]

# Get the "Value" sets for each category
Count_ML_values = set(Count_ML["Name"])
RA_ML_values = set(RA_ML["Name"])
CLR_ML_values = set(CLR_ML["Name"])
Count_LR_values = set(Count_LR["Name"])
RA_LR_values = set(RA_LR["Name"])
CLR_LR_values = set(CLR_LR["Name"])
ALDEx2_values = set(ALDEx2["Name"])
MaAslin2_values = set(MaAslin2["Name"])
ANCOMBC_values = set(ANCOMBC["Name"])
LinDA_values = set(LinDA["Name"])

# Option 1: Take into account all 10 possibilities of differential abundance analysis (DAA) tests:

#!pip uninstall -y networkx
#!pip install networkx==3.0

#import networkx as nx
#print(nx.__version__)

# You can also skip directly until this point and replace this with your actual summary DataFrame
df = pd.DataFrame(merged_complete)

# Add your renaming dictionary here
name_replacements = {
    'CountML': 'Count ML',
    'CountLR': 'Count LR',
    'relativeAbundanceML': 'RA ML',
    'relativeAbundanceLR': 'RA LR',
    'CLRML': 'CLR ML',
    'CLRLR': 'CLR LR',
    'aldex2': 'ALDEx2',
    'maaslin2': 'MaAslin2',
    'ancombc': 'ANCOM-BC',
    'LinDA': 'LinDA',
    # Add more name replacements as needed
}

# Replace the old names with the new names in the DataFrame
merged_complete['Method'] = merged_complete['Method'].replace(name_replacements)

# Create a network graph object
G = nx.Graph()

# Add nodes and edges to the graph based on the DataFrame
for index, row in df.iterrows():
    G.add_node(row['Name'], node_type='name')
    G.add_node(row['Method'], node_type='method')
    G.add_edge(row['Name'], row['Method'])

# Set figure resolution and size
plt.figure(dpi=600, figsize=(12, 12))

# Initialize the ForceAtlas2 object
forceatlas2 = ForceAtlas2(
    # Set your desired parameters for ForceAtlas2 here
)

# Apply the ForceAtlas2 layout to the graph
M = adjacency_matrix(G).tolil().astype('f')  # Add this line to create the adjacency matrix M
positions = forceatlas2.forceatlas2(M, pos=None, iterations=50)
positions = np.array(positions)

# Convert the list of positions to a dictionary
pos = {node: (x, y) for node, x, y in zip(G.nodes(), positions[:, 0], positions[:, 1])}

# Optional: Create a spring layout with custom parameters (increase iterations, change k value)
#pos = nx.spring_layout(G, seed=42, iterations=400, k=0.10)
#pos = nx.circular_layout(G, scale=2)
#pos = nx.kamada_kawai_layout(G, scale=4)

# Define colors (you need to check this in your summary table to understand what was the maximum number of tests you found differences)
name_color = "#77B5A8"
method_color = "#c9e3a2"
degree_colors = {
    4: "#9ea616",
    5: "#686cb7",
    6: "#afd4f2"
}    #,
#    7: "#eac1dd",
#}

# Highlight nodes connected to four or more other nodes
highlighted_nodes = [node for node, degree in G.degree() if degree >= 4 and G.nodes[node]['node_type'] == 'name']

# Create a dictionary for node sizes
node_sizes_dict = {
    4: 350,
    5: 400,
    6: 450
}
    #    ,
#    7: 600
#}

# Set node colors based on degree and degree_colors
node_colors = []
for node, data in G.nodes(data=True):
    if data['node_type'] == 'name':
        degree = G.degree(node)
        if node in highlighted_nodes:
            node_colors.append(degree_colors[degree])
        else:
            node_colors.append(name_color)
    else:
        node_colors.append(method_color)

# Initialize node_sizes list
node_sizes = []

# Set node sizes based on degree and node_sizes_dict
for node, data in G.nodes(data=True):
    if data['node_type'] == 'name':
        degree = G.degree(node)
        if node in highlighted_nodes:
            node_sizes.append(node_sizes_dict[degree])
        else:
            node_sizes.append(75)
    else:
        node_sizes.append(800)

# Create a dictionary to map nodes to their index
highlighted_nodes_sorted = sorted(highlighted_nodes, key=lambda x: G.degree(x), reverse=True)
node_index_dict = {node: index for index, node in enumerate(highlighted_nodes_sorted, 1)}

# Draw the edges with transparency and lower zorder
for (n1, n2) in G.edges():
    x_values = [pos[n1][0], pos[n2][0]]
    y_values = [pos[n1][1], pos[n2][1]]
    plt.plot(x_values, y_values, color='k', alpha=0.3, zorder=0)

# Draw the nodes with the specified colors and sizes, and a higher zorder
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=1)

font_family = "Arial"

# Draw the labels for the nodes connected to four or more other nodes
highlighted_node_labels = {node: f"{node_index_dict[node]}" for node in highlighted_nodes}
nx.draw_networkx_labels(G, pos, labels=highlighted_node_labels, font_size=12, font_weight='light', font_family=font_family)

# Draw the method labels
node_labels = {}
for node, data in G.nodes(data=True):
    if data['node_type'] == 'method':
        node_labels[node] = name_replacements.get(node, node)

# Draw the method labels
nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=10, font_weight='bold', font_family=font_family)

# Create a sorted dictionary of connected nodes
connected_nodes_dict = {}
for node in highlighted_nodes:
    degree = G.degree(node)
    if degree not in connected_nodes_dict:
        connected_nodes_dict[degree] = []
    connected_nodes_dict[degree].append(node)

# Display the list of nodes with 4 or more connections on the side of the graph
connected_nodes_str = ""
for degree, nodes in sorted(connected_nodes_dict.items(), reverse=True):
    connected_nodes_str += f"{degree} Tests:\n"
    connected_nodes_str += "\n".join([f"{node_index_dict[node]}. {node}" for node in nodes])
    connected_nodes_str += "\n\n"

plt.text(1.0b2, 0.5, f"Legend title\nlegend subtitle\n\n{connected_nodes_str}", transform=plt.gca().transAxes, fontsize=10, verticalalignment='center')

# Add a custom legend outside the graph
legend_elements = [
    plt.Line2D([0], [0], marker='o', color=method_color, label='DAA Method', linestyle='', markersize=10)
]

# Add the legend elements for the degrees with custom colors
for degree in sorted(degree_colors, reverse=True):
    legend_elements.append(
        plt.Line2D([0], [0], marker='o', color=degree_colors[degree], label=f'{degree} Tests', linestyle='', markersize=10)
    )

legend_elements.append(
    plt.Line2D([0], [0], marker='o', color=name_color, label='≤ 3 Tests', linestyle='', markersize=10)
)
    
plt.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=14)

# Remove axis and display the graph
plt.axis('off')
plt.show()

# Extract names of nodes with more than 3 connections
nodes_with_more_than_3_connections = []

for node in highlighted_nodes:
    nodes_with_more_than_3_connections.append(node)

print("Nodes with more than 3 connections:")
for node_name in nodes_with_more_than_3_connections:
    print(node_name)
