from collections import defaultdict
import itertools
from Bio.Data import CodonTable
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import json
import glob
import networkx as nx
import collections
import pandas as pd

bias_files = glob.glob('Biases/*.txt') # Glob collects all file paths matching the pattern

results = []

for file_path in bias_files:
    print("Processing file:", file_path)  # Print the current file being processed

    filename = os.path.splitext(os.path.basename(file_path))[0]

    def initialize_mutation_likelihoods(frequencies):
        # Update the "C" dictionary with an "A" entry and the likelihood equal to the "T" value in the "G" dictionary
        frequencies["C"]["A"] = frequencies["G"]["T"]

        # Update the "G" dictionary with an "A" entry and the likelihood equal to the "T" value in the "C" dictionary
        frequencies["G"]["A"] = frequencies["C"]["T"]

        # Update the "G" dictionary with a "C" entry and the likelihood equal to the "G" value in the "C" dictionary
        frequencies["G"]["C"] = frequencies["C"]["G"]

        # Update the "T" dictionary with a "A" entry and the likelihood equal to the "T" value in the "A" dictionary
        frequencies["T"]["A"] = frequencies["A"]["T"]

        # Update the "T" dictionary with a "C" entry and the likelihood equal to the "G" value in the "A" dictionary
        frequencies["T"]["C"] = frequencies["A"]["G"]

        # Update the "T" dictionary with a "G" entry and the likelihood equal to the "C" value in the "A" dictionary
        frequencies["T"]["G"] = frequencies["A"]["C"]

        return frequencies

    # Define the 6 unique frequencies
    def read_frequencies_from_file(file_path):
        frequencies = {}
        k_value = None
        with open(file_path, 'r') as file:
            data = json.loads(file.read().strip())
            k_value = data.pop('k', None)
            frequencies = data
        return frequencies, k_value
    
    frequencies, k = read_frequencies_from_file(file_path)
    mutation_likelihoods = initialize_mutation_likelihoods(frequencies)
    
    # Get the title and output filenames from the bias filename
    title = filename.replace("_Bias", "")
    output_filename = f"{filename}_Transition_Probabilities_Heatmap.tif"
    output_filename_directed_graph = f"{filename}_Amino_Acid_Conversion_Probabilities_Directed_Graph.tif"

    codon_table = {
        "C": ["TGT", "TGC"],
        "A": ["GCT", "GCC", "GCA", "GCG"],
        "G": ["GGT", "GGC", "GGA", "GGG"],
        "P": ["CCT", "CCC", "CCA", "CCG"],
        "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
        "T": ["ACT", "ACC", "ACA", "ACG"],
        "D": ["GAT", "GAC"],
        "E": ["GAA", "GAG"],
        "N": ["AAT", "AAC"],
        "Q": ["CAA", "CAG"],
        "H": ["CAT", "CAC"],
        "K": ["AAA", "AAG"],
        "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
        "I": ["ATT", "ATC", "ATA"],
        "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
        "V": ["GTT", "GTC", "GTA", "GTG"],
        "M": ["ATG"],
        "F": ["TTT", "TTC"],
        "W": ["TGG"],
        "Y": ["TAT", "TAC"],
        "*": ["TAA", "TAG", "TGA"]
    }

    codon_transitions = defaultdict(lambda: defaultdict(float))

    for aa1, codons1 in codon_table.items():
        for codon1 in codons1:
            for aa2, codons2 in codon_table.items():
                for codon2 in codons2:
                    differences = sum(1 for a, b in zip(codon1, codon2) if a != b)
                    if differences == 1:
                        nt1, nt2 = next((a, b) for a, b in zip(codon1, codon2) if a != b)
                        codon_transitions[codon1][codon2] += mutation_likelihoods[nt1][nt2]


    # Use codon labels
    codon_labels = list(codon_transitions.keys())
    transition_matrix = np.zeros((len(codon_labels), len(codon_labels)))

    for i in range(len(codon_labels)):
        for j in range(len(codon_labels)):
            codon1 = codon_labels[i]
            codon2 = codon_labels[j]
            transition_matrix[i, j] = codon_transitions[codon1][codon2]

    # Set the diagonal values (self-transitions) to 0
    np.fill_diagonal(transition_matrix, 0)

    # Normalise by the sum of the matrix. This makes each value 'the likelihood of this codon change occurring relative to all others per unit of mutational time'
    total_sum = transition_matrix.sum()
    normalized_matrix = transition_matrix / total_sum

    # Initialize the colormap for the different groups
    colors = sns.color_palette("husl", 7)

    # Define the groups and their corresponding colors
    group_colors = {
        "Sulfur": colors[0],
        "Small": colors[1],
        "Acid/Amide": colors[2],
        "Basic": colors[3],
        "Hydrophobic": colors[4],
        "Aromaticity": colors[5],
        "Stop": colors[6],
    }

    # Define the group boundaries
    group_boundaries = [0, 2, 24, 32, 42, 56, 61, 64]

    # Create a list to store the color of each codon
    codon_colors = []
    for i in range(len(group_boundaries) - 1):
        group_color = colors[i]
        for _ in range(group_boundaries[i + 1] - group_boundaries[i]):
            codon_colors.append(group_color)

    # Create a mask for the diagonal elements
    mask = np.zeros_like(normalized_matrix, dtype=bool)
    np.fill_diagonal(mask, True)

    # Create a custom colormap with grey color for the diagonal
    cmap = sns.color_palette("YlGnBu", as_cmap=True)
    cmap.set_bad(color='lightgrey')

    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(normalized_matrix, annot=False, cmap=cmap, mask=mask, xticklabels=codon_labels, yticklabels=codon_labels, vmax=0.01)

    # Add lines to separate groups
    for boundary in group_boundaries[1:-1]:
        plt.axhline(y=boundary, color='black', linewidth=0.5)
        plt.axvline(x=boundary, color='black', linewidth=0.5)

    # Color-code the x-axis and y-axis labels
    for i, (label, color) in enumerate(zip(ax.get_xticklabels(), codon_colors)):
        label.set_color(color)
        ax.get_yticklabels()[i].set_color(color)

    # Add a legend
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=key,
                              markerfacecolor=value, markersize=8)
                       for key, value in group_colors.items()]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.17, 1))

    plt.xlabel("Codon 2 (After)")
    plt.ylabel("Codon 1 (Before)")

    # Get the name of the current script without the .py extension
    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]

    # Create the Heatmap directory if it doesn't exist
    directory = "Heatmap"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Update the title and save paths
    plt.title(f"{title} Transition Probabilities Heatmap")
    plt.savefig(os.path.join(directory, output_filename), format="tif", dpi=300)
    plt.close()


    ####### directed graph #######
    title = filename.replace("_Bias", "")

    # Create a directed graph
    G = nx.DiGraph()

    # Invert the codon table
    codon_table_inverse = {codon: aa for aa, codons in codon_table.items() for codon in codons}

    # Add nodes to the graph
    for codon in codon_labels:
        aa = codon_table_inverse[codon]
        G.add_node(codon, label = aa)

    # Add weighted edges to the graph
    for i in range(len(codon_labels)):
        for j in range(len(codon_labels)):
            codon1 = codon_labels[i]
            codon2 = codon_labels[j]
            weight = normalized_matrix[i, j]
            if weight > 0.00001:
                G.add_edge(codon1, codon2, weight=weight)

    # Create color palettes with distinct colors for each amino acid
    aa_palette = sns.color_palette("hls", len(codon_table))

    # Map amino acids to colors
    aa_to_color = dict(zip(codon_table.keys(), aa_palette))

    # Create a mapping of codons to colors and labels
    codon_to_color = {}
    for aa, codons in codon_table.items():
        for codon in codons:
            codon_to_color[codon] = aa_to_color[aa]

    # Use the mapping to get the color of each node
    node_colors = [codon_to_color[node] for node in G.nodes()]
   
    # Visualize the graph
    pos = nx.spring_layout(G, seed=1, k=k)  # Use the stored k value

    # Draw the nodes
    nx.draw_networkx_nodes(G, pos, node_color=node_colors)

    # Draw the node labels
    labels = nx.get_node_attributes(G, 'label')
    nx.draw_networkx_labels(G, pos, labels=labels)

    # Draw the edges with varying transparency based on the weight
    edge_colors = [G[u][v]["weight"] for u, v in G.edges()]
    offset = 0.05

    for (u, v, weight) in G.edges(data='weight'):
        # Draw the edges with a small offset
        nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=1.5, alpha=0.5,
                               edge_color=[weight], edge_cmap=plt.cm.Blues, 
                               edge_vmin=0, edge_vmax=0.01, ax=None,
                               arrows=True, arrowstyle='-|>', arrowsize=15,
                               connectionstyle=f"arc3, rad={offset}")

        # Draw the reverse edges with a small offset in the opposite direction
        if G.has_edge(v, u):
            nx.draw_networkx_edges(G, pos, edgelist=[(v, u)], width=1.5, alpha=0.5,
                                   edge_color=[G[v][u]["weight"]], edge_cmap=plt.cm.Blues,
                                   edge_vmin=0, edge_vmax=0.01, ax=None,
                                   arrows=True, arrowstyle='-|>', arrowsize=15,
                                   connectionstyle=f"arc3, rad={-offset}")
    

    plt.axis("off")

    # Create the Directed_Graphs directory if it doesn't exist
    directory = "Directed_Graphs"
    if not os.path.exists(directory):
        os.makedirs(directory)

    plt.title(f"{title} Codon Conversion Probabilities Directed Graph")
    plt.savefig(os.path.join(directory, output_filename_directed_graph), format="tif", dpi=300)
    plt.close()

    # Compute centrality measures
    degree_centrality = nx.degree_centrality(G)
    closeness_centrality = nx.closeness_centrality(G, distance = "weight")
    
    # Create a DataFrame to hold the centrality measures
    df = pd.DataFrame({
        'amino_acid': [codon_table_inverse[codon] for codon in G.nodes()],
        'degree': [degree_centrality[codon] for codon in G.nodes()],
        'closeness': [closeness_centrality[codon] for codon in G.nodes()]
    })
    
    # Aggregate by amino acid
    df_aggregated = df.groupby('amino_acid').mean()
    
    # Add the filename and the aggregated values to the results
    for aa in df_aggregated.index:
        results.append({
            'file': filename,
            'amino_acid': aa,
            'degree': df_aggregated.loc[aa, 'degree'],
            'closeness': df_aggregated.loc[aa, 'closeness']
        })

##### Plotting Closeness of aggregated amino acids #####

# Convert the results to a DataFrame
df_results = pd.DataFrame(results)

# Normalize the measures by the number of codons for each amino acid
df_results['degree'] /= df_results['amino_acid'].map(lambda aa: len(codon_table[aa]))
df_results['closeness'] /= df_results['amino_acid'].map(lambda aa: len(codon_table[aa]))

# Remove "_Bias" from the file column
df_results['file'] = df_results['file'].str.replace('_Bias', '')

# Group data by file and amino_acid, and compute mean closeness
grouped_data = df_results.groupby(['file', 'amino_acid'])['closeness'].mean().reset_index()

# Pivot data to have amino acids as columns, files as index and closeness as values
pivoted_data = grouped_data.pivot(index='file', columns='amino_acid', values='closeness')

# Calculate total closeness for each system and sort based on it
pivoted_data['total_closeness'] = pivoted_data.sum(axis=1)
pivoted_data.sort_values('total_closeness', inplace=True)

# Normalize the closeness centrality for each system so that it sums to 1
pivoted_data = pivoted_data.div(pivoted_data['total_closeness'], axis=0)

# Drop the total_closeness column before plotting
pivoted_data.drop(columns='total_closeness', inplace=True)

# Plot the data as stacked bar plot
ax = pivoted_data.plot(kind='bar', stacked=True, figsize=(10,7))

# Set the title and axis labels
plt.title('Average relative closeness centrality values of codons by amino acid')
plt.ylabel('Relative Closeness Centrality')

# Get current legend
handles, labels = ax.get_legend_handles_labels()

# Reverse the order of legend items and adjust the location of the legend
ax.legend(reversed(handles), reversed(labels), title='Amino Acid', bbox_to_anchor=(1.0, 1), loc='upper left')

# Create the Properties directory if it doesn't exist
directory = "Properties"
if not os.path.exists(directory):
    os.makedirs(directory)

plt.savefig("Properties/System_Closenesses.tif", format="tif", dpi=300)
plt.close()

