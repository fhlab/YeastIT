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

for file_path in bias_files:
    print("Processing file:", file_path)  # Print the current file being processed

    filename = os.path.splitext(os.path.basename(file_path))[0]

    def initialize_mutation_likelihoods(frequencies):
        frequencies["C"]["A"] = frequencies["G"]["T"]
        frequencies["G"]["A"] = frequencies["C"]["T"]
        frequencies["G"]["C"] = frequencies["C"]["G"]
        frequencies["T"]["A"] = frequencies["A"]["T"]
        frequencies["T"]["C"] = frequencies["A"]["G"]
        frequencies["T"]["G"] = frequencies["A"]["C"]

        return frequencies

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
    
    title = filename.replace("_Bias", "")

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

    codon_labels = list(codon_transitions.keys())
    transition_matrix = np.zeros((len(codon_labels), len(codon_labels)))

    for i in range(len(codon_labels)):
        for j in range(len(codon_labels)):
            codon1 = codon_labels[i]
            codon2 = codon_labels[j]
            transition_matrix[i, j] = codon_transitions[codon1][codon2]

    np.fill_diagonal(transition_matrix, 0)

    # Normalise by the sum of the matrix. This makes each value 'the likelihood of this codon change occurring relative to all others per unit of mutational time'
    total_sum = transition_matrix.sum()
    normalized_matrix = transition_matrix / total_sum


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

    # Compute the closeness centrality for each codon
    closeness_centrality = nx.closeness_centrality(G, distance='weight')

    # Assign closeness centrality to each codon within the corresponding amino acid
    closeness_dict = defaultdict(lambda: defaultdict(float))

    for codon, centrality in closeness_centrality.items():
        aa = codon_table_inverse[codon]
        closeness_dict[aa][codon] = centrality

    # Normalize the closeness centrality values
    for aa, codons in closeness_dict.items():
        total_centrality = sum(closeness_dict[aa].values())
        for codon in codons:
            closeness_dict[aa][codon] /= total_centrality

    # Prepare data for plotting
    plot_data = defaultdict(dict)

    for aa, codons in closeness_dict.items():
        for codon, centrality in codons.items():
            plot_data[aa][codon] = centrality

    df = pd.DataFrame(plot_data).T

    # Plot the data
    fig, ax = plt.subplots(figsize=(15, 10))
    df.plot(kind='bar', stacked=True, ax=ax)

    # Iterate over each amino acid (or bar in the bar plot)
    for i, amino_acid in enumerate(df.index):
        y_offset = 0  # Initialize the vertical offset to adjust the labels
        # Iterate over each codon in this amino acid
        for codon in df.columns:
            value = df.loc[amino_acid, codon]
            if pd.notnull(value):  # Check for non-missing values
                # Place the label in the middle of the segment
                ax.text(i, y_offset + value / 2, codon, ha='center', va='center', fontsize=8)
                y_offset += value  # Add the height of the bar to the offset for the next label
  
    ax.get_legend().remove()

    plt.title('Closeness Centrality per Codon')
    plt.xlabel('Amino Acid')
    plt.ylabel('Closeness Centrality')

    # Save the plot in the 'Properties' folder with the appropriate .tif name
    plt.savefig(f'Properties/{title}_closeness_centrality.tif', format='tif', dpi=300)

    # Calculate in-degree and out-degree for each node considering edge weights
    in_degree = {node: sum([G.edges[predecessor, node]['weight'] for predecessor in G.predecessors(node)]) for node in G.nodes}
    out_degree = {node: sum([G.edges[node, successor]['weight'] for successor in G.successors(node)]) for node in G.nodes}

    # Calculate the ratio of out-degree to in-degree (with a small constant to avoid division by zero)
    ratio_dict = defaultdict(lambda: defaultdict(float))
    for node in G.nodes:
        aa = codon_table_inverse[node]
        ratio_dict[aa][node] = out_degree[node] / (in_degree[node] + 1e-6)  # Small constant to avoid division by zero

    # Normalize these values
    for aa, codons in ratio_dict.items():
        total_ratio = sum(ratio_dict[aa].values())
        for codon in codons:
            ratio_dict[aa][codon] /= total_ratio

    # Prepare data for plotting
    plot_data = defaultdict(dict)

    for aa, codons in ratio_dict.items():
        for codon, ratio in codons.items():
            plot_data[aa][codon] = ratio

    df = pd.DataFrame(plot_data).T

    # Fill NA/NaN values with 0
    df = df.fillna(0)

    # Plot the data
    fig, ax = plt.subplots(figsize=(20, 10))
    df.plot(kind='bar', stacked=True, ax=ax)

    # Calculate the ratio of out-degree to in-degree for each codon (with a small constant to avoid division by zero)
    ratio_dict = {codon: out_degree[codon] / (in_degree[codon] + 1e-6) for codon in df.columns}

    # Iterate over each amino acid (or bar in the bar plot)
    for i, amino_acid in enumerate(df.index):
        y_offset = 0  # Initialize the vertical offset to adjust the labels
        # Iterate over each codon in this amino acid
        for codon in df.columns:
            value = df.loc[amino_acid, codon]
            if pd.notnull(value):  # Check for non-missing values
                ratio = ratio_dict[codon]
                # Place the label in the middle of the segment
                ax.text(i, y_offset + value / 2, codon, ha='center', va='center', fontsize=8)
                y_offset += value  # Add the height of the bar to the offset for the next label
  
    ax.get_legend().remove()

    plt.title(f'Out-degree to In-degree Ratio per Codon')
    plt.xlabel('Amino Acid')
    plt.ylabel('Out/In Degree Ratio')


    # Save the plot in the 'Properties' folder with the appropriate .tif name
    plt.savefig(f'Properties/{title}_degree_ratio.tif', format='tif', dpi=300)







