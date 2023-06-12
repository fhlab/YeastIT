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
    # Calculate the closeness centrality for each node
    closeness_centralities = nx.closeness_centrality(G, distance='weight')

    # Normalize the closeness centrality to 1 per amino acid
    aa_closeness_sum = {aa: sum(closeness_centralities[codon] for codon in codons) for aa, codons in codon_table.items()}
    closeness_centralities_normalized = {codon: closeness_centralities[codon] / aa_closeness_sum[codon_table_inverse[codon]] for codon in codon_labels}

    # Codon data
    ecoli_usage = {
        "TTT": 0.58, "TTC": 0.42, "TTA": 0.14, "TTG": 0.13, 
        "TCT": 0.17, "TCC": 0.15, "TCA": 0.14, "TCG": 0.14, 
        "TAT": 0.59, "TAC": 0.41, "TAA": 0.61, "TAG": 0.09,
        "TGT": 0.46, "TGC": 0.54, "TGA": 0.30, "TGG": 1.00, 
        "CTT": 0.12, "CTC": 0.10, "CTA": 0.04, "CTG": 0.47, 
        "CCT": 0.18, "CCC": 0.13, "CCA": 0.20, "CCG": 0.49,
        "CAT": 0.57, "CAC": 0.43, "CAA": 0.34, "CAG": 0.66, 
        "CGT": 0.36, "CGC": 0.36, "CGA": 0.07, "CGG": 0.11, 
        "ATT": 0.49, "ATC": 0.39, "ATA": 0.11, "ATG": 1.00,
        "ACT": 0.19, "ACC": 0.40, "ACA": 0.17, "ACG": 0.25, 
        "AAT": 0.49, "AAC": 0.51, "AAA": 0.74, "AAG": 0.26,
        "AGT": 0.16, "AGC": 0.25, "AGA": 0.07, "AGG": 0.04, 
        "GTT": 0.28, "GTC": 0.20, "GTA": 0.17, "GTG": 0.35, 
        "GCT": 0.18, "GCC": 0.26, "GCA": 0.23, "GCG": 0.33, 
        "GAT": 0.63, "GAC": 0.37, "GAA": 0.68, "GAG": 0.32,
        "GGT": 0.35, "GGC": 0.37, "GGA": 0.13, "GGG": 0.15
    }

    # Create DataFrame
    df = pd.DataFrame({
        "Codon": list(closeness_centralities_normalized.keys()),
        "Closeness Centrality": list(closeness_centralities_normalized.values()),
        "E. coli Usage": [ecoli_usage[codon] for codon in closeness_centralities_normalized.keys()]
    })

    # Create scatter plot
    plt.figure(figsize=(10,6))
    sns.scatterplot(data=df, x='E. coli Usage', y='Closeness Centrality', hue='Codon', legend=False)
    plt.title(title)
    # Save the plot in the 'Properties' folder with the appropriate .tif name
    plt.savefig(f'Properties/{title}_EColi_Bias.tif', format='tif', dpi=300)







