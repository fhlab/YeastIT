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

bias_files = glob.glob('Biases/*.txt') # Glob collects all file paths matching the pattern

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
    output_filename = f"{filename}_Transition_Probabilities_Heatmap_Dayhoff_categories.tif"
    output_filename_graph = f"{filename}_Amino_Acid_Conversion_Probabilities_Graph_Dayhoff_categories.tif"
    output_filename_directed_graph = f"{filename}_Amino_Acid_Conversion_Probabilities_Directed_Graph_Dayhoff_categories.tif"
    output_filename_interactive = f"{filename}_interactive.html"
    output_filename_pagerank = f"{filename}_PageRank.tif"



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

    aa_transitions = defaultdict(lambda: defaultdict(float))

    for aa1, codons1 in codon_table.items():
        for codon1 in codons1:
            for aa2, codons2 in codon_table.items():
                for codon2 in codons2:
                    differences = sum(1 for a, b in zip(codon1, codon2) if a != b)
                    if differences == 1:
                        nt1, nt2 = next((a, b) for a, b in zip(codon1, codon2) if a != b)
                        aa_transitions[aa1][aa2] += mutation_likelihoods[nt1][nt2]

    # Get a list of amino acid labels including the stop codon label
    aa_labels = list(aa_transitions.keys())

    # Create an empty NumPy array with the same dimensions as the number of amino acid labels
    transition_matrix = np.zeros((len(aa_labels), len(aa_labels)))

    # Assign the values from the aa_transitions dictionary to the transition_matrix array
    for i in range(len(aa_labels)):
        for j in range(len(aa_labels)):
            aa1 = aa_labels[i]
            aa2 = aa_labels[j]
            transition_matrix[i, j] = aa_transitions[aa1][aa2]


    # Set the diagonal values (self-transitions) to 0
    np.fill_diagonal(transition_matrix, 0)

    # Define the small constant to add to all the zero values
    small_constant = 0.0001

    # Create a boolean mask for the elements of the transition_matrix array that are zero
    zero_mask = (transition_matrix == 0)

    # Add the small constant to the elements of the transition_matrix array that are zero
    transition_matrix[zero_mask] += small_constant

    # Normalize the rows of the transition_matrix array so that they sum to 1
    row_sums = transition_matrix.sum(axis=1, keepdims=True)
    normalized_matrix = transition_matrix / row_sums

    # Set rows to 0 if all the values in the row are approximately equal
    tolerance = 1e-3
    for i in range(normalized_matrix.shape[0]):
        if np.allclose(normalized_matrix[i, :], normalized_matrix[i, 0], rtol=tolerance):
            normalized_matrix[i, :] = np.zeros(normalized_matrix.shape[1])

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
    group_boundaries = [0, 1, 6, 10, 13, 17, 20, 21]

    # Create a list to store the color of each amino acid
    aa_colors = []
    for i in range(len(group_boundaries) - 1):
        group_color = colors[i]
        for _ in range(group_boundaries[i + 1] - group_boundaries[i]):
            aa_colors.append(group_color)

    # Create a mask for the diagonal elements
    mask = np.zeros_like(normalized_matrix, dtype=bool)
    np.fill_diagonal(mask, True)

    # Create a custom colormap with grey color for the diagonal
    cmap = sns.color_palette("YlGnBu", as_cmap=True)
    cmap.set_bad(color='lightgrey')

    # Plot the heatmap
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(normalized_matrix, annot=False, cmap=cmap, mask=mask, xticklabels=aa_labels, yticklabels=aa_labels, vmax=1)

    # Add lines to separate groups
    for boundary in group_boundaries[1:-1]:
        plt.axhline(y=boundary, color='black', linewidth=0.5)
        plt.axvline(x=boundary, color='black', linewidth=0.5)

    # Color-code the x-axis and y-axis labels
    for i, (label, color) in enumerate(zip(ax.get_xticklabels(), aa_colors)):
        label.set_color(color)
        ax.get_yticklabels()[i].set_color(color)

    # Add a legend
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=key,
                              markerfacecolor=value, markersize=8)
                       for key, value in group_colors.items()]
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.17, 1))

    plt.xlabel("Amino Acid 2 (After)")
    plt.ylabel("Amino Acid 1 (Before)")

    # Get the name of the current script without the .py extension
    script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]

    # Create the Heatmap directory if it doesn't exist
    directory = "Heatmap"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Update the title and save paths
    plt.title(f"{title} Transition Probabilities Heatmap (excluding self-transitions) - Dayhoff categories")
    plt.savefig(os.path.join(directory, output_filename), format="tif", dpi=300)
    plt.close()

    import networkx as nx

    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes to the graph
    for aa in aa_labels:
        G.add_node(aa)

    # Add weighted edges to the graph
    for i in range(len(aa_labels)):
        for j in range(len(aa_labels)):
            aa1 = aa_labels[i]
            aa2 = aa_labels[j]
            weight = normalized_matrix[i, j]
            if weight > 0.01:
                G.add_edge(aa1, aa2, weight=weight)

    # Calculate the number of nodes
    node_count = len(G.nodes)

    # Calculate the total number of edges where weight > 0.01
    edge_count = sum(1 for u, v, w in G.edges(data='weight') if w >= 0.01)

    # Calculate the average number of edges per node
    avg_edges_per_node = edge_count / node_count

    print(f"Average number of edges with weight >= 0.01 per node: {avg_edges_per_node}")

    # Visualize the graph
    pos = nx.spring_layout(G, seed=42, k=k)  # Use the stored k value

    # Draw the nodes
    nx.draw_networkx_nodes(G, pos, node_color=aa_colors)

    # Draw the edges with varying transparency based on the weight
    edge_colors = [G[u][v]["weight"] for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, edge_cmap=plt.cm.Blues, alpha=0.5, width=1.5)

    # Draw node labels
    nx.draw_networkx_labels(G, pos, font_size=10)

    plt.axis("off")

    # Create the Graphs directory if it doesn't exist
    directory = "Graphs"
    if not os.path.exists(directory):
        os.makedirs(directory)

    plt.title(f"{title} Amino Acid Conversion Probabilities Graph + Dayhoff categories")
    plt.savefig(os.path.join("Graphs", output_filename_graph), format="tif", dpi=300)
    plt.close()


    from pyvis.network import Network

    # Initialize the pyvis Network
    net = Network(notebook=False, width="100%", height="100%")
    net.force_atlas_2based(gravity=-100, central_gravity=0.01, spring_length=200, spring_strength=0.08)

    # Add the nodes to the network
    for aa_label, color in zip(aa_labels, aa_colors):
        net.add_node(aa_label, color=color)

    # Add the edges to the network
    threshold = 0.01
    for i in range(len(aa_labels)):
        for j in range(len(aa_labels)):
            aa1 = aa_labels[i]
            aa2 = aa_labels[j]
            weight_forward = normalized_matrix[i, j]
            weight_reverse = normalized_matrix[j, i]

            if weight_forward > threshold or weight_reverse > threshold:
                title = f"{aa1} -> {aa2}: {weight_forward:.4f}, {aa2} -> {aa1}: {weight_reverse:.4f}"
                net.add_edge(aa1, aa2, value=(weight_forward + weight_reverse) / 2, title=title)

    # Create the Interactive_Graphs directory if it doesn't exist
    directory = "Interactive_Graphs"
    if not os.path.exists(directory):
        os.makedirs(directory)

    net.save_graph(os.path.join("Interactive_Graphs", output_filename_interactive))
    plt.close()

    ####### directed #######
    title = filename.replace("_Bias", "")

    # Create a directed graph
    G = nx.DiGraph()

    # Add nodes to the graph
    for aa in aa_labels:
        G.add_node(aa)

    # Add weighted edges to the graph
    for i in range(len(aa_labels)):
        for j in range(len(aa_labels)):
            aa1 = aa_labels[i]
            aa2 = aa_labels[j]
            weight = normalized_matrix[i, j]
            if weight > 0.01:
                G.add_edge(aa1, aa2, weight=weight)

    # Visualize the graph
    pos = nx.spring_layout(G, seed=42, k=k)  # Use the stored k value

    # Draw the nodes
    nx.draw_networkx_nodes(G, pos, node_color=aa_colors)

    # Draw the edges with varying transparency based on the weight
    edge_colors = [G[u][v]["weight"] for u, v in G.edges()]
    offset = 0.05

    for (u, v, weight) in G.edges(data='weight'):
        # Draw the edges with a small offset
        nx.draw_networkx_edges(G, pos, edgelist=[(u, v)], width=1.5, alpha=0.5,
                               edge_color=[weight], edge_cmap=plt.cm.Blues, 
                               edge_vmin=0, edge_vmax=1, ax=None,
                               arrows=True, arrowstyle='-|>', arrowsize=15,
                               connectionstyle=f"arc3, rad={offset}")

        # Draw the reverse edges with a small offset in the opposite direction
        if G.has_edge(v, u):
            nx.draw_networkx_edges(G, pos, edgelist=[(v, u)], width=1.5, alpha=0.5,
                                   edge_color=[G[v][u]["weight"]], edge_cmap=plt.cm.Blues,
                                   edge_vmin=0, edge_vmax=1, ax=None,
                                   arrows=True, arrowstyle='-|>', arrowsize=15,
                                   connectionstyle=f"arc3, rad={-offset}")

    # Draw node labels
    nx.draw_networkx_labels(G, pos, font_size=10)

    plt.axis("off")

    # Create the Directed_Graphs directory if it doesn't exist
    directory = "Directed_Graphs"
    if not os.path.exists(directory):
        os.makedirs(directory)

    plt.title(f"{title} Amino Acid Conversion Probabilities Directed Graph + Dayhoff categories")
    plt.savefig(os.path.join(directory, output_filename_directed_graph), format="tif", dpi=300)
    plt.close()


    # Flatten the matrix into a 1D array
    flat_matrix = normalized_matrix.flatten()

    # Create a mask for the values that are greater than 0.05
    values_above_threshold_mask = flat_matrix > 0.05

    # Use the mask to get the values that are greater than 0.05
    values_above_threshold = flat_matrix[values_above_threshold_mask]

    # Print the number of values that are greater than 0.05
    print(f"Number of values in the heatmap that are above 0.05: {len(values_above_threshold)}")

    average_connectivity = nx.average_node_connectivity(G)
    print(f"Average Node Connectivity: {average_connectivity}")

    # Normalize the weights
    norm_weights = nx.get_edge_attributes(G, "weight")
    norm_weights = {edge: weight / sum(norm_weights.values()) for edge, weight in norm_weights.items()}
    nx.set_edge_attributes(G, norm_weights, "weight")

    # PageRank
    pr = nx.pagerank(G, alpha=0.9)

    page_rank_dict = nx.pagerank(G, alpha=0.85)

    # We create lists of keys (nodes) and values (PageRank)
    nodes = list(page_rank_dict.keys())
    page_rank_values = list(page_rank_dict.values())

    # Then we create a bar chart
    plt.figure(figsize=(10,8))  # Set the figure size
    plt.bar(nodes, page_rank_values)  # Create the bar chart
    plt.title('PageRank Values')  # Set the title
    plt.xlabel('Nodes')  # Set the x-label
    plt.ylabel('PageRank')  # Set the y-label
    plt.xticks(rotation=90)  # Rotate x-axis labels if they overlap
    
    # Create the PageRank directory if it doesn't exist
    directory = "PageRank"
    if not os.path.exists(directory):
        os.makedirs(directory)

    plt.title(f"{title} PageRank values for each amino acid")
    plt.savefig(os.path.join(directory, output_filename_pagerank), format="tif", dpi=300)
    plt.close()
