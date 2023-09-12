from collections import defaultdict
import itertools
from Bio.Data import CodonTable
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import sys
import glob
import json
import pandas as pd

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

# Initialize a list to store data for each filename
data_list = []

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

# Initialize an empty dictionary to store frequencies from different files
all_frequencies = {}

# Find all files matching the pattern in the current directory
frequency_files = glob.glob('Biases/*.txt')

# Loop through each frequency file
for filename in frequency_files:
    with open(filename, 'r') as file:
        # Read the entire file as JSON-like data
        frequency_data = json.load(file)

        # Extract frequencies from the nested dictionaries
        frequencies = {}
        for base1, data in frequency_data.items():
            if base1 != 'k':
                frequencies[base1] = {}
                for base2, value in data.items():
                    frequencies[base1][base2] = float(value)

    # Initialize mutation likelihoods
    mutation_likelihoods = initialize_mutation_likelihoods(frequencies)

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
    cmap = plt.get_cmap("YlGnBu")
    cmap.set_bad(color='lightgrey')

    # Plot the amino acid heatmap
    plt.figure(figsize=(8, 8))
    plt.rcParams.update({'font.size': 20})
    ax = sns.heatmap(normalized_matrix, annot=False, cbar=False, cmap=cmap, mask=mask, xticklabels=aa_labels, yticklabels=aa_labels, vmax=1)

    # Add lines to separate groups
    for boundary in group_boundaries[1:-1]:
        plt.axhline(y=boundary, color='black', linewidth=0.5)
        plt.axvline(x=boundary, color='black', linewidth=0.5)

    # Color-code the x-axis and y-axis labels
    for i, (label, color) in enumerate(zip(ax.get_xticklabels(), aa_colors)):
        label.set_color(color)
        ax.get_yticklabels()[i].set_color(color)

    # Get the name of the frequency file
    newfilename = os.path.splitext(os.path.basename(filename))[0]

    # Create a directory for saving the SVG files if it doesn't exist
    directory = "Results_AAmatrices"
    os.makedirs(directory, exist_ok=True)

    # Save the heatmap as an SVG file with a name corresponding to the frequency file
    amino_acid_heatmap_filename = f"{os.path.splitext(os.path.basename(filename))[0]}_amino_acid_heatmap.svg"
    amino_acid_heatmap_full_path = os.path.join(directory, amino_acid_heatmap_filename)
    plt.savefig(amino_acid_heatmap_full_path, format="svg", dpi=600, transparent=True)
    plt.close()

    #Save csvs with probabilities
    flat_matrix = normalized_matrix.flatten()
    # Create a mask for the values that are greater than 0.05
    values_above_threshold_mask = flat_matrix >= 0.05
    # Use the mask to get the values that are greater than or equal to 0.05
    values_above_threshold = flat_matrix[values_above_threshold_mask]
    directory = "Parameters_AAmatrices"
    os.makedirs(directory, exist_ok=True)
    filename= f"{newfilename}_probs05.csv"
    full_path = os.path.join(directory, filename)
    np.savetxt(full_path, values_above_threshold, delimiter=",")
    filename= f"{newfilename}_probsall.csv"
    full_path = os.path.join(directory, filename)
    np.savetxt(full_path, flat_matrix, delimiter=",")
    filename= f"{newfilename}_normmatrix.csv"
    full_path = os.path.join(directory, filename)
    np.savetxt(full_path, normalized_matrix, delimiter=",")

    # Define the categories for amino acids
    amino_acid_categories = {
        "P": ["S", "T", "Y", "N", "Q", "C"], #polar
        "NP": ["A", "V", "I", "L", "M", "F", "W", "G", "P"], #non-polar
        "B": ["R", "H", "K"], #positively charged
        "A": ["D", "E"], #negatively charged
        "*": ["*"] #stopcodon
    }

    # Create a dictionary to map amino acids to their categories
    amino_acid_to_category = {}
    for category, acids in amino_acid_categories.items():
        for acid in acids:
            amino_acid_to_category[acid] = category

    # Get a list of category labels
    category_labels = [category for category in amino_acid_categories.keys() if category != "*"]

    threshold = 0.05
    probable_transitions_counts = np.zeros((len(category_labels), len(category_labels)), dtype=int)

    for i, category1 in enumerate(category_labels):
        for j, category2 in enumerate(category_labels):
            count = 0
            for aa1, codons1 in codon_table.items():
                if amino_acid_to_category[aa1] == category1:
                    for codon1 in codons1:
                        for aa2, codons2 in codon_table.items():
                            if amino_acid_to_category[aa2] == category2:
                                for codon2 in codons2:
                                    differences = sum(1 for a, b in zip(codon1, codon2) if a != b)
                                    if differences == 1:
                                        aa_idx1 = aa_labels.index(aa1)
                                        aa_idx2 = aa_labels.index(aa2)
                                        if normalized_matrix[aa_idx1, aa_idx2] > threshold:
                                            count += 1
            probable_transitions_counts[i, j] = count

    # Plot the category heatmap
    plt.figure(figsize=(2, 2))
    plt.rcParams.update({'font.size': 15})
    cmap = plt.get_cmap("Blues")

    # Calculate the appropriate cell size for the heatmap
    ax = sns.heatmap(probable_transitions_counts, annot=True, fmt="d", cmap=cmap, cbar=False,
                    xticklabels=category_labels, yticklabels=category_labels, annot_kws={"size": 15},
                    square=True, linewidths=0.5, linecolor='black', cbar_kws={"shrink": 0.7})

    # Set the cell size by adjusting the position of the annotations
    ax.set_yticklabels(ax.get_yticklabels(), rotation=45)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    # Create a directory for saving the SVG files if it doesn't exist
    directory = "Results_AAcategories"
    os.makedirs(directory, exist_ok=True)

    # Save the heatmap as an SVG file with a name corresponding to the frequency file
    category_heatmap_filename = f"{os.path.splitext(os.path.basename(filename))[0]}_category_heatmap.svg"
    category_heatmap_full_path = os.path.join(directory, category_heatmap_filename)
    plt.savefig(category_heatmap_full_path, format="svg", dpi=600, transparent=True)
    plt.close()

    #Save csvs with probabilities of aa categories
    directory = "Parameters_AAcategories"
    os.makedirs(directory, exist_ok=True)
    filename= f"{newfilename}_probsall.csv"
    full_path = os.path.join(directory, filename)
    np.savetxt(full_path, probable_transitions_counts, delimiter=",")


    # Extract the filename without "_Bias_probsall" part
    filename_without_bias = os.path.splitext(os.path.basename(filename))[0].replace("_Bias_probsall", "")

    # Create a dictionary to store the data for this filename
    file_data = {
        'Filename': filename_without_bias,
    }

        
    # Get the labels for rows and columns from the heatmap
    row_labels = category_labels
    column_labels = category_labels

    # Add the values from probable_transitions_counts to the dictionary
    for i in range(len(row_labels)):
        for j in range(len(column_labels)):
            field_name = f'{row_labels[i]} to {column_labels[j]}'
            value = probable_transitions_counts[i, j]
            file_data[field_name] = value

    # # Flatten the probable_transitions_counts array and add its values to the dictionary
    # flattened_counts = probable_transitions_counts.flatten()
    # for i, value in enumerate(flattened_counts):
    #     field_name = f'Field_{i + 1}'
    #     file_data[field_name] = value

    # Append the data for this filename to the data list
    data_list.append(file_data)
    
# Create a DataFrame from the list of dictionaries
df = pd.DataFrame(data_list)

# Define the output CSV file path
output_csv_path = 'stats.csv'

# Save the DataFrame as a CSV file
df.to_csv(output_csv_path, index=False)

print(f'Saved data to {output_csv_path}')

