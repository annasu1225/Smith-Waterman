#!/usr/bin/python
__author__ = "Anna Su"
__email__ = "anna.su@yale.edu"
__copyright__ = "Copyright 2023"
__license__ = "GPL"
__version__ = "1.0.0"

'''
Usage: python SW.py -i <input file> -s <score file>
Example: python SW.py -i input.txt -s blosum62.txt
Note: Smith-Waterman Algorithm
Acknowledgement: received advice and inspirations from Keyi Li and Jason Liu
References: 
    (1) https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf 
    (2) https://github.com/slavianap/Smith-Waterman-Algorithm/blob/master/Script.py
'''

import numpy as np
import pandas as pd
import argparse

### This is one way to read in arguments in Python.
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

### Set the variables for pointers
STOP = 0
LEFT = 1
UP = 2
DIAGONAL = 3


### Get the sequences from the input file
def readInput(input):
    # Open the file containing the sequences
    with open(input, 'r') as file:
        # Read the lines of the file
        sequences = [line.strip() for line in file]
    return sequences


### Get the score matrix and store it as a dictionary
def readMatrix(score):
    # Open the file containing the score matrix
    with open(score) as f:
        # Read the lines of the file into a list
        lines = f.readlines()
    # Create an empty dictionary to hold the score matrix
    scores = {}
    # Extract the column labels from the first line of the file
    cols = lines[0].split()
    # Loop through the remaining lines of the file
    for line in lines[1:len(lines)-1]:
        # Split the line into a list of strings
        items = line.split()
        # Extract the row label from the first item
        row = items[0][0]
        # Create an empty list to hold the scores for this row
        row_scores = []
        # Loop through the remaining items in the line
        for item in items[1:]:
            # Convert the item to an integer and append it to the row scores list
            row_scores.append(int(item))
        # Add the row to the score matrix dictionary
        scores[row] = dict(zip(cols, row_scores))
    return scores


### Get the output score matrices
def outputScores(seq1, seq2, openGap, extGap, scores):
    # Create four matrices
    # X records the best alignment score of seq1 and seq2 ending with a gap in seq2
    # Y records the best alignment score of seq1 and seq2 ending with a gap in seq1
    # M records the best alignment score of seq1 and seq2 ending with a character match or mismatch
    # P records the pointers
    M = pd.DataFrame(0, index=range(len(seq2) + 1), columns=range(len(seq1) + 1))
    X = pd.DataFrame(0, index=range(len(seq2) + 1), columns=range(len(seq1) + 1))
    Y = pd.DataFrame(0, index=range(len(seq2) + 1), columns=range(len(seq1) + 1))
    P = pd.DataFrame(0, index=range(len(seq2) + 1), columns=range(len(seq1) + 1))

    # Loop through the sequences to get the scores
    # Notes: row: j, seq2; col: i, seq1
    for i in range(1, len(seq1)+1):
        for j in range(1, len(seq2)+1):
            X.iloc[j][i] = max(openGap + M.iloc[j][i - 1], extGap + X.iloc[j][i - 1])
            Y.iloc[j][i] = max(openGap + M.iloc[j - 1][i], extGap + Y.iloc[j - 1][i])
            M.iloc[j][i] = max(scores[seq1[i - 1]][seq2[j - 1]] + M.iloc[j - 1][i - 1], X.iloc[j][i], Y.iloc[j][i], 0)

            # Track where the cell's value is coming from
            if M.iloc[j][i] == 0:
                P.iloc[j][i] = STOP
            elif M.iloc[j][i] == Y.iloc[j][i]:
                P.iloc[j][i] = LEFT
            elif M.iloc[j][i] == X.iloc[j][i]:
                P.iloc[j][i] = UP
            elif M.iloc[j][i] == scores[seq1[i - 1]][seq2[j - 1]] + M.iloc[j - 1][i - 1]:
                P.iloc[j][i] = DIAGONAL

    return X, Y, M, P


### Get alignment from traceback
def traceback(M, P, seq1, seq2):
    # Initialize variables
    aligned_seq1 = ")"
    aligned_seq2 = ")"
    current_aligned_seq1 = ""
    current_aligned_seq2 = ""
    match = ""

    # Find the best alignment score and its position in the M dataframe
    max_value = max(M.max())
    max_j, max_i = np.where(M == max_value)
    max_j, max_i = int(max_j), int(max_i)

    # Post-alignment sequence
    post1 = ""
    post2 = ""
    if max_i < len(M.columns):
        post1 = seq1[max_i:]
    if max_j < len(M):
        post2 = seq2[max_j:]

    # Trace and compute the pathway with the local alignment
    while P.iloc[max_j][max_i] != STOP:
        if P.iloc[max_j][max_i] == DIAGONAL:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = seq2[max_j - 1]
            if current_aligned_seq1 == current_aligned_seq2:
                match += "|"
            else:
                match += " "
            max_i -= 1
            max_j -= 1
        elif P.iloc[max_j][max_i] == UP:
            current_aligned_seq1 = seq1[max_i - 1]
            current_aligned_seq2 = '-'
            match += " "
            max_i -= 1
        elif P.iloc[max_j][max_i] == LEFT:
            current_aligned_seq1 = '-'
            current_aligned_seq2 = seq2[max_j - 1]
            match += " "
            max_j -= 1

        aligned_seq1 += current_aligned_seq1
        aligned_seq2 += current_aligned_seq2

    # Add open bracket for the aligned sequences
    aligned_seq1 += "("
    aligned_seq2 += "("

    # Pre-alignment sequence
    pre1 = ""
    pre2 = ""
    if max_i > 0:
        pre1 = seq1[:max_i]
    if max_j > 0:
        pre2 = seq2[:max_j]

    # Add spacing for pre-alignment sequence and match symbol string
    if max_j > max_i:
        pre1 = " " * (max_j - max_i) + pre1
        match = " " * (max_j + 1) + match[::-1]
    elif max_j < max_i:
        pre2 = " " * (max_i - max_j) + pre2
        match = " " * (max_i + 1) + match[::-1]
    else:
        match = match[::-1]

    # Reverse the order of the aligned sequences and concatenate with pre and post alignment sequences
    aligned_seq1 = pre1 + aligned_seq1[::-1] + post1
    aligned_seq2 = pre2 + aligned_seq2[::-1] + post2

    return aligned_seq1, aligned_seq2, match


### Write the output file
def output_file(aligned_seq1, aligned_seq2 , match, M, seq1 ,seq2):
    with open('output.txt', 'w') as f:
        # Format Sequences section
        f.write("-----------\n|Sequences|\n-----------\n")
        f.write("sequence1\n" + seq1 + "\nsequence2\n" + seq2 + "\n")

        # Format Score Matrix section
        f.write("--------------\n|Score Matrix|\n--------------\n")
        f.write("\t".join("  " + seq1) + "\n")
        f.write("\t".join([" "] + [str(x) for x in list(M.iloc[0])]) + "\n")
        for i in range(len(seq2)):
            f.write("\t".join([seq2[i]] + [str(x) for x in list(M.iloc[i+1])]) + "\n")

        # Format Best Local Alignment section
        f.write("----------------------\n|Best Local Alignment|\n----------------------\n")
        f.write('Alignment Score:' + str(max(M.max())) + "\n")
        f.write("Alignment Results:\n" + aligned_seq1 + "\n" + match + "\n" + aligned_seq2 + "\n")


### Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
    # Implement the readInput function to get the sequences from the input file
    sequences = readInput(inputFile)
    seq1 = sequences[0]
    seq2 = sequences[1]

    # Implement the readMatrix function to get the score matrix
    scores = readMatrix(scoreFile)

    ### Implement the outputScores function to get the output score and pointer matrices
    X, Y, M, P = outputScores(seq1, seq2, openGap, extGap, scores)

    ### Implement the traceback function to get the aligned sequences and match symbol string
    aligned_seq1, aligned_seq2, match = traceback(M, P, seq1, seq2)

    ### Implement the output file function
    output_file(aligned_seq1, aligned_seq2, match, M, seq1, seq2)


### Run your Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)