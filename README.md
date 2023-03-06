# Smith-Waterman

This package implements a Smith-Waterman local alignment algorithm. The default gap penalties are: opening gap -2, extension gap = -1.

Input: 
1. A txt file containing two FASTA format sequences
2. A txt file containing blosum62 similarity scoring matrix

Ouput:
1. A txt file containing the input sequences, score matrix, and best local alignment results. 

Usage: 
python SW.py -i <input file> -s <score file>

Example: 
python SW.py -i input.txt -s blosum62.txt
