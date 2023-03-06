# Smith-Waterman

### Descriptions:
This package implements a Smith-Waterman local alignment algorithm. The default gap penalties are: opening gap -2, extension gap = -1.

### Input: 
1. A txt file containing two FASTA format sequences
2. A txt file containing blosum62 similarity scoring matrix

### Ouput:
1. A txt file containing the input sequences, score matrix, and best local alignment results. 

### Installation:
pip install "git+https://github.com/annasu1225/Smith-Waterman.git"

### Usage: 
python SW.py -i inputfile -s scorefile

### Example: 
python SW.py -i input.txt -s blosum62.txt

### Contact:
anna.su@yale.edu

### References:
1. [CMU CMSC423 Gap Penalties Lecture Slides](https://www.cs.cmu.edu/~ckingsf/bioinfo-lectures/gaps.pdf)
2. [Smith-Waterman Algorithm by Slavianap](https://github.com/slavianap/Smith-Waterman-Algorithm/blob/master/Script.py)
