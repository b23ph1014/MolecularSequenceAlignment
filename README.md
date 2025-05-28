# MolecularSequenceAlignment
This project implements three fundamental algorithms for sequence alignment in computational biology: Needleman-Wunsch, Smith-Waterman, and Hirschberg’s Algorithm. These algorithms are used to align biological sequences (like DNA, RNA, or proteins) for similarity analysis.

# Features
Global Alignment using Needleman-Wunsch,Local Alignment using Smith-Waterman,Space-Efficient Global Alignment using Hirschberg’s algorithm,Runtime and memory usage tracking,Visual representation of alignments,Customizable scoring parameters: match, mismatch, and gap penalties

# Files
Molecular_sequence_alignment.cpp – Main implementation file
Sequence input files – Plain text files containing the sequences (one per line)

# Input Format
Each sequence should be in a separate text file (e.g., seq1.txt, seq2.txt), with one line containing the sequence (e.g., ACGTAG).

# Sample Output
```
Needleman-Wunsch Alignment (Global)
Alignment Score: 5

Visual Alignment:
Sequence 1: ACGTAG
            ..*.. 
Sequence 2: ACGGAG

Runtime: 0.00123 seconds
Memory Usage: 42 KB
```
# Algorithms Used
Needleman-Wunsch,Global,space complexity of O(m×n),Best for full-length alignments
Smith-Waterman,Local,space complexity of O(m×n),Best for finding matching subsequences
Hirschberg,Global,space complexity of O(n),Best for long sequences (low memory)
