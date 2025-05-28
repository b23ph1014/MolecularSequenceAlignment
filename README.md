# MolecularSequenceAlignment
This project implements three fundamental algorithms for sequence alignment in computational biology: Needleman-Wunsch, Smith-Waterman, and Hirschberg’s Algorithm. These algorithms are used to align biological sequences (like DNA, RNA, or proteins) for similarity analysis.

# Features
Global Alignment using Needleman-Wunsch,Local Alignment using Smith-Waterman,Space-Efficient Global Alignment using Hirschberg’s algorithm,Runtime and memory usage tracking,Visual representation of alignments,Customizable scoring parameters: match, mismatch, and gap penalties

# Files
Molecular_sequence_alignment.cpp – Main implementation file
Sequence input files – Plain text files containing the sequences (one per line)

⚙️ How to Compile and Run
🛠️ Compilation
Use g++ or any modern C++ compiler:

bash
Copy
Edit
g++ -std=c++11 -O2 Molecular_sequence_alignment.cpp -o align
🚀 Execution
bash
Copy
Edit
./align <sequence_file1> <sequence_file2> [match=1] [mismatch=-1] [gap=-2]
Example:

bash
Copy
Edit
./align seq1.txt seq2.txt 2 -1 -2
📥 Input Format
Each sequence should be in a separate text file (e.g., seq1.txt, seq2.txt), with one line containing the sequence (e.g., ACGTAG).

🖥️ Sample Output
text
Copy
Edit
Needleman-Wunsch Alignment (Global)
Alignment Score: 5

Visual Alignment:
Sequence 1: ACGTAG
            ..*.. 
Sequence 2: ACGGAG

Runtime: 0.00123 seconds
Memory Usage: 42 KB
📈 Algorithms Used
Algorithm	Type	Space Complexity	When to Use
Needleman-Wunsch	Global	O(m×n)	Best for full-length alignments
Smith-Waterman	Local	O(m×n)	Best for finding matching subsequences
Hirschberg	Global	O(n)	Best for long sequences (low memory)
