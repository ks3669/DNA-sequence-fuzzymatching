# DNA-sequence-fuzzymatching
This Repo hosts my DNA sequence fuzzymatching's C++ code. This program uses a Trie to efficiently store and search DNA sequences, supports mutation analysis, and performs large-scale genomic queries. It handles billions of DNA sequences and executes fast searches with mismatch tolerance, making it useful for genomic research and mutation detection

I have given a quick rundown of the code below

1️⃣ Trie-Based DNA Sequence Storage
The program implements a Trie (prefix tree) structure to efficiently store and search DNA sequences from a dataset.
It uses a custom TrieNode class with four child nodes corresponding to A, C, G, and T bases.
The Trie enables fast insertion and lookup of genome sequences with potential mismatches.

2️⃣ DNA Sequence Manipulation
The program includes functions to convert, mutate, and handle DNA bases (A, C, G, T).
It replaces invalid characters (e.g., 'N') with 'A', ensuring clean genomic data processing.
Mutation functions randomly alter bases to simulate genetic variations for ProblemB.

3️⃣ Handling Large Genome Datasets
The program processes a confidential dataset containing 1 billion DNA sequences, each 32 characters long.
It extracts a genome segment of a specified size and loads it into memory for querying and mutation analysis.
It supports three modes: ProblemA (perfect match), ProblemB (mutated match), and ProblemC (large-scale query search).

4️⃣ Efficient Search with Mismatch Tolerance
The Trie structure supports searching sequences with up to one mismatch for DNA pattern recognition.
It uses recursive backtracking to traverse the Trie and find close matches.
This is useful for detecting genetic mutations and analyzing similar DNA sequences efficiently.

5️⃣ Mutated DNA Sequence Generation
For ProblemB, the program generates randomly mutated DNA sequences from the dataset.
It inserts these sequences into the Trie to analyze their similarity to the original dataset.
This helps simulate biological mutations and test the accuracy of DNA sequence matching.

6️⃣ Large-Scale Genome Query Matching
For ProblemC, the program performs a large-scale search for 36-mers in 100K, 10M, and 100M-character segments.
It times the execution speed for querying the dataset and reports the number of matching sequences.
This evaluates the scalability of the Trie structure in handling massive genomic data.

7️⃣ Command-Line Interface & Execution Control
The program is CLI-driven, taking input arguments for query path, subject dataset, problem type, and number of sequences.
It supports different execution modes (ProblemA, ProblemB, ProblemC) based on user selection.
Error handling is included to validate file paths, dataset formats, and problem type selection
