#ifndef GENOMETRIE_H  // Start of the preprocessor directive to prevent multiple inclusions of this header file.
#define GENOMETRIE_H  // Definition that marks the header has been included.
#include <string>     // Include the string library to use the std::string type.
#include <iostream>   // Include the iostream library for input and output operations.
#include <fstream>    // Include the fstream library to work with files.
#define ALPHABET_SIZE 4  // Define a constant for the number of possible bases in a DNA sequence.
#define SEQUENCE_LENGTH 36  // Define the length of the sequences used in the trie, here 36 bases long.
#define MUTATION_RATE 0.05  // Define the mutation rate as 5%.

int baseToIndex(char base);  // Convert a DNA base (A, C, G, T) to an index (0-3).
char indexToBase(int index);  // Convert an index (0-3) back to a DNA base (A, C, G, T).
char mutateBase(char base);  // Mutate a given DNA base to a different base according to a mutation rate.
struct TrieNode 
{  // Define the structure of each node in the trie.
    TrieNode* children[ALPHABET_SIZE];  // An array of pointers to child nodes, one for each base in the alphabet.
    bool isEndOfWord;  // Boolean flag to mark the end of a sequence in the trie.
    TrieNode();  // Constructor for TrieNode, initializes the node.
    ~TrieNode();  // Destructor for TrieNode, cleans up resources.
};

class GenomeTrie 
{  // Define a class for managing the trie that stores genomic sequences.
private:
    TrieNode* root;  // Pointer to the root node of the trie.

public:
    GenomeTrie();  // Constructor for GenomeTrie, initializes the trie.
    ~GenomeTrie();  // Destructor for GenomeTrie, cleans up resources.
    void insertSequence(const std::string& sequence);  // Insert a DNA sequence into the trie.
    bool searchSequence(const std::string& sequence);  // Search for a DNA sequence in the trie.
    int countTrieNodes();  // Count the total number of nodes in the trie.

private:
    bool searchRecursive(TrieNode* node, const std::string& sequence, int depth, int mismatches);  // Helper method to search with mismatches allowed.
    int countNodes(TrieNode* node);  // Helper method to recursively count nodes.
};

void loadGenomeSegment(const std::string& filename, std::string& genomeSegment, int segmentSize);  // Load a segment of a genome from a file into a string.
void populateTrieWithQueries(const std::string& filename, GenomeTrie& trie);  // Read sequences from a file and insert them into the trie.
void generateMutatedSequences(const std::string& segment, int count, GenomeTrie& trie, bool mutate);  // Generate and insert mutated sequences into the trie.
int countMatchesWithMismatches(const std::string& segment, GenomeTrie& trie);  // Count how many sequences in the trie match a given segment with possible mismatches.
#endif // GENOMETRIE_H  // End of the preprocessor directive.
