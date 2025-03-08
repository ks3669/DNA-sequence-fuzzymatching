#include "GenomeTrie.h"
#include <cstdlib>
#include <ctime>
#include <chrono>

// Converts a DNA base to an integer index.
// Input: a character representing a DNA base (A, C, G, T).
// Returns: an integer index (0-3) corresponding to the DNA base.
int baseToIndex(char base) {
    switch (base) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 0; // Replace 'N' and other invalid characters with 'A'
    }
}

// Converts an integer index back to a DNA base.
// Input: an integer index (0-3).
// Returns: a character representing a DNA base (A, C, G, T).
char indexToBase(int index) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    return bases[index % 4];
}

// Mutates a given DNA base to a different random base.
// Input: a character representing a DNA base to mutate.
// Returns: a character representing the mutated DNA base.
char mutateBase(char base) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    int mutationIndex = rand() % 4;
    return base == bases[mutationIndex] ? bases[(mutationIndex + 1) % 4] : bases[mutationIndex];
}

// Constructor for TrieNode that initializes all children to nullptr and isEndOfWord to false.
// No input parameters.
// No return value (constructor).
TrieNode::TrieNode() : isEndOfWord(false) {
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        children[i] = nullptr;
    }
}

// Destructor for TrieNode that deallocates all the children nodes.
// No input parameters.
// No return value (destructor).
TrieNode::~TrieNode() {
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        delete children[i];
    }
}

// Constructor for GenomeTrie that initializes the trie with a new root node.
// No input parameters.
// No return value (constructor).
GenomeTrie::GenomeTrie() : root(new TrieNode()) {}

// Destructor for GenomeTrie that deallocates the root node.
// No input parameters.
// No return value (destructor).
GenomeTrie::~GenomeTrie() { delete root; }

// Inserts a sequence into the trie.
// Input: a string representing the DNA sequence to insert.
// No return value.
void GenomeTrie::insertSequence(const std::string &sequence) {
    TrieNode *crawler = root;
    for (char base : sequence) {
        int index = baseToIndex(base);
        if (!crawler->children[index]) {
            crawler->children[index] = new TrieNode();
        }
        crawler = crawler->children[index];
    }
    crawler->isEndOfWord = true;
}

// Searches for a sequence in the trie with up to one mismatch allowed.
// Input: a string representing the DNA sequence to search.
// Returns: true if the sequence is found with up to one mismatch, false otherwise.
bool GenomeTrie::searchSequence(const std::string &sequence) {
    return searchRecursive(root, sequence, 0, 0);
}

// Counts the total number of nodes in the trie.
// No input parameters.
// Returns: the total number of nodes in the trie.
int GenomeTrie::countTrieNodes() {
    return countNodes(root);
}

// Helper function to search the trie recursively, allowing for one mismatch.
// Input: a pointer to the current TrieNode, the sequence being searched, the current depth in the sequence, and the number of mismatches found so far.
// Returns: true if the sequence is found with up to one mismatch, false otherwise.
bool GenomeTrie::searchRecursive(TrieNode *node, const std::string &sequence, int depth, int mismatches) {
    if (!node) {
        return false;
    }
    if (depth == sequence.size()) {
        return node->isEndOfWord && mismatches <= 1;
    }

    int index = baseToIndex(sequence[depth]);
    bool found = false;

    // Recursive search allowing one mismatch
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        if (mismatches < 1 || i == index) {
            found = searchRecursive(node->children[i], sequence, depth + 1, mismatches + (i == index ? 0 : 1));
            if (found) break;
        }
    }
    return found;
}

// Helper function to count the number of nodes in the trie recursively.
// Input: a pointer to the TrieNode from which to start counting.
// Returns: the total number of nodes from the given node downward.
int GenomeTrie::countNodes(TrieNode *node) {
    if (!node) return 0;
    int count = 1; // Count this node
    for (int i = 0; i < ALPHABET_SIZE; i++) {
        count += countNodes(node->children[i]);
    }
    return count;
}

// Loads a genome segment from a file into a string, replacing any 'N' bases with 'A'.
// Input: the filename to read from, a string to store the genome segment, and the desired segment size.
// No return value.
void loadGenomeSegment(const std::string &filename, std::string &genomeSegment, int segmentSize) {
    std::ifstream file(filename);
    std::string fullGenome;
    if (!file.is_open()) {
        throw std::runtime_error("Could not open the file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() or line[0] == '>') continue;
        for (char &c : line) {
            if (c == 'N') c = 'A'; // Replace 'N' with 'A'
        }
        fullGenome += line;
    }
    file.close();

    if (fullGenome.size() > segmentSize) {
        int start = rand() % (fullGenome.size() - segmentSize);
        genomeSegment = fullGenome.substr(start, segmentSize);
    } else {
        genomeSegment = fullGenome; // Use the whole genome if smaller than requested segment size
    }
}

// Loads DNA sequences from a file into the trie.
// Input: the filename from which to read sequences and the trie into which to insert them.
// No return value.
void populateTrieWithQueries(const std::string &filename, GenomeTrie &trie) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open the query file: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() or line[0] == '>') continue; // Skip header lines
        trie.insertSequence(line);
    }
    file.close();
}

// Generates a specified number of mutated sequences from a genome segment and inserts them into the trie.
// Input: the segment from which to generate sequences, the number of sequences to generate, the trie into which to insert the sequences, and a boolean indicating whether to mutate the sequences.
// No return value.
void generateMutatedSequences(const std::string &segment, int count, GenomeTrie &trie, bool mutate) {
    for (int i = 0; i < count; i++) {
        int start = rand() % (segment.size() - SEQUENCE_LENGTH + 1);
        std::string sequence = segment.substr(start, SEQUENCE_LENGTH);
        if (mutate) {
            for (char &c : sequence) {
                if (rand() < RAND_MAX * MUTATION_RATE) {
                    c = mutateBase(c);
                }
            }
        }
        trie.insertSequence(sequence);
    }
}

// Counts the number of sequences in a genome segment that match sequences in the trie, allowing for up to one mismatch.
// Input: the segment to search and the trie containing the sequences to match against.
// Returns: the number of matching sequences found.
int countMatchesWithMismatches(const std::string &segment, GenomeTrie &trie) {
    int matches = 0;
    for (int i = 0; i <= segment.size() - SEQUENCE_LENGTH; i++) {
        if (trie.searchSequence(segment.substr(i, SEQUENCE_LENGTH))) {
            matches++;
        }
    }
    return matches;
}

// Main function to handle command-line arguments and execute functions based on the problem type.
// Input: command-line arguments specifying paths to data files, the problem type, and the number of sequences to generate or match.
// Returns: 0 if successful, 1 if there is an error with input arguments.
int main(int argc, char* argv[]) 
{
    if (argc != 5) 
	{
        std::cerr << "Usage: " << argv[0] << " <query dataset path> <subject dataset path> <problem type> <number of 36-mers>" << std::endl;
        return 1;
    }
    std::string queryPath = argv[1];
    std::string subjectPath = argv[2];
    std::string problemType = argv[3];
    int num36mers = atoi(argv[4]);
    std::srand(std::time(nullptr)); // Seed for random number generation
    GenomeTrie trie; //this contains the trie with characters taken from the subjectdataset
    std::string genomeSegment; //this contains the series of characters that we have selected as the segment
    if (problemType == "ProblemA" || problemType == "ProblemB") 
	{
        loadGenomeSegment(subjectPath, genomeSegment, 50000);  // Load a 50K segment from the subject dataset
        std::cout << "Loaded a 50K segment from the genome.\n";
        if (problemType == "ProblemA") //if the user chooses the ProblemA option we generate the nmers without using the mutation
		{
            generateMutatedSequences(genomeSegment, num36mers, trie, false);
            std::cout << "Trie size after insertion: " << trie.countTrieNodes() << std::endl;
            std::cout << "Matches with up to 1 mismatch: " << countMatchesWithMismatches(genomeSegment, trie) << std::endl;
        } 
		else if (problemType == "ProblemB") //if the user chooses the ProblemB option we generate the nmers with mutation error
		{
            generateMutatedSequences(genomeSegment, num36mers, trie, true);
            std::cout << "Trie size after insertion with mutations: " << trie.countTrieNodes() << std::endl;
            std::cout << "Matches with up to 1 mismatch: " << countMatchesWithMismatches(genomeSegment, trie) << std::endl;
        }
    } 
	else if (problemType == "ProblemC") 
	{
        populateTrieWithQueries(queryPath, trie); // Load all queries into the trie
        int segmentSizes[] = {100000, 10000000, 100000000}; // Sizes for ProblemC given in the question
        for (int i = 0; i < 3; i++) //we iterate through our segments sizes to conduct the perfect match search
		{
            loadGenomeSegment(subjectPath, genomeSegment, segmentSizes[i]);
            auto start = std::chrono::high_resolution_clock::now(); //start the clock to calculate the time it took for searching
            int hits = 0; //this variable holds the number of hits that we got
            for (int j = 0; j <= genomeSegment.length() - SEQUENCE_LENGTH; j++) 
			{
                std::string sequence = genomeSegment.substr(j, SEQUENCE_LENGTH);
                if (trie.searchSequence(sequence)) 
				{
                    hits++; //here we increment the hits when we find a hit
                }
            }
            auto end = std::chrono::high_resolution_clock::now(); //here we stop to clock to calculate the time taken
            std::chrono::duration<double> elapsed = end - start;
            std::cout << "Time taken to find all 36-mers in a " << segmentSizes[i] << " character segment: " << elapsed.count() << " seconds.\n"; //time taken to find
            std::cout << "Number of hits found for " << segmentSizes[i] << " segment: " << hits << std::endl; //number of hits we got
        }
    } 
	else 
	{
        std::cerr << "Unsupported problem specified." << std::endl; //display this when user enters a invalid problem type
        return 1;
    }
    return 0;
}
