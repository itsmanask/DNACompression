#include <iostream>
#include <map>
#include <vector>
#include <iomanip> // for precision in arithmetic coding
#include <chrono>  // for time calculation
using namespace std;

struct node {
    string data1;
    int data2;
    bool used;
    int flag;
    node* left;
    node* right;

    node(string nucleotide, int frequency)
        : data1(nucleotide), data2(frequency), used(false), left(nullptr), right(nullptr), flag(-1) {}
};

class HuffmanBST {
public:
    map<char, int> freq_map;           // Frequency map for nucleotides
    map<char, string> huffmanCodes;    // Map to store the Huffman codes
    map<string, char> reverseHuffman;  // Reverse map for decoding Huffman sequences
    node* root = nullptr;              // Root of Huffman tree
    vector<node*> nodes;               // Vector to store pointers to nodes
    map<char, double> arithmeticProb;  // Probability map for arithmetic coding
    map<char, pair<double, double>> intervals; // Interval map for each nucleotide in arithmetic coding

    HuffmanBST() {
        nodes.push_back(new node("A", 0));
        nodes.push_back(new node("C", 0));
        nodes.push_back(new node("G", 0));
        nodes.push_back(new node("T", 0));
    }

    void countfreq(string dnasequence) {
        freq_map.clear();  // Ensure this is reset for every run
        for (char nucleotide : dnasequence) {
            freq_map[nucleotide]++;
        }

        for (auto& n : nodes) {
            n->data2 = freq_map[n->data1[0]];
        }
    }

    void displayNodes() {
        cout << "Nucleotide Frequencies:" << endl;
        for (auto& n : nodes) {
            cout << n->data1 << " = " << n->data2;
            if (n->used) {
                cout << " (combined)";
            }
            cout << endl;
        }
    }

    void findLowestFrequencyNucleotides(node*& lowest1, node*& lowest2) {
        lowest1 = nullptr;
        lowest2 = nullptr;
        for (auto& n : nodes) {
            if (!n->used) {
                if (lowest1 == nullptr || n->data2 < lowest1->data2) {
                    lowest2 = lowest1;
                    lowest1 = n;
                } else if (lowest2 == nullptr || n->data2 < lowest2->data2) {
                    lowest2 = n;
                }
            }
        }
    }

    void combineLowestNodes(node*& parent) {
        node* lowest1 = nullptr;
        node* lowest2 = nullptr;

        findLowestFrequencyNucleotides(lowest1, lowest2);

        if (lowest1 && lowest2) {
            int combinedFrequency = lowest1->data2 + lowest2->data2;
            string combinedNucleotide = lowest1->data1 + "+" + lowest2->data1;

            parent = new node(combinedNucleotide, combinedFrequency);
            parent->left = lowest1;
            parent->right = lowest2;
            lowest1->flag = 0;
            lowest2->flag = 1;

            lowest1->used = true;
            lowest2->used = true;

            cout << "\nCombined " << lowest1->data1 << " and " << lowest2->data1
                 << " into " << combinedNucleotide << " = " << combinedFrequency << endl;
        }
    }

    void buildHuffmanTree() {
        while (true) {
            node* lowest1 = nullptr;
            node* lowest2 = nullptr;
            findLowestFrequencyNucleotides(lowest1, lowest2);
            if (lowest1 == nullptr || lowest2 == nullptr) {
                break;
            }
            combineLowestNodes(root);
            nodes.push_back(root);
        }
    }

    void storeHuffmanCodes(node* root, string code) {
        if (root == nullptr) return;

        if (!root->left && !root->right) {
            char nucleotide = root->data1[0];
            huffmanCodes[nucleotide] = code;
            reverseHuffman[code] = nucleotide;  // Store reverse mapping for decoding
        }

        storeHuffmanCodes(root->left, code + "0");
        storeHuffmanCodes(root->right, code + "1");
    }

    void displayHuffmanCodes() {
        cout << "\nHuffman Codes for Nucleotides:" << endl;
        for (auto& pair : huffmanCodes) {
            cout << pair.first << " = " << pair.second << endl;
        }
    }

    string encodeDNA(string dnasequence) {
        string encodedDNA = "";
        for (char nucleotide : dnasequence) {
            encodedDNA += huffmanCodes[nucleotide];
        }
        return encodedDNA;
    }

    string decodeHuffman(string encodedDNA) {
        string decodedDNA = "";
        string code = "";
        cout << "\nHuffman Decoding Steps:\n";
        for (char bit : encodedDNA) {
            code += bit;
            if (reverseHuffman.find(code) != reverseHuffman.end()) {
                char decodedChar = reverseHuffman[code];
                decodedDNA += decodedChar;
                cout << "Decoded " << code << " to " << decodedChar << endl;
                code = "";
            }
        }
        return decodedDNA;
    }

    // Arithmetic Coding Methods

    void calculateProbabilities(string dnasequence) {
        int total = dnasequence.length();
        double cumulative = 0.0;

        for (auto& pair : freq_map) {
            arithmeticProb[pair.first] = static_cast<double>(pair.second) / total;

            // Calculate the interval for each nucleotide
            intervals[pair.first] = {cumulative, cumulative + arithmeticProb[pair.first]};
            cumulative += arithmeticProb[pair.first];
        }
    }

    void displayProbabilities() {
        cout << "\nArithmetic Coding Probabilities and Intervals:" << endl;
        for (auto& pair : arithmeticProb) {
            cout << pair.first << " = " << fixed << setprecision(3) << pair.second
                 << " Interval: [" << intervals[pair.first].first << ", "
                 << intervals[pair.first].second << "]" << endl;
        }
    }

    double encodeArithmetic(string dnasequence) {
        double low = 0.0, high = 1.0, range;

        for (char nucleotide : dnasequence) {
            range = high - low;
            high = low + range * intervals[nucleotide].second;
            low = low + range * intervals[nucleotide].first;
        }

        cout << "\nArithmetic Encoded Value: " << fixed << setprecision(6) << (low + high) / 2 << endl;
        return (low + high) / 2;
    }

    string decodeArithmetic(double encodedValue, int sequenceLength) {
        string decodedDNA = "";
        double low = 0.0, high = 1.0;
        cout << "\nArithmetic Decoding Steps:\n";
        for (int i = 0; i < sequenceLength; i++) {
            double range = high - low;
            double value = (encodedValue - low) / range;

            for (auto& pair : intervals) {
                if (value >= pair.second.first && value < pair.second.second) {
                    decodedDNA += pair.first;
                    cout << "Decoded value " << value << " to " << pair.first
                         << " with interval [" << pair.second.first << ", "
                         << pair.second.second << "]\n";
                    high = low + range * pair.second.second;
                    low = low + range * pair.second.first;
                    break;
                }
            }
        }

        return decodedDNA;
    }
};

int main() {
    HuffmanBST huffmanBST;

    // Example DNA sequence
    string dna;
    cout<<"Enter the DNA Sequence: ";
    cin>>dna;
    cout<<"\nCompression: \n\n";

    // Compression
    huffmanBST.countfreq(dna);
    huffmanBST.displayNodes();

    huffmanBST.buildHuffmanTree();
    huffmanBST.storeHuffmanCodes(huffmanBST.root, "");
    huffmanBST.displayHuffmanCodes();

    string huffmanEncoded = huffmanBST.encodeDNA(dna);
    cout << "\nHuffman Encoded DNA Sequence: " << huffmanEncoded << endl;

    huffmanBST.calculateProbabilities(dna);
    huffmanBST.displayProbabilities();
    double arithmeticEncoded = huffmanBST.encodeArithmetic(dna);

    cout<<"\n------------------------------------\n";
    cout<< "\nDecompression: \n";
    string decodedHuffman = huffmanBST.decodeHuffman(huffmanEncoded);
    cout << "\nHuffman Decoded DNA Sequence: " << decodedHuffman << endl;

    string decodedArithmetic = huffmanBST.decodeArithmetic(arithmeticEncoded, dna.length());
    cout << "\nArithmetic Decoded DNA Sequence: " << decodedArithmetic << endl;

    return 0;
}
