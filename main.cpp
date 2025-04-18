#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream> // For file reading
#include <chrono>  // For measuring runtime 
#include <climits> // For INT_MIN

using namespace std;

// Needleman-Wunsch Algorithm (Global Alignment Scoring Matrix)
vector<vector<int> > needlemanWunsch(const string &seq1, const string &seq2, int match, int mismatch, int gap_penalty) {
    int m = seq1.size();
    int n = seq2.size();
    vector<vector<int> > score_matrix(m + 1, vector<int>(n + 1)); // Score matrix of shape (m+1, n+1)

    for (int i = 0; i <= m; ++i)
        score_matrix[i][0] = i * gap_penalty; // Initializing first col with gap penalties
    for (int j = 0; j <= n; ++j)
        score_matrix[0][j] = j * gap_penalty; // Initializing first row with gap penalties.

    // Filling the score matrix
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match_score = (seq1[i - 1] == seq2[j - 1]) ? match : mismatch;
            score_matrix[i][j] = max({score_matrix[i - 1][j - 1] + match_score,  // Diagonal
                                      score_matrix[i - 1][j] + gap_penalty, // Up
                                      score_matrix[i][j - 1] + gap_penalty}); // Left
        }
    }

    return score_matrix;
}

// Traceback function for Needleman-Wunsch
pair<string, string> traceback(const string &seq1, const string &seq2, 
                              const vector<vector<int> > &matrix, 
                              int match, int mismatch, int gap_penalty) {
    string aligned1, aligned2;
    int i = seq1.size(), j = seq2.size();
    
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && 
            matrix[i][j] == matrix[i-1][j-1] + 
            ((seq1[i-1] == seq2[j-1]) ? match : mismatch)) {
            aligned1 = seq1[i-1] + aligned1;
            aligned2 = seq2[j-1] + aligned2;
            i--; j--;
        }
        else if (i > 0 && matrix[i][j] == matrix[i-1][j] + gap_penalty) {
            aligned1 = seq1[i-1] + aligned1;
            aligned2 = '-' + aligned2;
            i--;
        }
        else {
            aligned1 = '-' + aligned1;
            aligned2 = seq2[j-1] + aligned2;
            j--;
        }
    }
    
    return {aligned1, aligned2};
}

void printAlignmentWithVisualization(const string &aligned1, const string &aligned2, ostream &out = cout) {
    string alignment_line;

    for (size_t i = 0; i < aligned1.size(); ++i) {
        if (aligned1[i] == aligned2[i])
            alignment_line += '.'; // Match
        else if (aligned1[i] != '-' && aligned2[i] != '-')
            alignment_line += '*'; // Mismatch
        else
            alignment_line += ' '; // Gap
    }

    out << "Aligned Sequences:\n";
    out << "Sequence 1: " << aligned1 << "\n"; // Aligned sequence 1
    out << "            " << alignment_line << "\n"; // Alignment visualization
    out << "Sequence 2: " << aligned2 << "\n\n"; // Aligned sequence 2
}

// Smith-Waterman Algorithm (Local Alignment with traceback)
pair<string, string> smithWaterman(const string &seq1, const string &seq2, int match, int mismatch, int gap_penalty, ofstream* out = nullptr) {
    int m = seq1.size();
    int n = seq2.size();
    vector<vector<int> > score_matrix(m + 1, vector<int>(n + 1, 0));
    int max_i = 0, max_j = 0, max_score = 0;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match_score = (seq1[i - 1] == seq2[j - 1]) ? match : mismatch;
            score_matrix[i][j] = max({0,
                score_matrix[i - 1][j - 1] + match_score,
                score_matrix[i - 1][j] + gap_penalty,
                score_matrix[i][j - 1] + gap_penalty
            });
            if (score_matrix[i][j] > max_score) {
                max_score = score_matrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    string aligned1 = "", aligned2 = "";
    int i = max_i, j = max_j;

    while (i > 0 && j > 0 && score_matrix[i][j] != 0) {
        int score = score_matrix[i][j];
        int match_score = (seq1[i - 1] == seq2[j - 1]) ? match : mismatch;
        if (score == score_matrix[i - 1][j - 1] + match_score) {
            aligned1 = seq1[i - 1] + aligned1;
            aligned2 = seq2[j - 1] + aligned2;
            i--; j--;
        } else if (score == score_matrix[i - 1][j] + gap_penalty) {
            aligned1 = seq1[i - 1] + aligned1;
            aligned2 = '-' + aligned2;
            i--;
        } else {
            aligned1 = '-' + aligned1;
            aligned2 = seq2[j - 1] + aligned2;
            j--;
        }
    }

    cout << "\nSmith-Waterman Alignment (Local):\n";
    cout << "Score: " << max_score << "\n";
    cout << "Visual Alignment:\n";
    printAlignmentWithVisualization(aligned1, aligned2);

    if (out) {
        *out << "Smith-Waterman Alignment\n";
        *out << "Score: " << max_score << "\n";
        *out << "Visual Alignment:\n";
        printAlignmentWithVisualization(aligned1, aligned2, *out);
    }

    return {aligned1, aligned2};
}

// Computes the last row of N-W matrix...used in Hirschberg's algorithm
vector<int> nwScore(const string &a, const string &b, int match, int mismatch, int gap) {
    int m = a.size(), n = b.size();
    vector<int> prev(n + 1), curr(n + 1);
    for (int j = 0; j <= n; ++j) prev[j] = j * gap;

    
    for (int i = 1; i <= m; ++i) {
        curr[0] = i * gap;
        for (int j = 1; j <= n; ++j) {
            int cost = (a[i - 1] == b[j - 1]) ? match : mismatch;
            curr[j] = max({prev[j - 1] + cost, prev[j] + gap, curr[j - 1] + gap});
        }
        prev = curr;
    }
    return prev;
}

// Hirschberg's Algorithm: memory-efficient global alignment
void hirschberg(const string &a, const string &b, string &res_a, string &res_b, int match, int mismatch, int gap) {
    int m = a.size(), n = b.size();
  
    // Base cases
    if (m == 0) {
        res_a += string(n, '-');
        res_b += b;
    } else if (n == 0) {
        res_a += a;
        res_b += string(m, '-');
    } else if (m == 1 || n == 1) {
        // Needleman-Wunsch for smaller sequences
        auto matrix = needlemanWunsch(a, b, match, mismatch, gap);
        auto align = traceback(a, b, matrix, match, mismatch, gap);
        res_a += align.first;
        res_b += align.second;
    } else {
        // Recursive divide & conquer
        int mid = m / 2;
        auto scoreL = nwScore(a.substr(0, mid), b, match, mismatch, gap);
        auto scoreR = nwScore(string(a.rbegin(), a.rbegin() + (m - mid)), string(b.rbegin(), b.rend()), match, mismatch, gap);

        // Finding partition point
        int max_j = 0, max_score = INT_MIN;
        for (int j = 0; j <= n; ++j) {
            int val = scoreL[j] + scoreR[n - j];
            if (val > max_score) {
                max_score = val;
                max_j = j;
            }
        }
      
        // Recursively align the halves
        string left_a, left_b, right_a, right_b;
        hirschberg(a.substr(0, mid), b.substr(0, max_j), left_a, left_b, match, mismatch, gap);
        hirschberg(a.substr(mid), b.substr(max_j), right_a, right_b, match, mismatch, gap);

        res_a += left_a + right_a;
        res_b += left_b + right_b;
    }
}


// Helper function to calculate alignment score (used in Hirschberg's algorithm)
int calculateAlignmentScore(const string &aligned1, const string &aligned2, int match, int mismatch, int gap) {
    int score = 0;
    for (int i = 0; i < aligned1.size(); ++i) {
        if (aligned1[i] == '-' || aligned2[i] == '-') {
            score += gap;
        } else if (aligned1[i] == aligned2[i]) {
            score += match;
        } else {
            score += mismatch;
        }
    }
    return score;
}

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 6) {
        cerr << "Required Format: " << argv[0] << " <sequence_file1> <sequence_file2> [match=1] [mismatch=-1] [gap=-1]\n";
        return 1;
    }

    string seq1, seq2;
    ifstream file1(argv[1]), file2(argv[2]);
    if (!file1 || !file2) {
        cerr << "Error reading input files.\n";
        return 1;
    }

    getline(file1, seq1);
    getline(file2, seq2);

    // Values for match, mismatch, & gap (default ones)
    int match = 1, mismatch = -1, gap = -1;

    // Values for match, mismatch, & gap (chosen by the user)
    if (argc > 3) match = stoi(argv[3]);
    if (argc > 4) mismatch = stoi(argv[4]);
    if (argc > 5) gap = stoi(argv[5]);

    ofstream outfile("alignment_output.txt");

    auto start = chrono::high_resolution_clock::now();
    auto matrix = needlemanWunsch(seq1, seq2, match, mismatch, gap);
    auto aligned = traceback(seq1, seq2, matrix, match, mismatch, gap);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    cout << "\nNeedleman-Wunsch Runtime: " << duration.count() << " seconds\n";
    cout << "Alignment Score: " << matrix[seq1.size()][seq2.size()] << "\n";

    cout << "Visual Alignment:\n";
    printAlignmentWithVisualization(aligned.first, aligned.second);

    outfile << "Needleman-Wunsch Alignment\n";
    outfile << "Score: " << matrix[seq1.size()][seq2.size()] << "\n";
    outfile << "Visual Alignment:\n";
    printAlignmentWithVisualization(aligned.first, aligned.second, outfile);

    smithWaterman(seq1, seq2, match, mismatch, gap, &outfile);

    string h_align1, h_align2;
    auto start_hirschberg = chrono::high_resolution_clock::now();
    hirschberg(seq1, seq2, h_align1, h_align2, match, mismatch, gap);
    int h_score = calculateAlignmentScore(h_align1, h_align2, match, mismatch, gap);
    auto end_hirschberg = chrono::high_resolution_clock::now();
    chrono::duration<double> duration_hirschberg = end_hirschberg - start_hirschberg;
    cout << "\nHirschberg Runtime: " << duration_hirschberg.count() << " seconds\n";
    cout << "Hirschberg Alignment Score: " << h_score << "\n";
    cout << "Visual Alignment:\n";
    printAlignmentWithVisualization(h_align1, h_align2);

    outfile << "Hirschberg Alignment\n";
    outfile << "Hirschberg Alignment Score: " << h_score << "\n\n";
    outfile << "Visual Alignment:\n";
    printAlignmentWithVisualization(h_align1, h_align2, outfile);

    return 0;
}
