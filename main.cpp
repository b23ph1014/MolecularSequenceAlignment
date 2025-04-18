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

