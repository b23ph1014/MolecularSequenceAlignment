#include <iostream>
#include <vector>
#include <unordered_map>
#include <stack>
#include <algorithm>
#include <map>
#include <fstream> // For file reading
#include <chrono>  // For measuring runtime 
#include <climits> // For INT_MIN

using namespace std;

// Needleman-Wunsch Algorithm (Global Alignment Scoring Matrix)
vector<vector<int> > needlemanWunsch(const string &seq1, const string &seq2, int match = 1, int mismatch = -1, int gap_penalty = -1) {
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

// Smith-Waterman Algorithm (Local Alignment)
vector<vector<int> > smithWaterman(const string &seq1, const string &seq2, int match = 1, int mismatch = -1, int gap_penalty = -1) {
    int m = seq1.size();
    int n = seq2.size();
    vector<vector<int> > score_matrix(m + 1, vector<int>(n + 1, 0)); // Initialize zeroes

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match_score = (seq1[i - 1] == seq2[j - 1]) ? match : mismatch;
            score_matrix[i][j] = max({0,
                                           score_matrix[i - 1][j - 1] + match_score, // Diagonal
                                           score_matrix[i - 1][j] + gap_penalty,     // Up
                                           score_matrix[i][j - 1] + gap_penalty});   // Left
        }
    }

    return score_matrix;
}
