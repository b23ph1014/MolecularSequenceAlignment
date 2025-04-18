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