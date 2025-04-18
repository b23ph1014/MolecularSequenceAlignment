#include <iostream>
#include <vector>
#include <algorithm>
#include <sys/resource.h>  // For getrusage
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
    
    //Beginning with bottom-right corner of the matrix
    while (i > 0 || j > 0) {
        // Move diagonally if characters match (or mismatch)
        if (i > 0 && j > 0 &&  
            matrix[i][j] == matrix[i-1][j-1] + 
            ((seq1[i-1] == seq2[j-1]) ? match : mismatch)) {
            aligned1 = seq1[i-1] + aligned1;
            aligned2 = seq2[j-1] + aligned2;
            i--; j--;
        }
        // Move up if gap in seq2
        else if (i > 0 && matrix[i][j] == matrix[i-1][j] + gap_penalty) {
            aligned1 = seq1[i-1] + aligned1;
            aligned2 = '-' + aligned2;
            i--;
        }
        // Move left if gap in seq1
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

    out << "Sequence 1: " << aligned1 << "\n"; // Aligned sequence 1
    out << "            " << alignment_line << "\n"; // Alignment visualization
    out << "Sequence 2: " << aligned2 << "\n"; // Aligned sequence 2
}

// Smith-Waterman Algorithm (Local Alignment with traceback)
pair<string, string> smithWaterman(const string &seq1, const string &seq2, int match, int mismatch, int gap_penalty) {
    int m = seq1.size();
    int n = seq2.size();

    vector<vector<int> > score_matrix(m + 1, vector<int>(n + 1, 0)); // Initializing matrix with 0s
    int max_i = 0, max_j = 0, max_score = 0;

    // Fill scoring matrix based on match/mismatch/gap rules
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int match_score = (seq1[i - 1] == seq2[j - 1]) ? match : mismatch;
            score_matrix[i][j] = max({0,
                score_matrix[i - 1][j - 1] + match_score, // Diagonal
                score_matrix[i - 1][j] + gap_penalty, // Up
                score_matrix[i][j - 1] + gap_penalty // Left
            });
            // Tracking highest score's position for traceback
            if (score_matrix[i][j] > max_score) {
                max_score = score_matrix[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // Start traceback from cell with maximum score
    string aligned1 = "", aligned2 = "";
    int i = max_i, j = max_j;

    while (i > 0 && j > 0 && score_matrix[i][j] != 0) {
        int score = score_matrix[i][j];
        int match_score = (seq1[i - 1] == seq2[j - 1]) ? match : mismatch;

        // If current cell came from a diagonal move (match/mismatch)
        if (score == score_matrix[i - 1][j - 1] + match_score) {
            aligned1 = seq1[i - 1] + aligned1;
            aligned2 = seq2[j - 1] + aligned2;
            i--; j--;
        }
        // If it came from the cell above (gap in seq2) 
        else if (score == score_matrix[i - 1][j] + gap_penalty) {
            aligned1 = seq1[i - 1] + aligned1;
            aligned2 = '-' + aligned2;
            i--;
        } 
        // Otherwise, it came from the left (gap in seq1)
        else {
            aligned1 = '-' + aligned1;
            aligned2 = seq2[j - 1] + aligned2;
            j--;
        }
    }

    cout << "\nSmith-Waterman Alignment\n";
    cout << "Alignment Score: " << max_score << "\n\n";
    cout << "Visual Alignment:\n";
    printAlignmentWithVisualization(aligned1, aligned2);

    return {aligned1, aligned2};
}

// Computes the last row of N-W matrix...used in Hirschberg's algorithm
vector<int> nwScore(const string &a, int a_start, int a_end,
    const string &b, int b_start, int b_end,
    int match, int mismatch, int gap) {
    int m = a_end - a_start;
    int n = b_end - b_start;

    vector<int> prev(n + 1), curr(n + 1);
    
    // Initializing first row
    for (int j = 0; j <= n; ++j)
        prev[j] = j * gap;
    
    // Fill matrix row-by-row, keeping only current and previous row
    for (int i = 1; i <= m; ++i) {
        curr[0] = i * gap;
        for (int j = 1; j <= n; ++j) {
            int cost = (a[a_start + i - 1] == b[b_start + j - 1]) ? match : mismatch;
            curr[j] = max({prev[j - 1] + cost,
                    prev[j] + gap,
                    curr[j - 1] + gap});
        }
        swap(prev, curr);
    }

    return prev; // Only last row will be returned
}

// Hirschberg's Algorithm: memory-efficient global alignment
void hirschberg(const string &a, int a_start, int a_end,
    const string &b, int b_start, int b_end,
    string &res_a, string &res_b,
    int match, int mismatch, int gap) {

    int m = a_end - a_start;
    int n = b_end - b_start;
    
    // Base case: one string is empty
    if (m == 0) {
        res_a += string(n, '-');
        res_b += b.substr(b_start, n);
    } else if (n == 0) {
        res_a += a.substr(a_start, m);
        res_b += string(m, '-');
    } 
    // Small enough to use standard alignment
    else if (m == 1 || n == 1) {
        // Fall back to standard NW alignment
        string sub_a = a.substr(a_start, m);
        string sub_b = b.substr(b_start, n);
        auto matrix = needlemanWunsch(sub_a, sub_b, match, mismatch, gap);
        auto aligned = traceback(sub_a, sub_b, matrix, match, mismatch, gap);
        res_a += aligned.first;
        res_b += aligned.second;
    } 
    // Recursive case: split problem in half
    else {
        int mid = a_start + m / 2;

        auto scoreL = nwScore(a, a_start, mid, b, b_start, b_end, match, mismatch, gap);
        auto scoreR = nwScore(a, mid, a_end, b, b_start, b_end, match, mismatch, gap);
        reverse(scoreR.begin(), scoreR.end());

        // Find partition point
        int max_j = 0, max_score = INT_MIN;
        for (int j = 0; j <= n; ++j) {
        int score = scoreL[j] + scoreR[n - j];
        if (score > max_score) {
            max_score = score;
            max_j = j;
            }
        }

        int b_mid = b_start + max_j;

        // Recursively align left and right halves
        hirschberg(a, a_start, mid, b, b_start, b_mid, res_a, res_b, match, mismatch, gap);
        hirschberg(a, mid, a_end, b, b_mid, b_end, res_a, res_b, match, mismatch, gap);
    }
}

// Helper function to calculate alignment score based on match, mismatch, gap penalty values (used in Hirschberg's algorithm)
int calculateAlignmentScore(const string &aligned1, const string &aligned2, int match, int mismatch, int gap) {
    int score = 0;

    // Loop through each aligned character in a sequence (seq1 here)
    for (int i = 0; i < aligned1.size(); ++i) {
        // If there's a gap in either sequence, apply gap penalty
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

// Returns the peak memory usage (in KB) of the current process
long long printMemoryUsage() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage); // Get emory usage details

    return usage.ru_maxrss / 1024; // Divided by 1024 to get the result in KB
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
    int match = 1, mismatch = -1, gap = -2;

    // Values for match, mismatch, & gap (chosen by the user)
    if (argc > 3) match = stoi(argv[3]);
    if (argc > 4) mismatch = stoi(argv[4]);
    if (argc > 5) gap = stoi(argv[5]);

    // Needleman
    auto start = chrono::high_resolution_clock::now();
    auto matrix = needlemanWunsch(seq1, seq2, match, mismatch, gap);
    auto aligned = traceback(seq1, seq2, matrix, match, mismatch, gap);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;

    cout << "\nNeedleman-Wunsch Alignment (Global)\n";
    cout << "Alignment Score: " << matrix[seq1.size()][seq2.size()] << "\n\n";

    cout << "Visual Alignment:\n";
    printAlignmentWithVisualization(aligned.first, aligned.second);

    cout << "\nRuntime: " << duration.count() << " seconds\n";
    long long NeedlemanMemoryUsage = printMemoryUsage();

    // Hirschberg
    string h_align1, h_align2;
    auto start_hirschberg = chrono::high_resolution_clock::now();
    hirschberg(seq1, 0, seq1.size(), seq2, 0, seq2.size(), h_align1, h_align2, match, mismatch, gap);
    int h_score = calculateAlignmentScore(h_align1, h_align2, match, mismatch, gap);
    auto end_hirschberg = chrono::high_resolution_clock::now();
    chrono::duration<double> duration_hirschberg = end_hirschberg - start_hirschberg;
    cout << "\nHirschberg Alignment" << "\n";
    cout << "Score: " << h_score << "\n\n";
    cout << "Visual Alignment:\n";
    printAlignmentWithVisualization(h_align1, h_align2);
    cout << "\nRuntime: " << duration_hirschberg.count() << " seconds\n";
    long long HirschbergMemoryUsage = printMemoryUsage();

    // Smith
    start = chrono::high_resolution_clock::now();
    smithWaterman(seq1, seq2, match, mismatch, gap);
    end = chrono::high_resolution_clock::now();
    duration = end - start;
    cout << "\nRuntime: " << duration.count() << " seconds\n";
    long long SmithMemoryUsage = printMemoryUsage();

    cout << "\nMemory Usage Comparison: \n";
    cout << "Needleman-Wunsch Memory Usage (in KB): " << NeedlemanMemoryUsage;
    cout << "\nSmith-Waterman Memory Usage (in KB): " << SmithMemoryUsage;
    cout << "\nHirschberg's Algorithm Memory Usage (in KB): " << HirschbergMemoryUsage << endl;

    return 0;
}
