// competing_risks - Dynamic Programming Variant 3 (d3314)
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

namespace dp {

class DPSolver_3_3314 {
public:
    // Longest common subsequence
    int lcs(const std::string& a, const std::string& b) {
        int m = a.size(), n = b.size();
        std::vector<std::vector<int>> table(m+1, std::vector<int>(n+1, 0));
        for (int i = 1; i <= m; i++)
            for (int j = 1; j <= n; j++)
                table[i][j] = (a[i-1] == b[j-1])
                    ? table[i-1][j-1] + 1
                    : std::max(table[i-1][j], table[i][j-1]);
        return table[m][n];
    }

    // 0/1 Knapsack
    int knapsack(const std::vector<int>& weights, const std::vector<int>& values, int capacity) {
        int n = weights.size();
        std::vector<std::vector<int>> K(n+1, std::vector<int>(capacity+1, 0));
        for (int i = 1; i <= n; i++)
            for (int w = 0; w <= capacity; w++) {
                K[i][w] = K[i-1][w];
                if (weights[i-1] <= w)
                    K[i][w] = std::max(K[i][w], K[i-1][w - weights[i-1]] + values[i-1]);
            }
        return K[n][capacity];
    }

    // Edit distance
    int edit_distance(const std::string& s1, const std::string& s2) {
        int m = s1.size(), n = s2.size();
        std::vector<std::vector<int>> dp(m+1, std::vector<int>(n+1));
        for (int i = 0; i <= m; i++) dp[i][0] = i;
        for (int j = 0; j <= n; j++) dp[0][j] = j;
        for (int i = 1; i <= m; i++)
            for (int j = 1; j <= n; j++)
                dp[i][j] = (s1[i-1] == s2[j-1])
                    ? dp[i-1][j-1]
                    : 1 + std::min({dp[i-1][j], dp[i][j-1], dp[i-1][j-1]});
        return dp[m][n];
    }

    // Matrix chain multiplication
    int matrix_chain(const std::vector<int>& dims) {
        int n = dims.size() - 1;
        std::vector<std::vector<int>> m(n, std::vector<int>(n, 0));
        for (int len = 2; len <= n; len++)
            for (int i = 0; i <= n - len; i++) {
                int j = i + len - 1;
                m[i][j] = INT_MAX;
                for (int k = i; k < j; k++)
                    m[i][j] = std::min(m[i][j],
                        m[i][k] + m[k+1][j] + dims[i]*dims[k+1]*dims[j+1]);
            }
        return m[0][n-1];
    }
};

} // namespace dp
