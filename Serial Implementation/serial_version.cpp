#include <bits/stdc++.h>
using namespace std;

// Serial construction of n-1 independent spanning trees (ISTs) of the bubble-sort network B_n.
// First we store each tree in memory as adjacency lists, then export to DOT files.

// Function prototypes
int factorial(int n);
vector<vector<int>> generatePermutations(int n);
string permToString(const vector<int>& p);
vector<int> swapAdjacent(const vector<int>& v, int pos);
vector<int> findPosition(const vector<int>& v, int t, int n);
vector<int> parent1(const vector<int>& v, int t, int n);

int main(int argc, char** argv) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <n (2<=n<=10)>\n";
        return 1;
    }
    int n = stoi(argv[1]);
    if (n < 2 || n > 10) {
        cerr << "n must be in range [2..10]\n";
        return 1;
    }

    // Generate all permutations of size n
    vector<vector<int>> perms = generatePermutations(n);
    int N = perms.size();

    // Map permutation string to index
    unordered_map<string, int> indexOf;
    indexOf.reserve(N);
    for (int i = 0; i < N; i++) {
        string key = permToString(perms[i]);
        indexOf[key] = i;
    }

    // Prepare storage: children[t][i] = list of child-indices of vertex i in tree T^n_{t+1}
    vector<vector<vector<int>>> children(n - 1, vector<vector<int>>(N));

    // Identity root label
    vector<int> rootVec(n);
    iota(rootVec.begin(), rootVec.end(), 1);
    string rootLabel = permToString(rootVec);
    int rootIndex = indexOf[rootLabel];

    // Build trees in memory
    for (int t = 1; t <= n - 1; t++) {
        for (int i = 0; i < N; i++) {
            string vLabel = permToString(perms[i]);
            if (i == rootIndex) continue;  // skip root
            vector<int> pVec = parent1(perms[i], t, n);
            string pLabel = permToString(pVec);
            int pIndex = indexOf[pLabel];
            children[t - 1][pIndex].push_back(i);
        }
    }

    // Export each tree to DOT
    for (int t = 1; t <= n - 1; t++) {
        string dotName = "Tn_" + to_string(t) + ".dot";
        ofstream dot(dotName);
        dot << "digraph T" << n << "_" << t << " {\n";
        dot << "  rankdir=TB;\n";
        // emit edges parent->child
        for (int p = 0; p < N; p++) {
            string pLabel = permToString(perms[p]);
            for (int cIdx = 0; cIdx < children[t - 1][p].size(); ++cIdx) {
                int c = children[t - 1][p][cIdx];
                string cLabel = permToString(perms[c]);
                dot << "  \"" << pLabel << "\" -> \"" << cLabel << "\";\n";
            }
        }
        dot << "}\n";
        dot.close();
        cout << "Wrote " << dotName << "\n";
    }
    return 0;
}

int factorial(int n) {
    int f = 1;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

vector<vector<int>> generatePermutations(int n) {
    vector<int> base(n);
    iota(base.begin(), base.end(), 1);
    vector<vector<int>> all;
    do { all.push_back(base); } while (next_permutation(base.begin(), base.end()));
    return all;
}

string permToString(const vector<int>& p) {
    string s;
    for (int i = 0; i < p.size(); ++i) {
        s += char('0' + p[i]);
    }
    return s;
}

vector<int> swapAdjacent(const vector<int>& v, int symbol) {
    // Find the position of the given symbol in the permutation, then swap it with its successor
    int n = v.size();
    vector<int> u = v;
    // locate index j where v[j] == symbol
    int j = -1;
    for (int idx = 0; idx < n; ++idx) {
        if (v[idx] == symbol) { j = idx; break; }
    }
    if (j < 0 || j >= n-1) {
        // symbol not found or already at last position—no change
        return u;
    }
    std::swap(u[j], u[j+1]);
    return u;
}

vector<int> findPosition(const vector<int>& v, int t, int n) {
    vector<int> u = swapAdjacent(v, t);
    vector<int> root(n);
    iota(root.begin(), root.end(), 1);
    if (t == 2 && u == root) return swapAdjacent(v, t - 1);
    int vn1 = v[n - 2];
    if (vn1 == t || vn1 == n - 1) {
        int r;
        for (r = n - 1; r >= 0; --r) if (v[r] != r+1) break;
        return swapAdjacent(v, r+1);
    }
    return swapAdjacent(v, t);
}

vector<int> parent1(const vector<int>& v, int t, int n) {
    int vn = v[n - 1];
    int vn1 = v[n - 2];
    vector<int> p;
    vector<int> root(n);
    iota(root.begin(), root.end(), 1);

    if (vn == n) {
        if (t != n - 1) p = findPosition(v, t, n);
        else p = swapAdjacent(v, vn1);
    } else {
        if (vn == n - 1 && vn1 == n && swapAdjacent(v, n) != root) {
            if (t == 1) p = swapAdjacent(v, n);
            else p = swapAdjacent(v, t - 1);
        } else {
            if (vn == t) p = swapAdjacent(v, n);
            else p = swapAdjacent(v, t);
        }
    }
    return p;
}
