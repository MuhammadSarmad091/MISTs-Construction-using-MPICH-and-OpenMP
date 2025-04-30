#include <bits/stdc++.h>
using namespace std;

// This program constructs the n-1 independent spanning trees (ISTs) of
a bubble-sort network B_n in serial, and exports each tree to a DOT file
// for visualization with Graphviz.

// We enumerate all n! permutations of {1..n} in lexicographic order.
// Each permutation is represented as vector<int> of size n, and also as a string label.

// Function prototypes
int factorial(int n);
vector<vector<int>> generatePermutations(int n);
string permToString(const vector<int>& p);
int findR(const vector<int>& v);
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
    auto perms = generatePermutations(n);
    int N = perms.size();

    // Build a map from permutation string to its index
    unordered_map<string,int> indexOf;
    indexOf.reserve(N);
    for (int i = 0; i < N; i++) {
        indexOf[ permToString(perms[i]) ] = i;
    }

    // For each tree t = 1 .. n-1, compute edges and write DOT
    for (int t = 1; t <= n-1; t++) {
        // Open DOT file
        string dotName = "Tn_" + to_string(t) + ".dot";
        ofstream dot(dotName);
        dot << "digraph T" << n << "_" << t << " {\n";
        dot << "  rankdir=TB;\n";  // top-to-bottom

        // For each vertex v != root, compute its parent p in tree T^n_t
        vector<int> root(n);
        iota(root.begin(), root.end(), 1);
        string rootLabel = permToString(root);

        for (int i = 0; i < N; i++) {
            const auto &v = perms[i];
            string vLabel = permToString(v);
            if (vLabel == rootLabel) continue;
            auto p = parent1(v, t, n);
            string pLabel = permToString(p);
            dot << "  \"" << pLabel << "\" -> \"" << vLabel << "\";\n";
        }
        dot << "}\n";
        dot.close();
        cout << "Wrote " << dotName << "\n";
    }
    return 0;
}

// factorial helper
int factorial(int n) {
    int f = 1;
    for (int i = 2; i <= n; i++) f *= i;
    return f;
}

// generate permutations in lex order
vector<vector<int>> generatePermutations(int n) {
    vector<int> base(n);
    iota(base.begin(), base.end(), 1);
    vector<vector<int>> all;
    do {
        all.push_back(base);
    } while (next_permutation(base.begin(), base.end()));
    return all;
}

string permToString(const vector<int>& p) {
    string s;
    for (int x : p) s += char('0' + x);
    return s;
}

// r(v): position of first symbol from right not in its identity place
int findR(const vector<int>& v) {
    int n = v.size();
    for (int i = n-1; i >= 0; i--) {
        if (v[i] != i+1) return i+1;  // positions are 1-based
    }
    return 1;
}

// SwapAdjacent: swap positions pos and pos+1 (1-based pos)
vector<int> swapAdjacent(const vector<int>& v, int pos) {
    vector<int> u = v;
    int i = pos-1;
    swap(u[i], u[i+1]);
    return u;
}

// FindPosition as in pseudocode
vector<int> findPosition(const vector<int>& v, int t, int n) {
    // if t=2 and swapping at t yields root
    auto u = swapAdjacent(v, t);
    vector<int> root(n);
    iota(root.begin(), root.end(), 1);
    if (t == 2 && u == root) {
        return swapAdjacent(v, t-1);
    }
    int vn1 = v[n-2];
    if (vn1 == t || vn1 == n-1) {
        int r = findR(v);
        return swapAdjacent(v, r-1);
    }
    return swapAdjacent(v, t);
}

// Parent1 implements Algorithm 1 for tree t
vector<int> parent1(const vector<int>& v, int t, int n) {
    int vn = v[n-1];
    vector<int> p;
    vector<int> root(n);
    iota(root.begin(), root.end(), 1);

    if (vn == n) {
        if (t != n-1) p = findPosition(v, t, n);
        else p = swapAdjacent(v, vn-1);
    } else {
        int vn1 = v[n-2];
        if (vn == n-1 && vn1 == n && swapAdjacent(v, n-1) != root) {
            if (t == 1) p = swapAdjacent(v, n-1);
            else p = swapAdjacent(v, t-1);
        } else {
            if (vn == t) p = swapAdjacent(v, n-1);
            else p = swapAdjacent(v, t);
        }
    }
    return p;
}
