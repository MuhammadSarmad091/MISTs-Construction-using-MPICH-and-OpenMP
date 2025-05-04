#include <iostream>
#include <bits/stdc++.h>
#include <cstdint>
using namespace std;

static vector< vector<uint8_t> > perms;      // all vertices, each permutation symbol in [1..n]
static vector< vector<uint8_t> > pos;        // pos[v][symbol] = index of symbol in perms[v]
static vector<uint8_t> firstWrong;           // firstWrong[v] = first wrong-from-right position, 1-based
static vector<uint8_t> root;                 // identity permutation [1..n]

vector< vector<uint8_t> > generatePermutations(int n) {
    vector<uint8_t> base(n);
    iota(base.begin(), base.end(), 1);
    vector<vector<uint8_t>> all;
    do all.push_back(base);
    while(next_permutation(base.begin(), base.end()));
    return all;
}

string permToString(const vector<uint8_t>& p) {
    string s;
    s.reserve(p.size());
    for (uint8_t x : p) s.push_back(char('0' + x));
    return s;
}



void preprocess(int n) {
    size_t N = perms.size();
    for (size_t v = 0; v < N; ++v) {
        auto &P = perms[v];
        for (int j = 0; j < n; ++j)
            pos[v][ P[j] ] = (uint8_t)j;
            int r;
            for (r = n - 1; r >= 0; --r) if (P[r] != r+1) break;
        firstWrong[v] = (r <= 0 ? 1 : (uint8_t)r);
    }
}

vector<uint8_t> swapAdjacent(size_t vIdx, uint8_t symbol) {
    auto v = perms[vIdx];
    int j = pos[vIdx][symbol];
    if (j < 0 || j + 1 >= (int)v.size()) return v;
    swap(v[j], v[j+1]);
    return v;
}

vector<uint8_t> findPosition(size_t vIdx, int t, int n) {
    auto v = perms[vIdx];
    auto u = swapAdjacent(vIdx, (uint8_t)t);
    if (t == 2 && u == root) return swapAdjacent(vIdx, (uint8_t)(t-1));
    uint8_t vn1 = v[n-2];
    if (vn1 == t || vn1 == n-1) return swapAdjacent(vIdx, (uint8_t)(firstWrong[vIdx] + 1));
    return u;
}

vector<uint8_t> parent1(size_t vIdx, int t, int n) {
    auto v = perms[vIdx];
    uint8_t vn = v[n-1], vn1 = v[n-2];
    if (vn == n) {
        return (t != n-1 ? findPosition(vIdx, t, n)
                         : swapAdjacent(vIdx, vn1));
    }
    if (vn == n-1 && vn1 == n && swapAdjacent(vIdx, (uint8_t)n) != root) {
        return (t == 1 ? swapAdjacent(vIdx, (uint8_t)n)
                       : swapAdjacent(vIdx, (uint8_t)(t-1)));
    }
    return (vn == t ? swapAdjacent(vIdx, (uint8_t)n)
                    : swapAdjacent(vIdx, (uint8_t)t));
}

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
    clock_t start = clock();
    // Generate and store all permutations of size n
    perms = generatePermutations(n);
    int N = perms.size();

    // Prepare identity root globally
    root.resize(n);
    iota(root.begin(), root.end(), static_cast<uint8_t>(1));

    // Preprocess per-vertex information
    pos.assign(N, vector<uint8_t>(n+1));
    firstWrong.assign(N, static_cast<uint8_t>(1));
    preprocess(n);

    // Map permutation string to index for lookup
    unordered_map<string, int> indexOf;
    indexOf.reserve(N);
    for (int i = 0; i < N; i++) {
        indexOf[ permToString(perms[i]) ] = i;
    }

    // Prepare children lists: children[t][pIdx] -> list of child indices
    vector<vector<vector<int>>> children(n-1, vector<vector<int>>(N));

    // Identify root index (identity permutation)
    int rootIndex = indexOf[ permToString(root) ];

    // Build all trees using precomputed structures
    for (uint8_t t = 1; t <= n-1; t++) {
        for (int vIdx = 0; vIdx < N; vIdx++) {
            if (vIdx == rootIndex) continue;
            vector<uint8_t> pVec = parent1(vIdx, t, n);
            int pIdx = indexOf[ permToString(pVec) ];
            children[t-1][pIdx].push_back(vIdx);
        }
    }

    clock_t end = clock();
    double time = (double)(end-start)/CLOCKS_PER_SEC;
    cout<<"Time : "<<time<<endl;

    /*

    // Export to DOT files
    for (uint8_t t = 1; t <= n-1; t++) {
        string fname = "Tn_" + to_string(t) + ".dot";
        ofstream dot(fname);
        dot << "digraph T" << n << "_" << static_cast<int>(t) << " {\n  rankdir=TB;\n";
        for (int pIdx = 0; pIdx < N; pIdx++) {
            string pLabel = permToString(perms[pIdx]);
            for (int c : children[t-1][pIdx]) {
                string cLabel = permToString(perms[c]);
                dot << "  \"" << pLabel << "\" -> \"" << cLabel << "\";\n";
            }
        }
        dot << "}\n";
    }
        */
    return 0;
}

