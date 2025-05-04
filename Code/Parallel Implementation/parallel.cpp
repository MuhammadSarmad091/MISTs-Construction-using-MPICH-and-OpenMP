#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>
using namespace std;

// Hybrid MPI+OpenMP construction of n-1 independent spanning trees (ISTs)
// of bubble-sort network B_n, using compact data types to reduce memory.

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
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0) cerr << "Usage: " << argv[0] << " <n (2<=n<=10)>\n";
        MPI_Finalize(); return 1;
    }
    int n = stoi(argv[1]);
    if (n < 2 || n > 10) {
        if (rank == 0) cerr << "n must be in range [2..10]\n";
        MPI_Finalize(); return 1;
    }

    // Generate permutations and setup
    perms = generatePermutations(n);
    size_t N = perms.size();
    root.resize(n);
    iota(root.begin(), root.end(), 1);
    pos.assign(N, vector<uint8_t>(n+1));
    firstWrong.assign(N, 1);
    preprocess(n);

    unordered_map<string,int> indexOf;
    indexOf.reserve(N);
    for (size_t i = 0; i < N; ++i)
        indexOf[ permToString(perms[i]) ] = (int)i;

    int T = n - 1;
    int per = T / size, rem = T % size;
    int start_t = (rank < rem ? rank*(per+1)+1 : rem*(per+1)+(rank-rem)*per+1);
    int end_t   = (rank < rem ? start_t+per : start_t+per-1);
    vector<int> assigned_t;
    for (int t = start_t; t <= end_t; ++t) assigned_t.push_back(t);

    size_t M = assigned_t.size() * N;
    vector<tuple<int,uint32_t,uint32_t>> edges;

    #pragma omp parallel
    {
        vector<tuple<int,uint32_t,uint32_t>> buf;
        #pragma omp for nowait
        for (size_t i = 0; i < M; ++i) {
            int li = i / N;
            uint32_t vIdx = i % N;
            if (vIdx == (uint32_t)indexOf[permToString(root)]) continue;
            int t = assigned_t[li];
            vector<uint8_t> p = parent1(vIdx, t, n);
            uint32_t pIdx = indexOf[ permToString(p) ];
            buf.emplace_back(li, pIdx, vIdx);
        }
        #pragma omp critical
        edges.insert(edges.end(), buf.begin(), buf.end());
    }

    // Build local and send to root
    if (rank == 0) {
        vector<vector<vector<uint32_t>>> children_global(T, vector<vector<uint32_t>>(N));
        // own
        for (auto &e : edges) {
            int li; uint32_t p,v;
            tie(li,p,v) = e;
            int t = assigned_t[li] - 1;
            children_global[t][p].push_back(v);
        }
        // receive others
        for (int src = 1; src < size; ++src) {
            int pst = (src < rem ? src*(per+1)+1 : rem*(per+1)+(src-rem)*per+1);
            int pen = (src < rem ? pst+per : pst+per-1);
            for (int t = pst; t <= pen; ++t) {
                int tIdx = t-1;
                int cnt;
                MPI_Recv(&cnt,1,MPI_INT,src,t, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                vector<uint32_t> buf(2*cnt);
                MPI_Recv(buf.data(),2*cnt,MPI_UINT32_T,src,t+T, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int k=0; k<cnt; ++k)
                    children_global[tIdx][ buf[2*k] ].push_back(buf[2*k+1]);
            }
        }
        // export
        for (int t=1; t<=T; ++t) {
            ofstream dot("Tn_"+to_string(t)+".dot");
            dot<<"digraph T"<<n<<"_"<<t<<" {\n  rankdir=TB;\n";
            for (uint32_t p=0; p<N; ++p)
                for (auto c: children_global[t-1][p])
                    dot<<"  \""<<permToString(perms[p])<<"\" -> \""<<permToString(perms[c])<<"\";\n";
            dot<<"}\n";
        }
    } else {
        // send edges
        vector<vector<uint32_t>> buf(assigned_t.size());
        vector<int> cnt(assigned_t.size());
        for (auto &e : edges) {
            int li; uint32_t p,v;
            tie(li,p,v) = e;
            buf[li].push_back(p);
            buf[li].push_back(v);
        }
        for (int i=0; i<assigned_t.size(); ++i) {
            int t = assigned_t[i];
            cnt[i] = buf[i].size()/2;
            MPI_Send(&cnt[i],1,MPI_INT,0,t, MPI_COMM_WORLD);
            MPI_Send(buf[i].data(),2*cnt[i],MPI_UINT32_T,0,t+T, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
