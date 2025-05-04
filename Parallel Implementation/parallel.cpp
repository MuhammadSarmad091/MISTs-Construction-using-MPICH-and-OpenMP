#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>
using namespace std;

// Hybrid MPI+OpenMP construction of n-1 independent spanning trees (ISTs)
// of bubble-sort network B_n.  Each MPI rank computes a subset of trees;
// within each rank, OpenMP parallelizes over (tree,vertex) pairs.

// Global storage for permutations, preprocessing, and identity root
static vector<vector<int>> perms;
static vector<vector<int>> pos;
static vector<int> firstWrong;
static vector<int> root;

// Function prototypes
vector<vector<int>> generatePermutations(int n);
string permToString(const vector<int>& p);
void preprocess(int n);
vector<int> swapAdjacent(int vIdx, int symbol);
vector<int> findPosition(int vIdx, int t, int n);
vector<int> parent1(int vIdx, int t, int n);

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0) cerr << "Usage: " << argv[0] << " <n (2<=n<=10)>\n";
        MPI_Finalize();
        return 1;
    }
    int n = stoi(argv[1]);
    if (n < 2 || n > 10) {
        if (rank == 0) cerr << "n must be in range [2..10]\n";
        MPI_Finalize();
        return 1;
    }

    // Generate and store all permutations of size n
    perms = generatePermutations(n);
    int N = perms.size();

    // Prepare identity root globally
    root.resize(n);
    iota(root.begin(), root.end(), 1);

    // Preprocess per-vertex information
    pos.assign(N, vector<int>(n+1));
    firstWrong.assign(N, 1);
    preprocess(n);

    // Map permutation string to index
    unordered_map<string, int> indexOf;
    indexOf.reserve(N);
    for (int i = 0; i < N; i++) {
        indexOf[ permToString(perms[i]) ] = i;
    }

    // Determine tree assignment per MPI rank
    int T = n - 1;
    int per = T / size, rem = T % size;
    int start_t, end_t;
    if (rank < rem) {
        start_t = rank * (per + 1) + 1;
        end_t = start_t + per;
    } else {
        start_t = rem * (per + 1) + (rank - rem) * per + 1;
        end_t = start_t + per - 1;
    }
    vector<int> assigned_t;
    for (int t = start_t; t <= end_t; t++) assigned_t.push_back(t);

    // Local children lists for assigned trees
    vector<vector<vector<int>>> children_local(assigned_t.size(), vector<vector<int>>(N));

    // OpenMP parallel over (tree,vertex)
    omp_set_num_threads(3);
    int M = assigned_t.size() * N;
    vector<tuple<int,int,int>> threadEdges;
    #pragma omp parallel
    {
        vector<tuple<int,int,int>> buf;
        #pragma omp for schedule(static)
        for (int i = 0; i < M; i++) {
            int localIdx = i / N;
            int vIdx     = i % N;
            if (vIdx == indexOf[ permToString(root) ]) continue;
            int t = assigned_t[localIdx];
            int pIdx = indexOf[ permToString(parent1(vIdx, t, n)) ];
            buf.emplace_back(localIdx, pIdx, vIdx);
        }
        #pragma omp critical
        threadEdges.insert(threadEdges.end(), buf.begin(), buf.end());
    }
    // Build local children from threadEdges
    for (auto &e : threadEdges) {
        int li,p,v; tie(li,p,v) = e;
        children_local[li][p].push_back(v);
    }

    // MPI gather to root
    if (rank == 0) {
        vector<vector<vector<int>>> children_global(T, vector<vector<int>>(N));
        // copy own
        for (int i = 0; i < assigned_t.size(); i++)
            children_global[assigned_t[i]-1] = children_local[i];
        // receive others
        for (int src = 1; src < size; src++) {
            int pstart, pend;
            if (src < rem) { pstart = src*(per+1)+1; pend = pstart+per; }
            else { pstart = rem*(per+1)+(src-rem)*per+1; pend = pstart+per-1; }
            for (int t = pstart; t <= pend; t++) {
                int tIdx = t-1;
                int count;
                MPI_Recv(&count,1,MPI_INT,src,t, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                vector<int> buf(2*count);
                MPI_Recv(buf.data(),2*count,MPI_INT,src,t+T, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int k=0; k<count; k++) {
                    children_global[tIdx][ buf[2*k] ].push_back(buf[2*k+1]);
                }
            }
        }
        // Export global to DOT
        for (int t=1; t<=T; t++) {
            ofstream dot("Tn_"+to_string(t)+".dot");
            dot<<"digraph T"<<n<<"_"<<t<<" {\n  rankdir=TB;\n";
            for (int p=0; p<N; p++) for (int c: children_global[t-1][p])
                dot<<"  \""<<permToString(perms[p])<<"\" -> \""<<permToString(perms[c])<<"\";\n";
            dot<<"}\n";
        }
    } else {
        // send local to root
        for (int i=0; i<assigned_t.size(); i++) {
            int t = assigned_t[i], tIdx = t-1;
            vector<int> buf;
            for (int p=0; p<N; p++) for (int v: children_local[i][p]) {
                buf.push_back(p); buf.push_back(v);
            }
            int cnt=buf.size()/2;
            MPI_Send(&cnt,1,MPI_INT,0,t, MPI_COMM_WORLD);
            MPI_Send(buf.data(),2*cnt,MPI_INT,0,t+T, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}

vector<vector<int>> generatePermutations(int n) {
    vector<int> base(n); iota(base.begin(),base.end(),1);
    vector<vector<int>> all;
    do all.push_back(base);
    while(next_permutation(base.begin(),base.end()));
    return all;
}

string permToString(const vector<int>& p) {
    string s;
    for(int x:p) s.push_back(char('0'+x));
    return s;
}

void preprocess(int n) {
    int N = perms.size();
    for (int vIdx = 0; vIdx < N; vIdx++) {
        const vector<int> &v = perms[vIdx];
        // positions of each symbol
        for (int j = 0; j < n; j++) {
            pos[vIdx][ v[j] ] = j;    // 0-based index
        }
        // first wrong-from-right (1-based)
        int r = n;
        for (r = n - 1; r >= 0; --r) if (v[r] != r+1) break;
        firstWrong[vIdx] = (r<=0?1:r);
    }
}

vector<int> swapAdjacent(int vIdx,int symbol){
    auto v=perms[vIdx]; int j=pos[vIdx][symbol];
    if(j<0||j+1>=v.size()) return v;
    swap(v[j],v[j+1]); return v;
}

vector<int> findPosition(int vIdx,int t,int n){
    auto v=perms[vIdx], u=swapAdjacent(vIdx,t);
    if(t==2&&u==root) return swapAdjacent(vIdx,t-1);
    int vn1=v[n-2];
    if(vn1==t||vn1==n-1) return swapAdjacent(vIdx, firstWrong[vIdx]+1);
    return u;
}

vector<int> parent1(int vIdx,int t,int n){
    auto v=perms[vIdx]; int vn=v[n-1], vn1=v[n-2];
    if(vn==n) return (t!=n-1?findPosition(vIdx,t,n):swapAdjacent(vIdx,vn1));
    if(vn==n-1&&vn1==n && swapAdjacent(vIdx,n)!=root)
        return (t==1?swapAdjacent(vIdx,n):swapAdjacent(vIdx,t-1));
    return (vn==t?swapAdjacent(vIdx,n):swapAdjacent(vIdx,t));
}
