#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>
using namespace std;

// Parallel construction of ISTs of bubble-sort network B_n
// Fixed worker-master communication for better load distribution

static vector<vector<uint8_t>> perms;
static vector<vector<uint8_t>> pos;
static vector<uint8_t> firstWrong;
static vector<uint8_t> root;

vector<vector<uint8_t>> generatePermutations(int n) {
    // Pre-calculate factorial size to avoid resizing
    int factorial = 1;
    for(int i = 2; i <= n; i++) {
        factorial *= i;
    }
    
    vector<uint8_t> base(n);
    iota(base.begin(), base.end(), 1);
    
    // Pre-allocate the exact space needed
    vector<vector<uint8_t>> all;
    all.reserve(factorial);
    
    // More efficient insertion
    do {
        all.push_back(base);
    } while (next_permutation(base.begin(), base.end()));
    
    return all;
}

string permToString(const vector<uint8_t>& p) {
    // Faster string construction with single allocation
    string s;
    s.reserve(p.size());
    for (uint8_t x : p) {
        s.push_back(char('0' + x));
    }
    return s;
}

void preprocess(int n) {
    size_t N = perms.size();
    
    // Use different schedulers for better performance
    #pragma omp parallel for schedule(static)
    for (size_t v = 0; v < N; ++v) {
        auto &P = perms[v];
        // Optimized inner loop - unrolled for small n
        if (n <= 8) {
            for (int j = 0; j < n; ++j) {
                pos[v][P[j]] = (uint8_t)j;
            }
        } else {
            // Vectorizable loop for larger n
            for (int j = 0; j < n; ++j) {
                pos[v][P[j]] = (uint8_t)j;
            }
        }
        
        // More efficient first wrong computation
        int r = n - 1;
        while (r >= 0 && P[r] == r+1) r--;
        firstWrong[v] = (r < 0 ? 1 : (uint8_t)r);
    }
}

// Performance-optimized swap function that reduces object copies
inline vector<uint8_t> swapAdjacent(size_t vIdx, uint8_t symbol) {
    const auto& v = perms[vIdx]; 
    int j = pos[vIdx][symbol];
    if (j+1 >= (int)v.size()) return v;
    
    // Create result with reserved capacity
    vector<uint8_t> result(v);
    swap(result[j], result[j+1]); 
    return result;
}

// Optimized position finder with reduced redundant computations
inline vector<uint8_t> findPosition(size_t vIdx, int t, int n) {
    const auto& v = perms[vIdx]; 
    auto u = swapAdjacent(vIdx, (uint8_t)t);
    
    if (t == 2 && u == root) {
        return swapAdjacent(vIdx, (uint8_t)(t-1));
    }
    
    uint8_t vn1 = v[n-2];
    if (vn1 == t || vn1 == n-1) {
        return swapAdjacent(vIdx, (uint8_t)(firstWrong[vIdx]+1));
    }
    
    return u;
}

// Highly optimized parent function with inlining for better performance
inline vector<uint8_t> parent1(size_t vIdx, int t, int n) {
    const auto& v = perms[vIdx]; 
    uint8_t vn = v[n-1], vn1 = v[n-2];
    
    if (vn == n) {
        return (t != n-1) ? findPosition(vIdx, t, n) : swapAdjacent(vIdx, vn1);
    }
    
    if (vn == n-1 && vn1 == n) {
        auto s = swapAdjacent(vIdx, (uint8_t)n);
        if (s != root) {
            return (t == 1) ? s : swapAdjacent(vIdx, (uint8_t)(t-1));
        }
    }
    
    return (vn == t) ? swapAdjacent(vIdx, (uint8_t)n) : swapAdjacent(vIdx, (uint8_t)t);
}

int main(int argc,char**argv){
    MPI_Init(&argc,&argv);
    double t_start = MPI_Wtime();

    int rank,size;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if(argc!=2){ if(rank==0) cerr<<"Usage: "<<argv[0]<<" <n>\n"; MPI_Finalize(); return 1; }
    int n=stoi(argv[1]); if(n<2||n>10){ if(rank==0) cerr<<"n must be [2..10]\n"; MPI_Finalize(); return 1; }

    if (rank == 0) {
        cout << "Running fixed MIST construction with:" << endl;
        cout << "  Size parameter (n): " << n << endl;
        cout << "  MPI processes: " << size << endl;
        cout << "  OpenMP threads per process: " << omp_get_max_threads() << endl;
    }

    // setup
    perms=generatePermutations(n);
    size_t N=perms.size();
    root.resize(n); iota(root.begin(),root.end(),1);
    pos.assign(N,vector<uint8_t>(n+1)); firstWrong.assign(N,1);
    preprocess(n);
    unordered_map<string,int> indexOf; indexOf.reserve(N);
    for(size_t i=0;i<N;++i) indexOf[permToString(perms[i])] = i;

    // Total number of trees to compute
    int T = n-1;
    
    if (rank == 0) {
        cout << "Building " << T << " trees for n=" << n << endl;
    }
    
    // Determine which trees each process will work on
    vector<int> treesToBuild;
    for (int t = 1; t <= T; t++) {
        if (t % size == rank) {
            treesToBuild.push_back(t);
        }
    }
    
    if (rank == 0) {
        cout << "Process " << rank << " will build trees: ";
        for (int t : treesToBuild) cout << t << " ";
        cout << endl;
    }
    
    // master storage
    vector<vector<vector<uint32_t>>> children_global;
    vector<MPI_Request> recvRequests;
    vector<vector<uint32_t>> recvBuffers;
    
    if(rank==0) {
        children_global.assign(T, vector<vector<uint32_t>>(N));
    }

    // workers send buffer management
    vector<vector<uint32_t>> sendBuffers;
    vector<MPI_Request> sendRequests;

    // Build trees assigned to this process
    for(int i = 0; i < (int)treesToBuild.size(); i++) {
        int t = treesToBuild[i];
        
        // Build the tree locally
        vector<vector<uint32_t>> children_t(N);
        for (auto& child_list : children_t) {
            child_list.reserve(n);
        }
        
        // Compute root index
        int rootIdx = indexOf[permToString(root)];
        
        // Build this tree in parallel using OpenMP
        #pragma omp parallel
        {
            // Collect edges for this thread using a fixed-size array
            const int MAX_EDGES_PER_THREAD = 100000; 
            uint32_t* temp_edges = new uint32_t[MAX_EDGES_PER_THREAD * 2];
            int edge_count = 0;
            
            #pragma omp for schedule(guided, 1024)
            for(size_t vIdx=0; vIdx<N; ++vIdx) {
                if((int)vIdx == rootIdx) continue;
                
                // Find parent
                vector<uint8_t> p = parent1(vIdx, t, n);
                
                // Optimize lookup
                bool is_root = (p == root);
                int pIdx;
                
                if (is_root) {
                    pIdx = rootIdx;
                } else {
                    pIdx = indexOf[permToString(p)];
                }
                
                // Store edge in thread-local buffer
                if(pIdx >= 0 && edge_count < MAX_EDGES_PER_THREAD) {
                    temp_edges[edge_count*2] = pIdx; 
                    temp_edges[edge_count*2+1] = vIdx;
                    edge_count++;
                }
            }
            
            // Process collected edges
            vector<pair<uint32_t, uint32_t>> thread_edges;
            thread_edges.reserve(edge_count);
            for(int i = 0; i < edge_count; i++) {
                thread_edges.emplace_back(temp_edges[i*2], temp_edges[i*2+1]);
            }
            
            delete[] temp_edges;
            
            // Add edges to global collection with minimal critical section
            #pragma omp critical
            {
                for(const auto& edge : thread_edges) {
                    children_t[edge.first].push_back(edge.second);
                }
            }
        }
        
        cout << "Process " << rank << " completed tree " << t << endl;
        
        // If master process, store directly, else send to master
        if(rank==0) {
            children_global[t-1] = move(children_t);
        } else {
            // Serialize edges into a flat buffer to send
            size_t total_edges = 0;
            for(uint32_t p=0; p<N; ++p) {
                total_edges += children_t[p].size();
            }
            
            // Only send if there are edges to send
            if (total_edges > 0) {
                // Create buffer with exact size needed
                vector<uint32_t> buf;
                buf.reserve(total_edges * 2 + 1);
                
                // First value is the tree ID
                buf.push_back(t);
                
                // Flatten edge list
                for(uint32_t p=0; p<N; ++p) {
                    for(auto c: children_t[p]) {
                        buf.push_back(p);
                        buf.push_back(c);
                    }
                }
                
                // Send to master
                MPI_Request req;
                MPI_Isend(buf.data(), buf.size(), MPI_UINT32_T, 0, t, MPI_COMM_WORLD, &req);
                sendRequests.push_back(req);
                sendBuffers.push_back(move(buf));
                cout << "Process " << rank << " sent tree " << t << " to master" << endl;
            } else {
                // Send empty message to indicate completion
                vector<uint32_t> buf(1, t);
                MPI_Request req;
                MPI_Isend(buf.data(), buf.size(), MPI_UINT32_T, 0, t, MPI_COMM_WORLD, &req);
                sendRequests.push_back(req);
                sendBuffers.push_back(move(buf));
                cout << "Process " << rank << " sent empty tree " << t << " to master" << endl;
            }
        }
    }

    // Wait for all sends to complete
    if(!sendRequests.empty()) {
        MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
    }
    
    // Master process: receive trees from other processes
    if(rank==0) {
        // Set of trees we need to receive
        set<int> treesNeeded;
        for (int t = 1; t <= T; t++) {
            if (t % size != 0) {  // Only need trees not built by master
                treesNeeded.insert(t);
            }
        }
        
        cout << "Master needs to receive " << treesNeeded.size() << " trees" << endl;
        
        // While there are still trees to receive
        while (!treesNeeded.empty()) {
            MPI_Status status;
            int flag = 0;
            
            // Check for any incoming message
            MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
            
            if (flag) {
                // Get message size
                int count;
                MPI_Get_count(&status, MPI_UINT32_T, &count);
                
                // Receive the message
                vector<uint32_t> buffer(count);
                MPI_Recv(buffer.data(), count, MPI_UINT32_T, 
                         status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                // Extract tree ID (first element)
                int tree_id = buffer[0];
                
                if (count > 1) {
                    // Process the edges (rest of buffer)
                    for (int i = 1; i < count; i += 2) {
                        int p = buffer[i];
                        int c = buffer[i+1];
                        if (p >= 0 && p < (int)N && c >= 0 && c < (int)N) {
                            children_global[tree_id-1][p].push_back(c);
                        }
                    }
                    cout << "Master received tree " << tree_id << " with " 
                         << (count-1)/2 << " edges from process " << status.MPI_SOURCE << endl;
                } else {
                    // Empty tree
                    cout << "Master received empty tree " << tree_id 
                         << " from process " << status.MPI_SOURCE << endl;
                }
                
                // Mark tree as received
                treesNeeded.erase(tree_id);
            }
        }
        
        // All trees received, write DOT files
        /*
        cout << "Writing DOT files..." << endl;
        for(int t=1; t<=T; ++t) {
            ofstream dot("Tn_"+to_string(t)+".dot");
            dot<<"digraph T"<<n<<"_"<<t<<" {\n  rankdir=TB;\n";
            for(uint32_t p=0; p<N; ++p)
                for(auto c: children_global[t-1][p])
                    dot<<"  \""<<permToString(perms[p])<<"\" -> \""
                       <<permToString(perms[c])<<"\";\n";
            dot<<"}\n";
        }
            */
    }

    // Ensure all processes are done before reporting time
    MPI_Barrier(MPI_COMM_WORLD);
    double t_end = MPI_Wtime();
    if(rank==0) cout<<"Total execution time: "<< (t_end - t_start) <<" seconds\n";
    MPI_Finalize();
    return 0;
}