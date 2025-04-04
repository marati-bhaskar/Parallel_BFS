#include <cstdint>
#include <vector>
#include <fstream>
#include <queue>
#include <limits>
#include <iostream>
#include <omp.h>
#include <cstdlib>
#include <cstring>
#include <string>
#include <chrono>
#include <sstream>
#include <iomanip>

using vidType = uint32_t;
using eidType = uint64_t;
using weight_type = uint32_t;

struct Graph {
  eidType* rowptr;
  vidType* col;
  uint64_t N, M;
  vidType* curr_frontier;
  vidType* next_frontier;
  bool is_small_diameter;
  bool* curr_frontier_bm;
  bool* next_frontier_bm;
  int max_threads;
  vidType** local_frontier;
  vidType mdv;

  Graph() : rowptr(nullptr), col(nullptr), N(0), M(0), mdv(0) {}

  ~Graph() {
    // Free allocated memory
    free(curr_frontier);
    free(next_frontier);

    #pragma omp parallel for
    for (int i = 0; i < max_threads; ++i) {
      free(local_frontier[i]);
    }
    delete[] local_frontier;

    if(is_small_diameter){
      free(curr_frontier_bm);
      free(next_frontier_bm);
    }
    delete[] rowptr;  
    delete[] col; 
  }

  void alloc_ds()
  {
    // Allocate memory for global frontier arrays used in top-down approach
    curr_frontier = (vidType*)aligned_alloc(64, N * sizeof(vidType));
    next_frontier = (vidType*)aligned_alloc(64, N * sizeof(vidType));

    //if average degree is high, it will be more likely that the graph has small diameter
    is_small_diameter = (M/N > 7);

    //Allocate memory for frontier bitmaps used in bottom-up approach
    if(is_small_diameter){
      curr_frontier_bm = (bool*)aligned_alloc(64, N * sizeof(bool));
      next_frontier_bm = (bool*)aligned_alloc(64, N * sizeof(bool));
    }

    //get maximum no. of threads
    max_threads = omp_get_max_threads();

    // Allocate memory for local frontiers for each thread
    local_frontier = new vidType*[max_threads];
    #pragma omp parallel for
    for (int i = 0; i < max_threads; ++i) {
      local_frontier[i] = (vidType*)aligned_alloc(64, N * sizeof(vidType));
    }
  }

  bool read_binary(const char* filename) 
  {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
      std::cerr << "Could not open file: " << filename << std::endl;
      return false;
    }
    file.read((char*)&N, sizeof(N));
    file.read((char*)&M, sizeof(M));

    alloc_ds();

    rowptr = new eidType[N + 1];
    col = new vidType[M];
    file.read((char*)rowptr, sizeof(eidType) * (N + 1));
    file.read((char*)col, sizeof(vidType) * M);
    if (!file) {
      std::cerr << "Error reading graph data" << std::endl;
      delete[] rowptr;
      delete[] col;
      rowptr = nullptr;
      col = nullptr;
      return false;
    }
    mdv = 0;
    return true;
  }

  bool read_matrix_market(const char* filename) 
  {
    std::ifstream file(filename);
    if (!file) {
      std::cerr << "Could not open file: " << filename << std::endl;
      return false;
    }

    std::string line;
    while (std::getline(file, line)) {
      if (line[0] == '%') continue; // Skip comments
      break; // Found the first non-comment line
    }

    std::istringstream iss(line);
    iss >> N >> N >> M; // Read the number of rows (N), columns (N), and non-zeros (M)

    alloc_ds();

    rowptr = new eidType[N + 1]();
    col = new vidType[M];

    std::vector<std::pair<vidType, vidType>> edges;
    edges.reserve(M);

    eidType edge_index = 0;

    while (std::getline(file, line)) {
      vidType u, v;
      std::istringstream edge_iss(line);
      edge_iss >> u >> v;
      edge_iss.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

      // Adjust from 1-based indexing to 0-based indexing
      u--;
      v--;

      edges.emplace_back(u, v);
      rowptr[u + 1]++;
      col[edge_index++] = v;
    }

    for (vidType i = 1; i <= N; i++) {
      rowptr[i] += rowptr[i - 1];
    }

    uint64_t max_degree = 0;
    vidType max_degree_vertex = 0;
    for (vidType i = 1; i <= N; i++) {
      uint64_t deg = rowptr[i] - rowptr[i-1];
      if (max_degree < deg) {
        max_degree = deg;
        max_degree_vertex = i-1;
      }
    }
    mdv = max_degree_vertex;
    std::cout << "no. of neighbors of source vertex are: " << rowptr[mdv + 1] - rowptr[mdv] << std::endl;

    return true;
  }



  void bfs_serial(vidType source, weight_type* distances) const {
    for (uint64_t i = 0; i < N; i++) {
      distances[i] = std::numeric_limits<weight_type>::max();
    }
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<vidType> frontier;
    distances[source] = 0;
    frontier.push_back(source);
    while (!frontier.empty()) {
      std::vector<vidType> next_frontier;
      for (const auto& src : frontier) {
        for (uint64_t i = rowptr[src]; i < rowptr[src + 1]; i++) {
          vidType dst = col[i];
          if (distances[src] + 1 < distances[dst]) {
            distances[dst] = distances[src] + 1;
            next_frontier.push_back(dst);
          }
        }
      }
      frontier = std::move(next_frontier);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Serial BFS Time: " << diff.count() * 1000 << " ms" << std::endl;
  }

  void bfs_parallel(vidType source, weight_type* distances) 
  {
    for (uint64_t i = 0; i < N; i++) {
      distances[i] = std::numeric_limits<weight_type>::max();
    }
    auto start = std::chrono::high_resolution_clock::now();
    uint32_t curr_frontier_size = 1;
    curr_frontier[0] = source;
    distances[source] = 0;
    uint32_t next_frontier_size = 0;

    //top-down approach
    while (curr_frontier_size != 0) {
      //start the parallel region
      #pragma omp parallel
      {                 
        uint32_t local_frontier_size = 0;
        uint32_t frontier_starting_index;
        int tid = omp_get_thread_num();

        //split work among the threads
        #pragma omp for schedule(static) nowait 
        for (uint32_t idx = 0; idx < curr_frontier_size; ++idx) {
          vidType src = curr_frontier[idx];
          for (uint32_t i = rowptr[src]; i < rowptr[src + 1]; ++i) {
            vidType dst = col[i];
            if (distances[src] + 1 < distances[dst]) {
              if (__sync_bool_compare_and_swap(&distances[dst], std::numeric_limits<weight_type>::max(), distances[src] + 1)) {
                local_frontier[tid][local_frontier_size++] = dst;
              }
              
            }
          }
        }
        //collect the local frontiers into global frontier
        if(local_frontier_size != 0) {
          //atomically get the index in the (global)  next_frontier  from where local frontier elements of each element can be copied 
          frontier_starting_index = __atomic_fetch_add(&next_frontier_size, local_frontier_size, __ATOMIC_SEQ_CST);
          // copy the elements in local frontier to global frontier
          for (uint32_t i = 0; i < local_frontier_size; ++i) {
            next_frontier[frontier_starting_index + i] = local_frontier[tid][i];
          }
        }
        
        local_frontier_size = 0;
      }
      //swap current and next (global) frontiers
      std::swap(curr_frontier, next_frontier);
      curr_frontier_size = next_frontier_size;
      next_frontier_size = 0;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Parallel BFS Time: " << diff.count() * 1000 << " ms" << std::endl;
    std::cout << "Edges processed per second: " 
              << std::fixed << std::setprecision(0) 
              << M / (diff.count()) 
              << std::endl;  
  }

  bool verify_bfs(vidType source, bool verify) {
    std::vector<weight_type> parallel_distances(N);
    bfs_parallel(source, parallel_distances.data());
    if (!verify) {
      return true;
    }
    std::vector<weight_type> serial_distances(N);
    bfs_serial(source, serial_distances.data());
    for (uint64_t i = 0; i < N; i++) {
      if (serial_distances[i] != parallel_distances[i]) {
        std::cout << "Mismatch at vertex " << i << ": serial=" << serial_distances[i]
                  << " parallel=" << parallel_distances[i] << std::endl;
        return false;
      }
    }
    return true;
  }
};

bool ends_with(const std::string& str, const std::string& suffix) 
{
  return str.size() >= suffix.size() &&
         str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

int main(int argc, char* argv[]) 
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <graph_file.[pbin|mtx]>" << " bool: verify" << std::endl;
    return 1;
  }

  std::string filename = argv[1];
  bool verify = std::stoi(argv[2]);

  Graph g;

  if (ends_with(filename, ".pbin")) {
    if (!g.read_binary(filename.c_str())) {
      std::cerr << "Failed to read binary graph (.pbin)" << std::endl;
      return 1;
    }
  } else if (ends_with(filename, ".mtx")) {
    if (!g.read_matrix_market(filename.c_str())) {
      std::cerr << "Failed to read Matrix Market graph (.mtx)" << std::endl;
      return 1;
    }
  } else {
    std::cerr << "Unsupported file format. Use .pbin or .mtx" << std::endl;
    return 1;
  }

  bool is_correct = g.verify_bfs(g.mdv, verify);

  return 0;
}