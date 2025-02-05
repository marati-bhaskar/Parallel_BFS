#!/bin/bash
#SBATCH --job-name=bfs_runtime      # Job name
#SBATCH --output=bfs_output_%j.out  # Output file
#SBATCH --error=bfs_error_%j.err    # Error file
#SBATCH --time=24:00:00             # Time limit
#SBATCH --partition=mi2104x         # Partition name
#SBATCH --cpus-per-task=32          # Maximum CPUs per task


# Run different versions of BFS

# List of graph files
graph_files=("collaboration_network.pbin" "social_network.pbin" "synthetic_dense.pbin" "web_graph.pbin" "kNN_graph.pbin" "road_network_1.pbin" "road_network_2.pbin" "synthetic_sparse.pbin")
# graph_files=("delaunay_n24.mtx" "hugetrace-00020.mtx" "europe_osm.mtx" "com-Friendster.mtx" "mawi_201512020130.mtx" "kmer_V1r.mtx" "GAP-kron.mtx")
# Number of times to run each experiment
N=5   # Total runs
DISCARD=1  # Number of initial runs to discard
VALID_RUNS=$((N - DISCARD))  # Number of runs to average

# Set OpenMP environment variables for thread pinning
export OMP_NUM_THREADS=32
export OMP_PLACES=cores
export OMP_PROC_BIND=TRUE

# Output CSV files
runtime_csv="bfs_32t_v1_runtime.csv"
throughput_csv="bfs_32t_v1_throughput.csv"

# Create CSV headers
echo -n "Version," > $runtime_csv
echo -n "Version," > $throughput_csv
for graph in "${graph_files[@]}"; do
  echo -n "$graph," >> $runtime_csv
  echo -n "$graph," >> $throughput_csv
done

echo "" >> $runtime_csv
echo "" >> $throughput_csv

# Loop through BFS versions
# for version in {1..5}; do
version=1
  bfs_source="bfs_v${version}.cc"
  bfs_executable="bfs_v${version}.out"

  # Compile the BFS version
  g++ -std=c++11 -fopenmp -mavx2 -O3 -o $bfs_executable $bfs_source

  echo -n "v${version}," >> $runtime_csv
  echo -n "v${version}," >> $throughput_csv

  # Loop through graph files
  for graph in "${graph_files[@]}"; do
    runtimes=()
    throughputs=()

    # Run the experiment N times
    for ((run=1; run<=N; run++)); do
      output=$(taskset -c 0-31 ./$bfs_executable $graph 0)
      runtime=$(echo "$output" | grep "Parallel BFS Time" | awk '{print $4}')
      throughput=$(echo "$output" | grep "Edges processed per second" | awk '{print $5}')
      runtimes+=($runtime)
      throughputs+=($throughput)
    done

    # Discard the first `DISCARD` runs and calculate the average of the rest
    sum_runtime=0
    sum_throughput=0
    for ((i=DISCARD; i<N; i++)); do
      sum_runtime=$(echo "$sum_runtime + ${runtimes[i]}" | bc)
      sum_throughput=$(echo "$sum_throughput + ${throughputs[i]}" | bc)
    done
    average_runtime=$(echo "scale=2; $sum_runtime / $VALID_RUNS" | bc)
    average_throughput=$(echo "scale=2; $sum_throughput / $VALID_RUNS" | bc)

    # Append the average runtime and throughput to respective CSVs
    echo -n "$average_runtime," >> $runtime_csv
    echo -n "$average_throughput," >> $throughput_csv
  done

  echo "" >> $runtime_csv
  echo "" >> $throughput_csv

done
