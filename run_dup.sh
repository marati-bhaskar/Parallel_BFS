#!/bin/bash
#SBATCH --job-name=bfs_dup_runtime      # Job name
#SBATCH --output=bfs_dup_output_%j.out  # Output file
#SBATCH --error=bfs_dup_error_%j.err    # Error file
#SBATCH --time=24:00:00             # Time limit
#SBATCH --partition=mi2104x         # Partition name
#SBATCH --cpus-per-task=64          # Maximum CPUs per task

rm bfs.out
g++  -fopenmp -O3 -o bfs.out bfs_dup_count.cc

# Set OpenMP environment variables for thread pinning
export OMP_NUM_THREADS=64
export OMP_PLACES=cores
export OMP_PROC_BIND=TRUE

# List of graph files
graph_files=("collaboration_network.pbin" "social_network.pbin"   "synthetic_dense.pbin" "web_graph.pbin" "kNN_graph.pbin" "road_network_1.pbin" "road_network_2.pbin" "synthetic_sparse.pbin" )
# graph_files=("social_network.pbin" "road_network_1.pbin" "collaboration_network.pbin" "synthetic_sparse.pbin" "kNN_graph.pbin" "synthetic_dense.pbin" "road_network_2.pbin" "web_graph.pbin"  "kmer_V1r.mtx" "com-Friendster.mtx" "GAP-kron.mtx"  "MOLIERE_2016.mtx" "AGATHA_2015.mtx")
# graph_files=("kmer_V1r.mtx" )

# Output CSV file
output_csv="duplicate_64t.csv"

# Loop through graph files
for graph in "${graph_files[@]}"; do
  # Append dataset name with a newline before and after it
#   echo "" >> $output_csv      # Add a newline before the dataset name
  echo "$graph" >> $output_csv
  echo "" >> $output_csv      # Add a newline after the dataset name

  # Run the BFS program
  taskset -c 0-63 ./bfs.out $graph 0
done
