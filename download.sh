#!/bin/bash

urls=(
  https://suitesparse-collection-website.herokuapp.com/MM/GenBank/kmer_V1r.tar.gz
  https://suitesparse-collection-website.herokuapp.com/MM/GAP/GAP-kron.tar.gz
  https://suitesparse-collection-website.herokuapp.com/MM/SNAP/com-Friendster.tar.gz
  https://suitesparse-collection-website.herokuapp.com/MM/Sybrandt/MOLIERE_2016.tar.gz
  
  # "https://speedcode-proto-problem-data.s3.us-east-1.amazonaws.com/ppopp-graphs-v3/south-america_wgh.pbin"
  # "https://speedcode-proto-problem-data.s3.us-east-1.amazonaws.com/ppopp-graphs-v3/soc-LiveJournal1_sym.pbin"
  # "https://speedcode-proto-problem-data.s3.us-east-1.amazonaws.com/ppopp-graphs-v3/hollywood_2009_sym.pbin"
  # "https://speedcode-proto-problem-data.s3.us-east-1.amazonaws.com/ppopp-graphs-v3/GeoLifeNoScale_5_sym_wgh.pbin"
  # "https://speedcode-proto-problem-data.s3.us-east-1.amazonaws.com/ppopp-graphs-v3/enwiki-2023_sym.pbin"
  # "https://speedcode-proto-problem-data.s3.us-east-1.amazonaws.com/ppopp-graphs-v3/north-america_wgh.pbin"
  # "https://speedcode-proto-problem-data.s3.us-east-1.amazonaws.com/ppopp-graphs-v3/er.pbin"
  # "https://speedcode-proto-problem-data.s3.us-east-1.amazonaws.com/ppopp-graphs-v3/grid_1000_10000_sym.pbin"
)

for url in "${urls[@]}"; do
  wget "$url"
done
