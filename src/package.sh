#!/bin/bash

rm -rf ~/Desktop/package 2>/dev/null
mkdir ~/Desktop/package
mips64el-linux-g++ make_ewd.cpp -static -o ~/Desktop/package/make_graph
mips64el-linux-g++ bellman_omp.cpp -fopenmp -static -o ~/Desktop/package/_r_bellman_omp
mips64el-linux-g++ dijkstra_openmp_adj.cpp -fopenmp -o ~/Desktop/package/_r_dijkstra_omp -static
mips64el-linux-g++ bellman.cpp -static -o ~/Desktop/package/_r_bellman
cp ./running.sh ~/Desktop/package/