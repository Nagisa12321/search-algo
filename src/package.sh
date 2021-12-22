#!/bin/bash

rm -rf ~/Desktop/package 2>/dev/null
mkdir ~/Desktop/package
mips64el-linux-g++ make_ewd.cpp -static -o ~/Desktop/package/make_graph
mips64el-linux-g++ bellman_omp_v2.cpp -fopenmp -static -o ~/Desktop/package/bellman_omp -g
mips64el-linux-g++ bellman.cpp -fopenmp -static -o ~/Desktop/package/bellman -g
cp ./running.sh ~/Desktop/package/