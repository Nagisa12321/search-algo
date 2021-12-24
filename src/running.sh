#!/bin/bash
for i in {1..100}
do
  v=$(( $i * 100 ))
  e=$(( $v * 10 ))
  echo ">> vertexes=$v, edges=$e"
  ./make_graph $v $e

  echo "run bellman now..."
  ./_r_bellman $v >/dev/null 2>&1
  echo "run bellman_omp now..."
  ./_r_bellman_omp $v >/dev/null 2>&1
  echo "run dijkstra and dijkstra omp now..."
  ./_r_dijkstra_omp $v >/dev/null 2>&1

  rm "$v"
done