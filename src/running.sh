#!/bin/bash
for i in {1..100}
do
  v=$(( $i * 1000 ))
  e=$(( $v * 2 ))
  echo ">> vertexes=$v, edges=$e"
  ./make_graph $v $e

  echo "run bellman now"
  ./bellman $v >/dev/null 2>&1
  echo "run bellman_omp now"
  ./bellman_omp $v >/dev/null 2>&1

  rm "$v"
done