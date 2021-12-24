#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <unordered_map>

using namespace std;

int main(int argc, char **argv) {
  if (argc != 2) {
    cout << "please input ./dijkstra_mpi [input file name]" << endl;
    exit(1);
  }

  string input_file = argv[1];
  int comm_sz;    // process num
  int my_rank;    // process id  
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank != 0) {

  } else {
    //
    // Read the input file 
    // and send it to every process. 
    //
    int nv;     // number of vertexes 
    int edges;  // number of edges 
    cout << "process 0 start to read file." << endl;

    ifstream fis(input_file);
    fis >> nv >> edges;

    cout << "   Number of vertex: " << nv << endl;
    cout << "   Number of edges: " << edges << endl;

    cout << "Now read the file .... ";
    vector<unordered_map<int, double>> graph = vector<unordered_map<int, double>>(nv);
    for (int i = 0; i < nv; ++i) {
      graph[i][i] = 0;
    }

    int from, to;
    double distance;
    for (int i = 0; i < edges; ++i) {
      fis >> from >> to >> distance;
      graph[from][to] = distance;
    }
    cout << "ok" << endl;

    cout << "Now start to send the data to the process. " << endl;
    //
    // send the different to different process. 
    // 
    for (int rank = 1; rank < comm_sz; ++rank) {
      int first = ((rank - 1) * nv) / (comm_sz - 1);
      int last = ((rank) * nv) / (comm_sz - 1) - 1;
      cout << "index " << first << " to " << last 
           << " will send to process " << rank << endl;
    }
  }

  MPI_Finalize();
}