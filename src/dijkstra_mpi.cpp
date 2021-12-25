#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <mpi.h>
#include <mpio.h>
#include <string.h>
#include <unordered_map>
#include <vector>

using namespace std;

// inf
const int inf = 0xffff;
const int lock_tag = 0xeeff;

double distance(const vector<unordered_map<int, double>> &ohd, int i, int j);
void cycle_criticalregion(int comm_sz, int my_rank, function<void()> fn);

int main(int argc, char **argv) {
  if (argc != 2) {
    cout << "please input ./dijkstra_mpi [input file name]" << endl;
    exit(1);
  }

  string input_file = argv[1];
  int comm_sz; // process num
  int my_rank; // process id
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  //
  // Open File
  //
  const char *file_connected = "./connected.txt";
  const char *file_dist_to = "./dist_to.txt";
  const char *file_variables = "./dist_to.txt";
  // Create two files.
  // connected.txt and dist_to.txt
  MPI_File fh_connected, fh_dist_to, fh_variables;

  MPI_File_open(MPI_COMM_WORLD, file_connected,
                MPI_MODE_RDWR | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_CREATE,
                MPI_INFO_NULL, &fh_connected);
  MPI_File_open(MPI_COMM_WORLD, file_dist_to,
                MPI_MODE_RDWR | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_CREATE,
                MPI_INFO_NULL, &fh_dist_to);
  MPI_File_open(MPI_COMM_WORLD, file_variables,
                MPI_MODE_RDWR | MPI_MODE_DELETE_ON_CLOSE | MPI_MODE_CREATE,
                MPI_INFO_NULL, &fh_variables);

  if (my_rank != 0) {
    int nv;
    MPI_Recv(&nv, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int my_first = ((my_rank - 1) * nv) / (comm_sz - 1);
    int my_last = ((my_rank)*nv) / (comm_sz - 1) - 1;

    vector<unordered_map<int, double>> my_graph(nv);
    //
    // Recvice my data...
    //
    for (int i = my_first; i <= my_last; ++i) {
      double *tmp = new double[nv];
      MPI_Recv(tmp, nv, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (int v = 0; v < nv; ++v) {
        if (tmp[v] != inf) {
          my_graph[i][v] = tmp[v];
        }
      }
      delete[] tmp;
    }

    //
    // Show my graph, and some message
    //
    cycle_criticalregion(comm_sz, my_rank, [&] {
      cout << " >>>>>>> This is proc " << my_rank;
      cout << ", my_first is" << my_first << ", my_last is " << my_last << ". ";
      int i, j;
      cout << "  matrix:\n";
      cout << "\t";
      for (i = 0; i < nv; i++)
        cout << "\t[" << i << "]";
      cout << endl;
      for (i = my_first; i <= my_last; i++) {
        cout << "\t[" << i << "]";
        for (j = 0; j < nv; j++) {
          if (!my_graph[i].count(j)) {
            cout << "\tInf";
          } else {
            cout << "\t" << my_graph[i][j];
          }
        }
        cout << "\n";
      }
      cout << "\n" << endl;
    });

    //
    // Start the real algo
    //
    double *dist_to = new double[nv];
    int *connected = new int[nv];
    int mv;
    double md;

    for (int my_step = 1; my_step < nv; ++my_step) {
      // Init the md and mv. 
      if (my_rank == 1) {
        md = inf;
        mv = -1;
        MPI_File_write_at(fh_variables, 0, &md, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_write_at(fh_variables, 1, &mv, 1, MPI_INT, MPI_STATUS_IGNORE);
      }
      
      double my_md = inf;
      int my_mv = -1;
      //
      //  Each process finds the nearest unconnected node in its part of the graph.
      //  Some threads might have no unconnected nodes left.
      //  get dist_to and connected. 
      //
      MPI_File_read(fh_connected, connected, nv, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_read(fh_dist_to, dist_to, nv, MPI_DOUBLE, MPI_STATUS_IGNORE);
      for (int i = my_first; i <= my_last; ++i) {
        if (!connected[i] && dist_to[i] < my_md) {
          my_md = dist_to[i];
          my_mv = i;
        }
      }

      //
      // Choose the best mv, md
      //
      cycle_criticalregion(comm_sz, my_rank, [&]{
        //
        // Read the mv and md. 
        //
        MPI_File_read_at(fh_variables, 0, &mv, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_read_at(fh_variables, 1, &md, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);

        // 
        // If I have greate mv and md 
        // renew it!
        // 
        if (my_md < md) {
          MPI_File_write_at(fh_variables, 0, &md, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
          MPI_File_write_at(fh_variables, 1, &mv, 1, MPI_INT, MPI_STATUS_IGNORE);
        }
      });

      //
      // make sure that md is the small path to a new node.
      //
      MPI_Barrier(MPI_COMM_WORLD);
      //
      // Read the mv and md. 
      // mv and md are now the latest and correct values
      //
      MPI_File_read_at(fh_variables, 0, &mv, 1, MPI_INT, MPI_STATUS_IGNORE);
      MPI_File_read_at(fh_variables, 1, &md, 1, MPI_DOUBLE, MPI_STATUS_IGNORE);
      if (mv != -1) {
        for (int i = my_first; i <= my_last; ++i) {
          if (!connected[i]) {
            //
            // Why we can do this ? 
            // Because mv is now in the tree. 
            // So dist_to[mv] is the correct value!
            // 
            if (dist_to[i] > dist_to[mv] + distance(my_graph, mv, i)) {
              dist_to[i] = dist_to[mv] + distance(my_graph, mv, i);
            }
          }
        }
      }

      // update the dist_to to the file . 

      //
      // Make sure that every thread have update the dist_to.
      // 
      MPI_Barrier(MPI_COMM_WORLD);
    }

    delete[] connected;
    delete[] dist_to;
  } else {
    //
    // Read the input file
    // and send it to every process.
    //
    int nv;    // number of vertexes
    int edges; // number of edges
    cout << "process 0 start to read file." << endl;

    ifstream fis(input_file);
    fis >> nv >> edges;

    cout << "   Number of vertex: " << nv << endl;
    cout << "   Number of edges: " << edges << endl;

    //
    // Send the number of vertex to every process
    //
    for (int rank = 1; rank < comm_sz; ++rank) {
      MPI_Send(&nv, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
    }

    cout << "Now read the file .... ";
    vector<unordered_map<int, double>> graph(nv);
    for (int i = 0; i < nv; ++i) {
      graph[i][i] = 0;
    }

    int from, to;
    double d;
    for (int i = 0; i < edges; ++i) {
      fis >> from >> to >> d;
      graph[from][to] = d;
    }
    cout << "ok" << endl;

    //
    // Print the matrix
    //
    int i, j;
    cout << "\n";
    cout << "  Distance matrix:\n";
    cout << "\n";
    cout << "\t";
    for (i = 0; i < nv; i++)
      cout << "\t[" << i << "]";
    cout << endl;
    for (i = 0; i < nv; i++) {
      cout << "\t[" << i << "]";
      for (j = 0; j < nv; j++) {
        if (!graph[i].count(j)) {
          cout << "\tInf";
        } else {
          cout << "\t" << graph[i][j];
        }
      }
      cout << "\n";
    }

    cout << "Now start to send the data to the process. " << endl;
    //
    // send the different to different process.
    //
    for (int rank = 1; rank < comm_sz; ++rank) {
      int first = ((rank - 1) * nv) / (comm_sz - 1);
      int last = ((rank)*nv) / (comm_sz - 1) - 1;
			
      for (int i = first; i <= last; ++i) {
        double *tmp = new double[nv];
        for (int v = 0; v < nv; ++v) {
          tmp[v] = distance(graph, i, v);
        }
        // send every line of graph to the proc
        MPI_Send(tmp, nv, MPI_DOUBLE, rank, 0, MPI_COMM_WORLD);
        delete[] tmp;
      }
    }

    //
    // Setup the connected and dist_to. 
    //
    int *connected = new int[nv];
    double *dist_to = new double[nv];
    memset(connected, 0, nv * sizeof(int));
    memset(dist_to, 0, nv * sizeof(double));

    // 
    // setup the connected and dist_to
    // 
    connected[0] = true;
    for (int i = 0; i < nv; ++i) 
      dist_to[i] = distance(graph, 0, i);
    //
    // Write the file
    // 1) init the connected/dist_to
    // 2) init the md, mv
    // 
    MPI_File_write(fh_connected, connected, nv, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_write(fh_dist_to, dist_to, nv, MPI_DOUBLE, MPI_STATUS_IGNORE);

    delete[] connected;
    delete[] dist_to;

    //
    // Make a criticalregion
    //
    MPI_Send(&nv, 1, MPI_INT, 1, lock_tag, MPI_COMM_WORLD);
  }

  //
  // Close the opened file.
  //
  MPI_File_close(&fh_connected);
  MPI_File_close(&fh_dist_to);

  MPI_Finalize();
}

double distance(const vector<unordered_map<int, double>> &ohd, int i, int j) {
  if (!ohd[i].count(j))
    return inf;
  else
    return ohd[i].at(j);
}

//
// !Warnning:
// !It's better not to use it.
// !Because it's too slow...
//
void cycle_criticalregion(int comm_sz, int my_rank, function<void()> fn) {
  int value;
  MPI_Recv(&value, 1, MPI_INT, my_rank - 1, lock_tag, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  fn();
  if (my_rank < comm_sz - 1)
    MPI_Send(&value, 1, MPI_INT, my_rank + 1, lock_tag, MPI_COMM_WORLD);
}