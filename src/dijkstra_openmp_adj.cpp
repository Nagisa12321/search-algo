#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <ctime>
#include <omp.h>
#include <sys/time.h>
#include <unordered_map>

using namespace std;

// config
#define PRINT_DISTANCE_MATRIX false
#define PRINT_RESULT          false
#define PRINT_FILE            false
#define CHECK_TIME            1
#define DEBUG                 false
#define TIMER_FILE 		true

string __input_file_name;
const string __output_file_name("__result_omp.txt");
const int __i4_huge = 0xffff;
int __nv;


void __init(vector<unordered_map<int, double>> &__ohd);
double *__dijkstra_distance(const vector<unordered_map<int, double>> &__ohd);
double *__dijkstra_distance_single(const vector<unordered_map<int, double>> &__ohd);
void __print_matrix(const vector<unordered_map<int, double>> &__ohd);
void __print_result(double *__dist_to);
void __print_file(double *__dist_to);
bool __compare_result(const double *__dt1, const double *__dt2);


int main(int argc, char **argv) {
  __input_file_name = argv[1];
  double *__dist_to_multi, *__dist_to;
  vector<unordered_map<int, double>> __ohd;
#if CHECK_TIME
  struct timeval __ts_multi, __te_multi, __ts, __te;
#endif

  // init the program data;
  __init(__ohd);

  // print the distance matrix
  if (PRINT_DISTANCE_MATRIX)
    __print_matrix(__ohd);


#if CHECK_TIME
  gettimeofday(&__ts_multi, nullptr);
#endif
  __dist_to_multi = __dijkstra_distance(__ohd);
#if CHECK_TIME
  gettimeofday(&__te_multi, nullptr);
  clock_t __spend_multi = (__te_multi.tv_sec - __ts_multi.tv_sec) * 1000 
      + (__te_multi.tv_usec - __ts_multi.tv_usec) / 1000;
#endif
  
#if CHECK_TIME
  gettimeofday(&__ts, nullptr);
#endif
  __dist_to = __dijkstra_distance_single(__ohd);
#if CHECK_TIME
  gettimeofday(&__te, nullptr);
  clock_t __spend = (__te.tv_sec - __ts.tv_sec) * 1000 
      + (__te.tv_usec - __ts.tv_usec) / 1000;
#endif  

  // print result
  if (PRINT_RESULT)
    __print_result(__dist_to_multi);

#if CHECK_TIME
  cout << "  multi thread spend time: " << __spend_multi << " ms." << endl;
  cout << "  single thread spend time: " << __spend << " ms." << endl;
#endif
  cout << "   Start to compare the result ... " << endl;
  // bool __res = __compare_result(__dist_to_multi, __dist_to);
  // if (!__res) {
  //   cout << " ---- Failed ----" << endl;
  // } else {
  //   cout << " ---- Successful ----" << endl;
  // }


  // print to file
  if (PRINT_FILE)
    __print_file(__dist_to_multi);

  if (TIMER_FILE) {
    ofstream __fis("time.txt", ios_base::app);
    __fis << __input_file_name << "\t" << __spend << " ms. " << "\t" << "dijkstra\t" << endl;
    __fis << __input_file_name << "\t" << __spend_multi << " ms. " << "\t" << "dijkstra_omp\t" << endl;
  }

  delete[] __dist_to_multi;
  delete[] __dist_to;
}

void __init(vector<unordered_map<int, double>> &__ohd) {
  int i;
  int j;
  int edges;
  int from, to;
  double distance;
  ifstream fis(__input_file_name);
  if (!fis.is_open()) {
    std::cout << "failed to open " << __input_file_name << '\n';
    exit(1);
  }

  fis >> __nv;
  fis >> edges;
  cout << " Start to read the input file ..... " << endl;
  __ohd = vector<unordered_map<int, double>>(__nv);
  for (int i = 0; i < __nv; ++i)
    { __ohd[i][i] = 0; }
  for (i = 0; i < edges; ++i) {
    fis >> from >> to;
    fis >> distance;

    __ohd[from][to] = distance;
  }
  cout << "  Ok to read the input file .... " << endl;
}

static double __get_val(const vector<unordered_map<int, double>> &ohd, int i, int j) {
  if (!ohd[i].count(j)) return __i4_huge;
  else return ohd[i].at(j);
}

double *__dijkstra_distance(const vector<unordered_map<int, double>> &ohd) {
  vector<bool> connected(__nv, false);
  double *dist_to;
  int step;
  int mv; // the index of the nearest unconnected node.
  double md; // the distance from node 0 to the nearest unconnected node.
  int i;
  connected[0] = true;
  dist_to = new double[__nv];

  // init dist_to to one_step distance;
  for (i = 0; i < __nv; ++i)
    dist_to[i] = __get_val(ohd, 0, i);

#pragma omp parallel num_threads(12) \
  shared (connected, md, dist_to, mv, ohd)
  {
    int i;
    int my_step;
    double my_md; 
    int my_mv;
    int my_id = omp_get_thread_num();
    int nth = omp_get_num_threads();
    int my_first = (my_id * __nv) / nth;
    int my_last = ((my_id + 1) * __nv) / nth - 1;


#pragma omp critical
    { 
      cout << "  t" << my_id << ": " << my_first << " - " << my_last << endl;
    }
# pragma omp barrier

    for (my_step = 1; my_step < __nv; ++my_step) {
// 
// Only one thread have to do this.  
//
#pragma omp single
      {
        md = __i4_huge;
        mv = -1;
      }
//
//  Each thread finds the nearest unconnected node in its part of the graph.
//  Some threads might have no unconnected nodes left.
//
      my_md = __i4_huge;
      my_mv = -1;
      for (i = my_first; i <= my_last; ++i) {
        if (!connected[i] && dist_to[i] < my_md) {
          my_md = dist_to[i];
          my_mv = i;
        }
      }
//
// Choose the best mv, md
//
#pragma omp critical
      {
        if (my_md < md) {
          md = my_md;
          mv = my_mv;
        }
      }
//
// make sure that md is the small path to a new node.
//
# pragma omp barrier
//
// OpenMP does not like to BREAK out of a parallel region.
//
# pragma omp single 
      {
        if (mv != - 1) {
          connected[mv] = true;
          if (DEBUG)
            cout << "  t" << my_id
                << ": Connecting node " << mv << "\n";;
        }
      }  
//
// This barrier is make sure the connect is upload.
//
# pragma omp barrier

//
// now we are going to update the path_to vector (relax the node. )
//
      if (mv != -1) {
        for (i = my_first; i <= my_last; ++i) {
          if (!connected[i]) {
            //
            // Why we can do this ? 
            // Because mv is now in the tree. 
            // So dist_to[mv] is the correct value!
            // 
            if (dist_to[i] > dist_to[mv] + __get_val(ohd, mv, i)) {
              dist_to[i] = dist_to[mv] + __get_val(ohd, mv, i);
            }
          }
        }
      }
//
// Make sure that every thread have update the dist_to.
// 
# pragma omp barrier
    } // for
#pragma omp critical 
    {
      cout << "  t" << my_id << " end" << endl;
    }
  } // omp parallel

  return dist_to;
}

void __print_result(double *__dist_to) {
  int __i;
  // print results.
  cout << "\n";
  cout << "  Minimum distances from node 0:\n";
  cout << "\n";
  for (__i = 0; __i < __nv; __i++) {
    cout << "  " << setw(2) << __i << "  " << setw(2) << __dist_to[__i] << "\n";
  }
}

void __print_file(double *__dist_to) {
  // output to the file
  std::ofstream fos(__output_file_name);
  for (int i = 0; i < __nv; ++i) {
    fos << i << ": " << __dist_to[i] << endl;
  }
}

void __print_matrix(const vector<unordered_map<int, double>> &__ohd) {
  int __i, __j;
  cout << "\n";
  cout << "  Distance matrix:\n";
  cout << "\n";
  cout << "\t";
  for (__i = 0; __i < __nv; __i++)
    cout << "\t[" << __i << "]";
  cout << endl;
  for (__i = 0; __i < __nv; __i++) {
    cout << "\t[" << __i << "]";
    for (__j = 0; __j < __nv; __j++) {
      if (!__ohd[__i].count(__j)) {
        cout << "\tInf";
      } else {
        cout << "\t" << __ohd[__i].at(__j);
      }
    }
    cout << "\n";
  }
}

double *__dijkstra_distance_single(const vector<unordered_map<int, double>> &ohd) {
  vector<bool> connected(__nv, false);
  double *dist_to;
  int i;
  int step;
  int mv; // the index of the nearest unconnected node.
  double md; // the distance from node 0 to the nearest unconnected node.

  connected[0] = true;
  dist_to = new double[__nv];

  // init dist_to to one_step distance;
  for (i = 0; i < __nv; ++i)
    dist_to[i] = __get_val(ohd, 0, i);

  // every step will find a new node connect to the tree.
  for (step = 1; step < __nv; ++step) {
    md = __i4_huge;
    mv = -1;
    for (i = 0; i < __nv; ++i) {
      if (!connected[i] && dist_to[i] < md) {
        md = dist_to[i];
        mv = i;
      }
    }

    // can not find any vertex
    if (mv == -1) {
      cout << "\n";
      cout << "DIJKSTRA_DISTANCE - Warning!\n";
      cout << "  Search terminated early.\n";
      cout << "  Graph might not be connected.\n";
      break;
    }

    // find a good vertex, so mark it connected.
    connected[mv] = true;

    // update the dist_to (relax the finded vertex)
    for (i = 0; i < __nv; ++i) {
      // the connected can be the best
      // TODO: prove here
      if (!connected[i]) {
        if (dist_to[i] > dist_to[mv] + __get_val(ohd, mv, i)) {
          dist_to[i] = dist_to[mv] + __get_val(ohd, mv, i);
        }
      }
    }
  }
  return dist_to;
}


bool __compare_result(const double *__dt1, const double *__dt2) {
  bool __res = true;
  for (int i = 0; i < __nv; ++i) {
    if (__dt1[i] != __dt2[i]) {
      cout << "  !!!Wawrnning: " << i 
           << " -> " 
           << __dt1[i] 
           << ", " 
           << __dt2[i];
      __res = false;
    }
  }
  return __res;
}