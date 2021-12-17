#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

// config
#define PRINT_DISTANCE_MATRIX false
#define PRINT_RESULT true
#define PRINT_FILE true

string __input_file_name;
const string __output_file_name("__result_adj.txt");
const int __i4_huge = 0xffff;
int __nv;

void __init(vector<unordered_map<int, double>> &__ohd);
double *__dijkstra_distance(const vector<unordered_map<int, double>> &__ohd);
void __print_matrix(const vector<unordered_map<int, double>> &__ohd);
void __print_result(double *__dist_to);
void __print_file(double *__dist_to);

int main(int argc, char **argv) {
  __input_file_name = argv[1];
  double *__dist_to;
  vector<unordered_map<int, double>> __ohd;

  // init the program data;
  __init(__ohd);

  // print the distance matrix
  if (PRINT_DISTANCE_MATRIX)
    __print_matrix(__ohd);

  // use the algo
  __dist_to = __dijkstra_distance(__ohd);

  // print result
  if (PRINT_RESULT)
    __print_result(__dist_to);

  // print to file
  if (PRINT_FILE)
    __print_file(__dist_to);

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
  __ohd = vector<unordered_map<int, double>>(__nv);
  for (int i = 0; i < __nv; ++i)
    { __ohd[i][i] = 0; }
  for (i = 0; i < edges; ++i) {
    fis >> from >> to;
    fis >> distance;

    __ohd[from][to] = distance;
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

static double __get_val(const vector<unordered_map<int, double>> &ohd, int i, int j) {
  if (!ohd[i].count(j)) return __i4_huge;
  else return ohd[i].at(j);
}

double *__dijkstra_distance(const vector<unordered_map<int, double>> &ohd) {
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