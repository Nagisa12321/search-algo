#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <queue>
#include <string>
#include <sys/time.h>
#include <unordered_map>
#include <utility>
#include <vector>
#include <omp.h>

using namespace std;

#define DEBUG_MODE true
#define PRINT_MATRIX false
#define TIMER true
#define PRINT_FILE true
#define TIMER_FILE false

void __abort(const string &__message);
void __dbg(const string &__message);
void __info(const string &__message);
double __distance(const vector<unordered_map<int, double>> &__g, int __lhs,
                  int __rhs);
void init(vector<unordered_map<int, double>> &__g, const string &__file_path);
void show_graph(const vector<unordered_map<int, double>> &__g);
void floyd_sp(const vector<unordered_map<int, double>> &__g,
              double **__dist_to);
void show_result(double **__dist_to);
void print_file(double **__dist_to);
void timer_start();
void timer_end();
void timer_print(const string &name);

int __nv; // number of vertex.
struct timeval __ts, __te;
const int __inf = 0xffff;
const string __out_name("floyd_omp.txt");
string __input_file;

int main(int __argc, char **__argv) {
  if (__argc < 2) {
    __abort("Use ./floyd_omp [input file]");
  }
  //
  // Get the Input file
  //
  __input_file = __argv[1];

  //
  // The graph of the input file.
  // __graph[i][j] means the distance vetex i to vetex j.
  //
  vector<unordered_map<int, double>> __graph;

  //
  // The short path to every node.
  //
  double **__dist_to;
  //
  // Init the graph
  //
  init(__graph, __input_file);
  //
  // Show the graph
  //
  show_graph(__graph);
  //
  // Run the algo
  //
  __dist_to = new double *[__nv];
  for (int i = 0; i < __nv; ++i)
    __dist_to[i] = new double[__nv];
  timer_start();
  floyd_sp(__graph, __dist_to);
  timer_end();
  //
  // Show the result
  //
  show_result(__dist_to);
  //
  // Show the time spend
  //
  timer_print("flyod");
  //
  // Print to the file
  //
  print_file(__dist_to);
  for (int i = 0; i < __nv; ++i)
    delete[] __dist_to[i];
  delete[] __dist_to;
}

//
// abort with a message
// will print the error message and exit this process.
//
void __abort(const string &__message) {
  cout << " >> Exit with error: " << __message << endl;
  exit(1);
}

void __dbg(const string &__message) {
  if (DEBUG_MODE)
    cout << " (debug) " << __message << endl;
}

void __info(const string &__message) {
  cout << " (info) " << __message << endl;
}

//
// init the graph use the input file
// and store the number of vetexes in __nv
//
void init(vector<unordered_map<int, double>> &__g, const string &__file_path) {
  fstream __fis(__file_path);
  int __e;
  int __offset;
  int __from, __to;
  double __distance;

  if (!__fis.is_open()) {
    __abort("Can not open the file.");
  }

  __info("start to read the file...");
  __fis >> __nv >> __e;

  __dbg("start to create the vector...");
  __g = vector<unordered_map<int, double>>(__nv);
  __dbg("ok to create vector");
  for (__offset = 0; __offset < __nv; ++__offset) {
    __g[__offset][__offset] = 0;
  }

  for (__offset = 0; __offset < __e; ++__offset) {
    __fis >> __from >> __to >> __distance;
    __g[__from][__to] = __distance;
  }

  __info("Ok to read the file.");
}

//
// Show the __graph
//
void show_graph(const vector<unordered_map<int, double>> &__g) {
  if (!PRINT_MATRIX)
    return;
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
      if (!__g[__i].count(__j)) {
        cout << "\tInf";
      } else {
        cout << "\t" << __g[__i].at(__j);
      }
    }
    cout << "\n";
  }
}

//
// Get the distance from __lhs to __rhs
// the map if __lhs is not exist __rhs, then return __inf.
//
double __distance(const vector<unordered_map<int, double>> &__g, int __lhs,
                  int __rhs) {
  if (!__g[__lhs].count(__rhs)) {
    return __inf;
  } else {
    return __g[__lhs].at(__rhs);
  }
}

//
// use the bellman-ford algo
//
void floyd_sp(const vector<unordered_map<int, double>> &__g,
              double **__dist_to) {
  __info("flyod start now!");
#pragma omp parallel num_threads(12)
  {
    int my_id = omp_get_thread_num();
    int nth = omp_get_num_threads();
    int my_first = (my_id * __nv) / nth;
    int my_last = ((my_id + 1) * __nv) / nth - 1;
  // 
  // Init the dist to
  //
  for (int i = my_first; i <= my_last; ++i) {
    for (int v = 0; v < __nv; ++v) {
      __dist_to[i][v] = __inf;
    }
  }
  // 
  // v to v is zero!
  //
  for (int i = my_first; i < my_last; ++i) {
    __dist_to[i][i] = 0;
  }
  //
  // set the weight
  //
  for (int __from = my_first; __from <= my_last; ++__from) {
    for (const pair<int, double> &__entry : __g[__from]) {
			int __to = __entry.first;
			double __weight = __entry.second;
			__dist_to[__from][__to] = __weight;
		}
  }
	// 1 let dist be a |V| × |V| array of minimum distances initialized to ∞ (infinity)
	// 2 for each vertex v
	// 3    dist[v][v] ← 0
	// 4 for each edge (u,v)
	// 5    dist[u][v] ← w(u,v)  // the weight of the edge (u,v)
	// 6 for k from 1 to |V|
	// 7    for i from 1 to |V|
	// 8       for j from 1 to |V|
	// 9          if dist[i][j] > dist[i][k] + dist[k][j] 
	// 10             dist[i][j] ← dist[i][k] + dist[k][j]
	// 11         end if
// #pragma omp barrier
	for (int k = 0; k < __nv; ++k) {
		for (int i = 0; i < __nv; ++i) {
#pragma omp barrier
			for (int j = my_first; j <= my_last; ++j) {
				if (__dist_to[i][j] > __dist_to[i][k] + __dist_to[k][j])
					{ __dist_to[i][j] = __dist_to[i][k] + __dist_to[k][j]; } 
			}
		}
	}

  }
  __info("flyod end now.");
}

void show_result(double **__dist_to) {
  if (!PRINT_MATRIX)
    return;
  int __i, __j;
  cout << "\n";
  cout << "  Minimum distances from node 0:\n";
  cout << "\n";
  for (__i = 0; __i < __nv; __i++) {
    for (__j = 0; __j < __nv; ++__j)
      cout << "  " << __i << " -> " << __j << ": " 
					 << ((__dist_to[__i][__j] == __inf) ?  "INF" : to_string( __dist_to[__i][__j]))
           << endl;
  }
}

void timer_start() {
  if (!TIMER)
    return;
  gettimeofday(&__ts, nullptr);
}
void timer_end() {
  if (!TIMER)
    return;
  gettimeofday(&__te, nullptr);
}
void timer_print(const string &name) {
  if (!TIMER)
    return;
  clock_t __spend =
      (__te.tv_sec - __ts.tv_sec) * 1000 + (__te.tv_usec - __ts.tv_usec) / 1000;
  cout << name << "'s time is " << __spend << " ms. " << endl;

  if (TIMER_FILE) {
    ofstream __fis("time.txt", ios_base::app);
    __fis << __input_file << "\t" << __spend << " ms. "
          << "\t"
          << "flyod\t" << endl;
  }
}

void print_file(double **__dist_to) {
  if (!PRINT_FILE)
    return;
  // output to the file
  std::ofstream fos(__out_name);
  for (int i = 0; i < __nv; ++i) {
    for (int j = 0; j < __nv; ++j) {
      fos << i << " -> " << j << ": " << __dist_to[i][j] << endl;
    }
  }
}
