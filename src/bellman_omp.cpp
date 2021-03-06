#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <omp.h>
#include <ostream>
#include <queue>
#include <string>
#include <sys/time.h>
#include <unordered_map>
#include <vector>

using namespace std;

#define DEBUG_MODE 		false
#define PRINT_MATRIX	false
#define TIMER		    	true
#define PRINT_FILE		false
#define TIMER_FILE 		true

void __abort(const string &__message);
void __dbg(const string &__message);
void __info(const string &__message);
double __distance(const vector<unordered_map<int, double>> &__g, int __lhs,
                  int __rhs);
void init(vector<unordered_map<int, double>> &__g, const string &__file_path);
void show_graph(const vector<unordered_map<int, double>> &__g);
void bellman_ford_sp(const vector<unordered_map<int, double>> &__g,
                     double *__dist_to);
void show_result(double *__dist_to);
void print_file(double *__dist_to);
void timer_start();
void timer_end();
void timer_print(const string &name);

long __nv; // number of vertex.
long __ne;
struct timeval __ts, __te;
const int __inf = 0xffff;
const string __out_name("bellman_omp.txt");
string __input_file;

int main(int __argc, char **__argv) {
  if (__argc < 2) {
    __abort("Use ./bellman [input file]");
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
  double *__dist_to;
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
  __dist_to = new double[__nv];
  timer_start();
  bellman_ford_sp(__graph, __dist_to);
  timer_end();
  //
  // Show the result
  //
  show_result(__dist_to);
  //
  // Show the time spend
  //
  timer_print("bellman-omp");
  //
  // Print to the file
  //
  print_file(__dist_to);
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
  int __offset;
  int __from, __to;
  double __distance;

  if (!__fis.is_open()) {
    __abort("Can not open the file.");
  }

  __info("start to read the file...");
  __fis >> __nv >> __ne;

  __dbg("start to create the vector...");
  __g = vector<unordered_map<int, double>>(__nv);
  __dbg("ok to create vector");
  for (__offset = 0; __offset < __nv; ++__offset) {
    __g[__offset][__offset] = 0;
  }

  for (__offset = 0; __offset < __ne; ++__offset) {
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
void bellman_ford_sp(const vector<unordered_map<int, double>> &__g,
                     double *__dist_to) {
  //
  // is the vetex on the queue.
  //
  vector<bool> __on_q(__nv);
  //
  // The queue to Support the entire algorithm.
  //
  queue<int> __q;
  //
  // Init the __dist_to.
  // The distance to ohter vetex is endless
  // And 0 to itself is 0.
  //
  int __poll;
  for (int __v = 0; __v < __nv; ++__v)
    __dist_to[__v] = __inf;
  __dist_to[0] = 0;

  __info("bellman-ford start now!");
  //
  // put the 0 to the q.
  //
  __q.push(0);
  __on_q[0] = true;

//
// Open 12 threads.
//
#pragma omp parallel num_threads(12)
  {
    int my_id = omp_get_thread_num();
    int nth = omp_get_num_threads();
    int my_first = (my_id * __nv) / nth;
    int my_last = ((my_id + 1) * __nv) / nth - 1;
    int __v; // thread local
#pragma omp critical
    {
      __info(to_string(my_id) + " start to run!, first is " +
             to_string(my_first) + ", last is " + to_string(my_last));
    }
    for (;;) {
      //
      // One thread can take the head of the queue
      //
#pragma omp single
      {
        if (__q.empty())
          __poll = -1;
        else {
          __poll = __q.front();
          __q.pop();
          __on_q[__poll] = false;
        }
        __dbg("relax [" + to_string(__poll) + "]. ");
      }
      //
      // and then relax the __poll.
      // from is __poll and then to is __v;
      //
      if (__poll != -1) {
        for (__v = my_first; __v <= my_last; ++__v) {
          if (__dist_to[__v] >
              __dist_to[__poll] + __distance(__g, __poll, __v)) {
            __dist_to[__v] = __dist_to[__poll] + __distance(__g, __poll, __v);
            //
            // if the __v is not int the queue
            // then put it to the queue
            //
#pragma omp critical
            {
              if (!__on_q[__v]) {
                __q.push(__v);
                __on_q[__v] = true;
              }
            }
          }
        }
      } 
      if (__poll == -1)
        goto end;
        //
        // This barrier...
        //
#pragma omp barrier
    }

end: {}
  }

  __info("bellman-ford end now.");
}

void show_result(double *__dist_to) {
  int __i;
  cout << "\n";
  cout << "  Minimum distances from node 0:\n";
  cout << "\n";
  for (__i = 0; __i < __nv; __i++) {
    cout << "  " << setw(2) << __i << "  " << setw(2) << __dist_to[__i] << "\n";
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
    __fis << __input_file << "\t" << __spend << " ms. " << "\t" << "bellman-openmp\t" << endl;
  }
}

void print_file(double *__dist_to) {
  if (!PRINT_FILE)
    return;
  // output to the file
  std::ofstream fos(__out_name);
  for (int i = 0; i < __nv; ++i) {
    fos << i << ": " << __dist_to[i] << endl;
  }
}
