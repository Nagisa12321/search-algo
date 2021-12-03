#include "dijkstra_short_path.h"
#include "directed_edge.h"
#include "edge_weighted_digraph.h"
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>

int main(int __argc, char *__argv[]) {
    if (__argc != 2) {
        std::cout << "./dijkstra [source]" << std::endl;
        return 0;
    }
    int __source = atoi(__argv[1]);
    std::ifstream __fis("../resource/tinyEWD.txt");
    sp::edge_weighted_digraph __ewd;
    __fis >> __ewd;

    sp::dijkstra_short_path __dijkstra_sp(__ewd, __source);
    for (int __i = 0; __i < __ewd.vertex(); ++__i) {
        std::cout << __source << " to " << __i;
        printf(" (%4.2f): ", __dijkstra_sp.dist_to(__i));
        if (__dijkstra_sp.has_path_to(__i)) {
            for (const sp::directed_edge &__e : __dijkstra_sp.path_to(__i)) {
                std::cout << __e << ", ";
            }
            std::cout << std::endl;
        }
    }
}