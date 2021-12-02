#include "edge_weighted_digraph.h"
#include <fstream>
#include <iostream>

int main() {
    std::ifstream __fis("../resource/tinyEWD.txt");
    edge_weighted_digraph __ewd;
    __fis >> __ewd;
    std::cout << __ewd;
}