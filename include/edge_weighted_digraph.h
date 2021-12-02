#ifndef __EDGE_WEIGHTED_DIGRAPH_H
#define __EDGE_WEIGHTED_DIGRAPH_H
#include <istream>
#include <ostream>
#include <vector>
#include <iostream>
#include "directed_edge.h"

class edge_weighted_digraph {
    friend std::istream &operator>>(std::istream & __in, edge_weighted_digraph &__ewd) {
        int __vertex, __edges;
        __in >> __vertex >> __edges;
        for (int __v = 0; __v < __vertex; ++__v)
            __ewd._M_adj.push_back(std::vector<directed_edge>()); 
        for (int __e = 0; __e < __edges; ++__e) {
            int __from, __to;
            double __weight;
            __in >> __from >> __to >> __weight;
            __ewd._M_adj[__from].emplace_back(__from, __to, __weight);
        }
        __ewd._M_edges = __edges;
        __ewd._M_vertex = __vertex;

        return __in;
    }

    friend std::ostream &operator<<(std::ostream & __out, const edge_weighted_digraph &__ewd) {
        for (int __v = 0; __v < __ewd._M_vertex; ++__v) {
            __out << __v << ": ";
            for (int __e = 0; __e < __ewd._M_adj[__v].size(); ++__e) {
                __out << __ewd._M_adj[__v][__e] << ", ";
            }
            __out << std::endl;
        }
        return __out;
    }
public:

    edge_weighted_digraph() 
        : _M_vertex(0),
          _M_edges(0),
          _M_adj()
    {

    }
    edge_weighted_digraph(int __vertex)
        : _M_vertex(__vertex),
          _M_edges(0),
          _M_adj(__vertex)
    {

    }

    int vertex() const { return _M_vertex; }
    int edges() const { return _M_edges; }
    void add_edge(const directed_edge &__edge) {
        _M_adj[__edge.from()].push_back(__edge);
        ++_M_edges;
    }

    std::vector<directed_edge> adj(int __vertex) {
        return _M_adj[__vertex];
    }

    std::vector<directed_edge> edges() {
        std::vector<directed_edge> __res;
        for (int __v = 0; __v < _M_vertex; ++__v) {
            for (const directed_edge &__edge : _M_adj[__v]) 
                __res.push_back(__edge);
        }
        return __res;
    }


private:
    int     _M_vertex;
    int     _M_edges;
    std::vector<std::vector<directed_edge>> _M_adj;
};

#endif // __EDGE_WEIGHTED_DIGRAPH_H
