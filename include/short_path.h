#ifndef __SHORT_PATH_H
#define __SHORT_PATH_H

#include "edge_weighted_digraph.h"
#include "directed_edge.h"
#include <cstdint>
#include <vector>
#include <algorithm>
namespace sp {

class short_path {
public:
    short_path(const edge_weighted_digraph &__ewd, const int &__s) 
        : _M_ewd(__ewd),
          _M_s(__s),
          _M_edge_to(__ewd.vertex()),
          _M_dist_to(__ewd.vertex(), INT32_MAX)
    {
        _M_dist_to[__s] = 0;
    }

    double dist_to(int __v) { return _M_dist_to[__v]; }
    bool has_path_to(int __v) { return _M_dist_to[__v] < INT32_MAX; }
    std::vector<directed_edge> path_to(int __v) {
        std::vector<directed_edge> __res;
        if (!has_path_to(__v)) return __res;
        for (directed_edge *__e = _M_edge_to[__v]; __e; __e = _M_edge_to[__e->from()])
            __res.push_back(*__e);
        std::reverse(__res.begin(), __res.end());
        return __res;
    }

private:
    edge_weighted_digraph _M_ewd;
    int _M_s;
    std::vector<directed_edge *> _M_edge_to;
    std::vector<double> _M_dist_to;

    void relax(directed_edge &__e) {
        int __from = __e.from(), __to = __e.to();
        if (_M_dist_to[__to] > _M_dist_to[__from] + __e.weight()) {
            _M_dist_to[__to] = _M_dist_to[__from] + __e.weight();
            _M_edge_to[__to] = &__e;
        }
    }

    void relax(int __from) {
        for (directed_edge *__e : _M_ewd.adj(__from)) {
            int __to = __e->to();
            if (_M_dist_to[__to] > _M_dist_to[__from] + __e->weight()) {
                _M_dist_to[__to] = _M_dist_to[__from] + __e->weight();
                _M_edge_to[__to] = __e;
            }
        }
    }
};

}

#endif // __SHORT_PATH_H
