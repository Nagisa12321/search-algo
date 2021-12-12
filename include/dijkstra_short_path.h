#ifndef __DIJKSTRA_SHORT_PATH_H
#define __DIJKSTRA_SHORT_PATH_H

#include "edge_weighted_digraph.h"
#include "directed_edge.h"
#include <cstdint>
#include <vector>
#include <algorithm>
#include <queue>
namespace sp {

class dijkstra_short_path {
public:
    dijkstra_short_path(const edge_weighted_digraph &__ewd, const int &__s) 
        : _M_ewd(__ewd),
          _M_s(__s),
          _M_edge_to(__ewd.vertex()),
          _M_dist_to(__ewd.vertex(), INT32_MAX)
    {
        _M_dist_to[__s] = 0;
        
        auto cmp = [&](const int &__v1, const int &__v2) {
            return _M_dist_to[__v1] > _M_dist_to[__v2];
        };
        std::priority_queue<int, std::vector<int>, decltype(cmp)> __pq(cmp);
        std::vector<bool> __relaxed(_M_ewd.vertex(), false);

        __pq.push(__s);
        while (!__pq.empty()) {
            int __from = __pq.top();
            __pq.pop();

            if (__relaxed[__from]) continue;
            __relaxed[__from] = true;

            for (directed_edge *__e : _M_ewd.adj(__from)) {
                int __to = __e->to();
                if (_M_dist_to[__to] > _M_dist_to[__from] + __e->weight()) {
                    _M_dist_to[__to] = _M_dist_to[__from] + __e->weight();
                    _M_edge_to[__to] = __e;
                    __pq.push(__to);
                }
            }

        }
    }

    double dist_to(int __v) { return _M_dist_to[__v]; }
    bool has_path_to(int __v) { return _M_dist_to[__v] < INT32_MAX; }
    std::vector<directed_edge> path_to(int __v) {
        std::vector<directed_edge> __res;
        if (!has_path_to(__v)) return __res;
        for (directed_edge *__e = _M_edge_to[__v]; __e != nullptr; __e = _M_edge_to[__e->from()])
            __res.push_back(*__e);
        std::reverse(__res.begin(), __res.end());
        return __res;
    }

private:
    const edge_weighted_digraph &_M_ewd;
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
};

}

#endif // __DIJKSTRA_SHORT_PATH_H
