#ifndef __DIRECTED_EDGE_H
#define __DIRECTED_EDGE_H
#include <iostream>

namespace sp {

class directed_edge {
    friend std::ostream &operator<<(std::ostream & __out, const directed_edge &__de) {
        __out << "[" << __de._M_v 
              << " -> " << __de._M_w 
              << ", " << __de._M_weight << "]";
        return __out;
    }
public:
    directed_edge(int __v, int __w, double __weight) 
        : _M_v(__v),
          _M_w(__w),
          _M_weight(__weight)
    {

    }

    double weight() const { return _M_weight; }
    int from() const { return _M_v; }
    int to() const { return _M_w; }
private:
    int         _M_v;
    int         _M_w;
    double      _M_weight;
};

}

#endif // __DIRECTED_EDGE_H
