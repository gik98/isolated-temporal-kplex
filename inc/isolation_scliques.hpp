#ifndef ISOLATION_SCLIQUES_H_
#define ISOLATION_SCLIQUES_H_

#include <Graph.hpp>

NodeSetSet min_vertex_cover_bounded(SGraph &g, int c);

NodeSetSet max_c_isolated_clique(SGraph &g, int c);

NodeSetSet min_c_isolated_clique(SGraph &g, int c);

NodeSetSet avg_c_isolated_clique(SGraph &g, int c);

#endif