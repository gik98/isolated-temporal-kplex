#ifndef ISOLATION_SPLEXES_H
#define ISOLATION_SPLEXES_H

#include <Graph.hpp>

NodeSetSet min_bdd_d_set(SGraph &g, int k, int d, NodeSet &candidate);

void foreach_kplex_pivot(int k, NodeSet &pivot_candidates, function<void(NodeSet &)> callback);

NodeSetSet min_c_isolated_kplex(SGraph &g, int c, int k);

NodeSetSet max_c_isolated_kplex(SGraph &g, int c, int k);

NodeSetSet avg_c_isolated_kplex(SGraph &g, int c, int k);

NodeSetSet max_c_isolated_kplex_restricted(SGraph &g, int c, int k, const NodeSet &restriction);

NodeSetSet avg_c_isolated_kplex_restricted(SGraph &g, int c, int k, const NodeSet &restriction);

#endif