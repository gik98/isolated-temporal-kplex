#ifndef ISOLATION_TPLEXES_HPP_
#define ISOLATION_TPLEXES_HPP_

#include <Graph.hpp>

NodeSetIntervalSet c_isolated_temporal_kplex(TGraph &g, int k, int c, TemporalIsolationType isolation);

NodeSetSet alltime_max_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start, NodeTime t_stop, int delta);

NodeSetSet max_usually_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start, NodeTime t_stop, int delta);

NodeSetSet avg_alltime_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start, NodeTime t_stop, int delta);

NodeSetSet usually_avg_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start, NodeTime t_stop, int delta);

NodeSetSet alltime_avg_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start, NodeTime t_stop, int delta);

NodeSetSet usually_max_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start, NodeTime t_stop, int delta);

#endif