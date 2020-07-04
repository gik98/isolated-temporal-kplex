#ifndef GRAPH_HPP_
#define GRAPH_HPP_

#include <boost/functional/hash.hpp>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

using std::function;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::unordered_set;
using std::vector;

typedef int NodeId;
typedef int NodeTime;
typedef set<NodeId> NodeSet;
typedef pair<NodeSet, pair<NodeTime, NodeTime>> NodeSetInterval;

typedef pair<NodeTime, NodeTime> Interval;

typedef unordered_set<NodeSet, boost::hash<NodeSet>> NodeSetSet;
typedef unordered_set<NodeSetInterval, boost::hash<NodeSetInterval>> NodeSetIntervalSet;

#define NODETIME_MAX INT_MAX
#define NODETIME_MIN INT_MIN

enum TemporalIsolationType { ALLTIME_MAX, MAX_USUALLY, AVG_ALLTIME, USUALLY_AVG, ALLTIME_AVG, USUALLY_MAX };

struct SEdge {
	NodeId nodeFrom, nodeTo;

	SEdge() {
	}

	SEdge(NodeId nodeFrom, NodeId nodeTo) : nodeFrom(nodeFrom), nodeTo(nodeTo) {
	}

	bool operator==(const SEdge &rhs) const {
		return nodeFrom == rhs.nodeFrom && nodeTo == rhs.nodeTo;
	}
};

struct TEdge {
	NodeId nodeFrom, nodeTo;
	NodeTime tStart, tStop;

	TEdge() {
	}

	TEdge(NodeId nodeFrom, NodeId nodeTo, NodeTime tStart, NodeTime tStop)
	    : nodeFrom(nodeFrom), nodeTo(nodeTo), tStart(tStart), tStop(tStop) {
	}

	bool operator==(const TEdge &rhs) const {
		return nodeFrom == rhs.nodeFrom && nodeTo == rhs.nodeTo && tStart == rhs.tStart && tStop == rhs.tStop;
	}
};

std::size_t hash_value(SEdge const& e);

std::size_t hash_value(TEdge const& e);

class SGraph {
	map<NodeId, unordered_set<SEdge, boost::hash<SEdge>>> adj_list;

      public:
	SGraph();

	SGraph(vector<SEdge> &edgeList);

	SGraph(unordered_set<SEdge, boost::hash<SEdge>> &edgeList);

	SGraph(const SGraph &g2);

	void addNode(NodeId n);

	void addEdge(SEdge e);

	int getNodesCount();

	int getEdgesCount();

	int degree(NodeId node);

	int degree(NodeId node, const NodeSet &restriction);

	int outdegree(NodeId node, const NodeSet &in);

	int mindegree(const NodeSet &restriction);

	NodeSet neighbourhood(NodeId node);

	NodeSet neighbourhood(NodeId node, const NodeSet &restriction);

	NodeSet neighbourhoodGreaterDeg(NodeId node);

	NodeSet neighbourhoodSmallerDeg(NodeId node);

	SGraph buildComplement(const NodeSet &restriction);

	void forallNodes(function<void(NodeId)> callback, bool parallel);

	void forallNeighbours(NodeId node, function<void(NodeId)> callback, bool parallel);

	NodeSet getNodes();

	bool hasEdge(NodeId u, NodeId v);

	bool isConnected(const NodeSet &subset);

	NodeSet getReachableNodes(NodeId u);

	bool isKplex(const NodeSet &plex, int k);

	int outdegree_sum(const NodeSet &restriction);
};

class TGraph {
	map<NodeId, unordered_set<TEdge, boost::hash<TEdge>>> adj_list;
	NodeTime lifetime_begin, lifetime_end;

      public:
	TGraph();

	TGraph(vector<TEdge> &edgeList);

	TGraph(const TGraph &g2);

	void addNode(NodeId n);

	void addEdge(TEdge e);

	void forallNeighbours(NodeId node, NodeTime t, function<void(NodeId &)> callback, bool parallel);

	int getNodesCount();

	int getEdgesCount();

	int getEdgesInstantsCount();

	NodeTime getLifetimeBegin();

	NodeTime getLifetimeEnd();

	NodeSet getNodes();

	SGraph buildIntersectionGraph(const NodeSet &restriction, NodeTime t_start, NodeTime t_stop);

	SGraph buildAuxGraph(const NodeSet &restriction, NodeTime t_start, NodeTime t_stop, NodeTime t_crit);

	int outdegree(NodeId node, const NodeSet &in, NodeTime t);

	int outdegree_max(NodeId node, const NodeSet &in, NodeTime t_start, NodeTime t_stop);

	int outdegree_time_sum(NodeId node, const NodeSet &restriction, NodeTime t_start, NodeTime t_stop);

	int degree_sum(NodeId node, NodeTime t_start, NodeTime t_stop);

	void forallNodes(function<void(NodeId)> callback, bool parallel);

	bool isKplex(const NodeSet &plex, int k, NodeTime t_start, NodeTime t_stop);

	int degree(NodeId node, NodeTime t, const NodeSet &restriction);

	bool isIsolated(const NodeSet &plex, NodeTime t_start, NodeTime t_stop, TemporalIsolationType isolation, int c);

	map<NodeId, unordered_set<TEdge, boost::hash<TEdge>>> getAdjacencyList();
};

string nodeset_to_string(const NodeSet &nodeset);

string nodesetset_to_string(const NodeSetSet &nodesetset);

string nodesetinterval_to_string(const NodeSetInterval &nodeset);

#endif