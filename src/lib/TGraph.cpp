#include <Graph.hpp>

#include <algorithm>
#include <stack>

#include <spdlog/spdlog.h>

using std::max;
using std::min;

std::size_t hash_value(TEdge const &e) {
	size_t seed = 0;

	boost::hash_combine(seed, e.nodeFrom);
	boost::hash_combine(seed, e.nodeTo);
	boost::hash_combine(seed, e.tStart);
	boost::hash_combine(seed, e.tStop);
	return seed;
}

TGraph::TGraph() {
	this->lifetime_begin = NODETIME_MAX;
	this->lifetime_end = NODETIME_MIN;
}

TGraph::TGraph(vector<TEdge> &edgeList) {
	this->lifetime_begin = NODETIME_MAX;
	this->lifetime_end = NODETIME_MIN;

	for (TEdge &e : edgeList) {
		(this->adj_list[e.nodeFrom]).insert(e);
		(this->adj_list[e.nodeTo]).insert(TEdge(e.nodeTo, e.nodeFrom, e.tStart, e.tStop));

		this->lifetime_begin = min(this->lifetime_begin, e.tStart);
		this->lifetime_end = max(this->lifetime_end, e.tStop);
	}
}

TGraph::TGraph(const TGraph &g2) {
	this->adj_list = g2.adj_list;
	this->lifetime_begin = g2.lifetime_begin;
	this->lifetime_end = g2.lifetime_end;
}

void TGraph::addNode(NodeId n) {
	if (this->adj_list.find(n) == this->adj_list.end()) {
		this->adj_list[n] = unordered_set<TEdge, boost::hash<TEdge>>();
	}
}

void TGraph::addEdge(TEdge e) {
	(this->adj_list[e.nodeFrom]).insert(e);
	(this->adj_list[e.nodeTo]).insert(TEdge(e.nodeTo, e.nodeFrom, e.tStart, e.tStop));

	this->lifetime_begin = min(this->lifetime_begin, e.tStart);
	this->lifetime_end = max(this->lifetime_end, e.tStop);
}

int TGraph::getNodesCount() {
	return (int)(this->adj_list.size());
}

int TGraph::getEdgesCount() {
	int cnt = 0;
	for (const auto &l : this->adj_list) {
		cnt += (int)(l.second.size());
	}

	return cnt / 2;
}

int TGraph::getEdgesInstantsCount() {
	int cnt = 0;
	for (const auto &l : this->adj_list) {
		for (const TEdge e : l.second) {
			cnt += e.tStop - e.tStart + 1;
		}
	}

	return cnt / 2;
}

NodeTime TGraph::getLifetimeBegin() {
	return this->lifetime_begin;
}

NodeTime TGraph::getLifetimeEnd() {
	return this->lifetime_end;
}

void TGraph::forallNeighbours(NodeId node, NodeTime t, function<void(NodeId &)> callback, bool parallel) {
#pragma omp parallel if (parallel)
	{
#pragma omp single
		{
			for (TEdge n : this->adj_list[node]) {
				if (n.tStart <= t && t <= n.tStop) {
#pragma omp task if (parallel)
					callback(n.nodeTo);
				}
			}
		}
	}
}

NodeSet TGraph::getNodes() {
	NodeSet ret;
	for (const auto &entry : this->adj_list) {
		ret.insert(entry.first);
	}

	return ret;
}

/* TODO: note here, we do not want "fragmented" edges */
SGraph TGraph::buildIntersectionGraph(const NodeSet &restriction, NodeTime t_start, NodeTime t_stop) {
	unordered_set<SEdge, boost::hash<SEdge>> adj;
	for (NodeId u : restriction) {
		for (TEdge e : this->adj_list[u]) {
			if (e.tStart <= t_start && t_stop <= e.tStop) {
				adj.insert(SEdge(e.nodeFrom, e.nodeTo));
			}
		}
	}

	return SGraph(adj);
}

SGraph TGraph::buildAuxGraph(const NodeSet &restriction, NodeTime t_start, NodeTime t_stop, NodeTime t_crit) {
	/* First, get all edges in the intersection graph of V \ restriction */
	NodeSet outer = this->getNodes();
	for (NodeId u : restriction) {
		outer.erase(u);
	}

	SGraph intersect = this->buildIntersectionGraph(outer, t_start, t_stop);

	/* Then, add edges internal to the restriction subset which are valid at the crit time */
	for (NodeId u : restriction) {
		this->forallNeighbours(u, t_crit,
				       [&](NodeId v) {
					       if (restriction.find(v) != restriction.end()) {
						       intersect.addEdge(SEdge(u, v));
					       }
				       },
				       false);
	}

	return intersect;
}

int TGraph::outdegree(NodeId node, const NodeSet &in, NodeTime t) {
	int outdeg = 0;
	this->forallNeighbours(node, t,
			       [&](NodeId v) {
				       if (in.find(v) == in.end()) {
					       outdeg++;
				       }
			       },
			       false);

	return outdeg;
}

int TGraph::outdegree_max(NodeId node, const NodeSet &in, NodeTime t_start, NodeTime t_stop) {
	int outdeg_max = -1;
	for (NodeTime i = t_start; i <= t_stop; i++) {
		int outdegree = 0;
		this->forallNeighbours(node, i,
				       [&](NodeId v) {
					       if (in.find(v) == in.end()) {
						       outdegree++;
					       }
				       },
				       false);
		outdeg_max = max(outdeg_max, outdegree);
	}

	return outdeg_max;
}

int TGraph::degree_sum(NodeId node, NodeTime t_start, NodeTime t_stop) {
	int res = 0;

	for (TEdge e : this->adj_list[node]) {
		NodeTime start = max(t_start, e.tStart);
		NodeTime stop = min(t_stop, e.tStop);

		if (stop - start >= 0) {
			res += stop - start + 1;
		}
	}

	return res;
}

int TGraph::outdegree_time_sum(NodeId node, const NodeSet &restriction, NodeTime t_start, NodeTime t_stop) {
	int res = 0;

	for (TEdge e : this->adj_list[node]) {
		if (restriction.find(e.nodeTo) == restriction.end()) {
			NodeTime start = max(t_start, e.tStart);
			NodeTime stop = min(t_stop, e.tStop);

			if (stop - start >= 0) {
				res += stop - start + 1;
			}
		}
	}

	return res;
}

void TGraph::forallNodes(function<void(NodeId)> callback, bool parallel) {
#pragma omp parallel if (parallel)
	{
#pragma omp single
		{
			for (const auto &entry : this->adj_list) {
#pragma omp task if (parallel)
				callback(entry.first);
			}
		}
	}
}

int TGraph::degree(NodeId node, NodeTime t, const NodeSet &restriction) {
	int deg = 0;

	this->forallNeighbours(node, t,
			       [&](NodeId v) {
				       if (restriction.find(v) != restriction.end()) {
					       deg++;
				       }
			       },
			       false);

	return deg;
}

bool TGraph::isKplex(const NodeSet &plex, int k, NodeTime t_start, NodeTime t_stop) {
	for (NodeTime t = t_start; t <= t_stop; t++) {
		for (NodeId v : plex) {
			if (this->degree(v, t, plex) < (int)plex.size() - k) {
				return false;
			}
		}
	}
	return true;
}

bool TGraph::isIsolated(const NodeSet &plex, NodeTime t_start, NodeTime t_stop, TemporalIsolationType isolation,
			int c) {
	switch (isolation) {
	case ALLTIME_MAX: {
		for (NodeId v : plex) {
			if (this->outdegree_max(v, plex, t_start, t_stop) >= c) {
				return false;
			}
		}

		return true;
	}
	case USUALLY_AVG: {
		long val = 0;
		for (NodeId v : plex) {
			val += this->outdegree_time_sum(v, plex, t_start, t_stop);
		}

		return val < c * (long)(t_stop - t_start + 1) * (long)plex.size();
	}
	case USUALLY_MAX: {
		long val = 0;
		for (NodeTime t = t_start; t <= t_stop; t++) {
			int max_outdeg = 0;
			for (NodeId v : plex) {
				int outdeg = this->outdegree(v, plex, t);
				if (outdeg > max_outdeg) {
					max_outdeg = outdeg;
				}
			}
		}

		return val < (long)(t_stop - t_start + 1) * c;
	}

	case ALLTIME_AVG: {
		for (NodeTime t = t_start; t <= t_stop; t++) {
			int condition = 0;

			for (NodeId v : plex) {
				condition += this->outdegree(v, plex, t);
			}

			if (condition >= c * (int)plex.size()) {
				return false;
			}
		}

		return true;
	}
	case MAX_USUALLY:
	case AVG_ALLTIME: { /* TODO */
		return false;
	}
	}

	return false;
}

map<NodeId, unordered_set<TEdge, boost::hash<TEdge>>> TGraph::getAdjacencyList() {
	return this->adj_list;
}

string nodesetinterval_to_string(const NodeSetInterval &nodeset) {
	std::stringstream ss;
	ss << "Set size: " << nodeset.first.size() << ", interval[" << nodeset.second.first << ","
	   << nodeset.second.second << "]; elems: ";
	for (NodeId i : nodeset.first) {
		ss << i << " ";
	}

	return ss.str();
}
