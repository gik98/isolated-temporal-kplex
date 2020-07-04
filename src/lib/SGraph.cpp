#include <Graph.hpp>

#include <algorithm>
#include <stack>

#include <spdlog/spdlog.h>

using std::find;
using std::stack;

std::size_t hash_value(SEdge const &e) {
	size_t seed = 0;

	boost::hash_combine(seed, e.nodeFrom);
	boost::hash_combine(seed, e.nodeTo);
	return seed;
}

SGraph::SGraph() {
}

SGraph::SGraph(vector<SEdge> &edgeList) {
	for (const SEdge &e : edgeList) {
		(this->adj_list[e.nodeFrom]).insert(e);
		(this->adj_list[e.nodeTo]).insert(SEdge(e.nodeTo, e.nodeFrom));
	}
}

SGraph::SGraph(const SGraph &g2) {
	this->adj_list = g2.adj_list;
}

SGraph::SGraph(unordered_set<SEdge, boost::hash<SEdge>> &edgeList) {
	for (const SEdge &e : edgeList) {
		(this->adj_list[e.nodeFrom]).insert(e);
		if (edgeList.find(SEdge(e.nodeTo, e.nodeFrom)) == edgeList.end()) {
			(this->adj_list[e.nodeTo]).insert(SEdge(e.nodeTo, e.nodeFrom));
		}
	}
}

void SGraph::addNode(NodeId n) {
	if (this->adj_list.find(n) == this->adj_list.end()) {
		this->adj_list[n] = unordered_set<SEdge, boost::hash<SEdge>>();
	}
}

void SGraph::addEdge(SEdge e) {
	(this->adj_list[e.nodeFrom]).insert(e);
	(this->adj_list[e.nodeTo]).insert(SEdge(e.nodeTo, e.nodeFrom));
}

int SGraph::getNodesCount() {
	return (int)(this->adj_list.size());
}

int SGraph::getEdgesCount() {
	int cnt = 0;
	for (const auto &l : this->adj_list) {
		cnt += (int)(l.second.size());
	}

	return cnt / 2;
}

int SGraph::degree(NodeId node) {
	return (int)((this->adj_list[node]).size());
}

int SGraph::degree(NodeId node, const NodeSet &restriction) {
	int deg = 0;

	this->forallNeighbours(node,
			       [&](NodeId n) {
				       if (restriction.find(n) != restriction.end()) {
					       deg++;
				       }
			       },
			       false);

	return deg;
}

int SGraph::mindegree(const NodeSet &restriction) {
	int mindeg = 0;
	for (NodeId u : restriction) {
		int deg = this->degree(u);
		if (deg < mindeg) {
			mindeg = deg;
		}
	}

	return mindeg;
}

int SGraph::outdegree(NodeId node, const NodeSet &in) {
	int deg = 0;

	this->forallNeighbours(node,
			       [&](NodeId n) {
				       if (in.find(n) == in.end()) {
					       deg++;
				       }
			       },
			       false);

	return deg;
}

NodeSet SGraph::neighbourhood(NodeId node) {
	NodeSet v;
	for (SEdge e : this->adj_list[node]) {
		v.insert(e.nodeTo);
	}
	return v;
}

NodeSet SGraph::neighbourhood(NodeId node, const NodeSet &restriction) {
	NodeSet v;
	for (SEdge e : this->adj_list[node]) {
		if (restriction.find(e.nodeTo) != restriction.end()) {
			v.insert(e.nodeTo);
		}
	}
	return v;
}

NodeSet SGraph::neighbourhoodGreaterDeg(NodeId node) {
	NodeSet v;

	this->forallNeighbours(node,
			       [&](NodeId n) {
				       if (this->degree(n) >= this->degree(node)) {
					       v.insert(n);
				       }
			       },
			       false);

	return v;
}

NodeSet SGraph::neighbourhoodSmallerDeg(NodeId node) {
	NodeSet v;

	this->forallNeighbours(node,
			       [&](NodeId n) {
				       if (this->degree(n) <= this->degree(node)) {
					       v.insert(n);
				       }
			       },
			       false);

	return v;
}

SGraph SGraph::buildComplement(const NodeSet &restriction) {
	SGraph ret;

	for (NodeId u : restriction) {
		ret.addNode(u);
		for (NodeId v : restriction) {
			if (u != v && !this->hasEdge(u, v)) {
				ret.addEdge(SEdge(u, v));
			}
		}
	}

	return ret;
}

void SGraph::forallNodes(function<void(NodeId)> callback, bool parallel) {
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

void SGraph::forallNeighbours(NodeId node, function<void(NodeId)> callback, bool parallel) {
#pragma omp parallel if (parallel)
	{
#pragma omp single
		{
			for (SEdge n : this->adj_list[node]) {
#pragma omp task if (parallel)
				callback(n.nodeTo);
			}
		}
	}
}

NodeSet SGraph::getNodes() {
	NodeSet ret;
	for (const auto &entry : this->adj_list) {
		ret.insert(entry.first);
	}

	return ret;
}

bool SGraph::hasEdge(NodeId u, NodeId v) {
	return find(this->adj_list[u].begin(), this->adj_list[u].end(), SEdge(u, v)) != this->adj_list[u].end();
}

bool SGraph::isConnected(const NodeSet &subset) {
	if (subset.empty()) {
		return true;
	}

	set<NodeId> visited;
	stack<NodeId> dfs;

	dfs.push(*subset.begin());
	visited.insert(*subset.begin());

	while (!dfs.empty()) {
		NodeId u = dfs.top();
		dfs.pop();

		this->forallNeighbours(u,
				       [&](NodeId v) {
					       if (subset.find(v) != subset.end() && visited.find(v) == subset.end()) {
						       visited.insert(v);
						       dfs.push(v);
					       }
				       },
				       false);
	}

	return visited.size() == subset.size();
}

NodeSet SGraph::getReachableNodes(NodeId u) {
	NodeSet res;

	stack<NodeId> dfs;
	dfs.push(u);
	res.insert(u);

	while (!dfs.empty()) {
		NodeId u = dfs.top();
		dfs.pop();

		this->forallNeighbours(u,
				       [&](NodeId v) {
					       if (res.find(v) == res.end()) {
						       res.insert(v);
						       dfs.push(v);
					       }
				       },
				       false);
	}

	return res;
}

bool SGraph::isKplex(const NodeSet &plex, int k) {
	for (NodeId u : plex) {
		if (this->degree(u, plex) < (int)(plex.size()) - k) {
			return false;
		}
	}

	return true;
}

int SGraph::outdegree_sum(const NodeSet &restriction) {
	int sum = 0;
	for (NodeId u : restriction) {
		sum += this->outdegree(u, restriction);
	}

	return sum;
}

string nodeset_to_string(const NodeSet &nodeset) {
	std::stringstream ss;
	ss << "Set size: " << nodeset.size() << "; elems: ";
	for (NodeId i : nodeset) {
		ss << i << " ";
	}

	return ss.str();
}

string nodesetset_to_string(const NodeSetSet &nodesetset) {
	std::stringstream ss;
	for (const NodeSet &s: nodesetset) {
		ss << nodeset_to_string(s) << std::endl;
	}

	return ss.str();
}