#include <isolation_scliques.hpp>

#include <cassert>

using std::shared_ptr;

bool is_vertex_cover(vector<NodeId> &elems, SGraph &g);

void min_vertex_cover_rec(int offset, int left, vector<NodeId> &elems, NodeSetSet &result, SGraph &g);

bool is_vertex_cover(vector<NodeId> &elems, SGraph &g) {
	NodeSet covered;
	covered.insert(elems.begin(), elems.end());
	for (NodeId u : elems) {
		g.forallNeighbours(u, [&](NodeId v) { covered.insert(v); }, false);
	}

#ifdef DEBUG_SCLIQUE
	assert((int)(covered.size()) <= g.getNodesCount());
#endif

	return (int)(covered.size()) == g.getNodesCount();
}

void min_vertex_cover_rec(int offset, int left, vector<NodeId> &stack, vector<NodeId> &elems, NodeSetSet &result,
			  SGraph &g) {
	if (is_vertex_cover(stack, g)) {
		/* Add to solution set and prune */
		result.insert(NodeSet(stack.begin(), stack.end()));
	} else if (left > 0) {
		for (int i = offset; i <= (int)(elems.size()) - left; i++) {
			stack.push_back(elems[i]);
			min_vertex_cover_rec(i + 1, left - 1, stack, elems, result, g);
			stack.pop_back();
		}
	}
}

NodeSetSet min_vertex_cover_bounded(SGraph &g, int c) {
	NodeSet nodes = g.getNodes();
	vector<NodeId> nodes_v;
	vector<NodeId> stack;
	NodeSetSet result;

	/* Add only nodes that have degree > 0 */
	for (const auto &u : nodes) {
		if (g.degree(u, nodes) > 0) {
			nodes_v.push_back(u);
		}
	}

	if (nodes_v.size() > 0) {
#ifdef DEBUG_SCLIQUE
		assert(nodes.size() > 0);
		assert(nodes_v.size() <= nodes.size());
		assert(c >= 0);
#endif

		min_vertex_cover_rec(0, c, stack, nodes_v, result, g);
	} else {
		result.insert(NodeSet());
	}

	return result;
}
