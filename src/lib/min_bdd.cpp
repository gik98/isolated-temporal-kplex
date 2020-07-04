#include <isolation_splexes.hpp>

#include <algorithm>
#include <memory>

#include <Graph.hpp>

using std::make_shared;
using std::set_difference;

bool is_min_bdd_d(SGraph &g, NodeSet &deletion, int d) {
	bool crit = true;
	g.forallNodes([&](NodeId u) {
		int deg = g.outdegree(u, deletion);
		if (deg > d) {
			crit = false;
		}
	}, false);

	return crit;
}

void min_bdd_search(SGraph &g, NodeSet &deletion, int d, int k, NodeSet &candidate_set, NodeSetSet &result) {
	if (is_min_bdd_d(g, deletion, d)) {
		/* Add to solution set and prune */
		result.insert(deletion);
	} else if (k > 0) {
		bool flag_1 = true, flag_2 = true;
		g.forallNodes([&](NodeId u) {
			if (flag_1) {
				if (deletion.find(u) != deletion.end()) {
					/* Skip nodes which are already part of the minbdd */
					return;
				}

				if (g.outdegree(u, deletion) > d) {
					/* Recursion */
					g.forallNeighbours(u, [&](NodeId v) {
						if (flag_2) {
							if (deletion.find(v) != deletion.end()) {
								return;
							}

							if (candidate_set.find(v) == candidate_set.end()) {
								/* We cannot delete vertices which do not belong to the
								 * candidate set */
								return;
							}

							/*
							 * We will need to work on a copy if we will go parallel
							 *
							 * NodeSet deletion_cpy = deletion;
							 * deletion_cpy.insert(v);
							 */
							deletion.insert(v);
							min_bdd_search(g, deletion, d, k - 1, candidate_set, result);
							deletion.erase(v);
						}
					}, false);

					deletion.insert(u);
					min_bdd_search(g, deletion, d, k - 1, candidate_set, result);
					deletion.erase(u);

					flag_1 = false;
				}
			}
		}, false);
	}
}

NodeSetSet min_bdd_d_set(SGraph &g, int max_del, int d, NodeSet &candidate_set) {
	NodeSetSet sol, ret;
	/*
	 * First step: remove all vertices with a degree greater than k + d.
	 * They necessairly are in the solution of every min bdd
	 */
	int curr_k = max_del;
	NodeSet kernel;
	bool skip = false;

	g.forallNodes([&](NodeId u) {
		if (g.degree(u) > max_del + d) {
			if (candidate_set.find(u) != candidate_set.end()) {
				kernel.insert(u);
				curr_k--;
			} else {
				/* We cannot remove vertices which do not belong to the candidate set */
				skip = true;
			}
		}
	}, false);

	if (curr_k < 0 || skip) {
		return sol;
	}

	/* Second step: start enumerating in bdd */
	min_bdd_search(g, kernel, d, max_del, candidate_set, sol);

	/* Third step: maximality check */
	for (const NodeSet &s : sol) {
		bool maximal = true;
		for (const NodeSet &t : sol) {
			if (s == t || s.size() > t.size()) {
				continue;
			}
			if (s.size() < t.size()) {
				maximal = false;
				break;
			}

			NodeSet intr;
			set_difference(s.begin(), s.end(), t.begin(), t.end(), std::inserter(intr, intr.begin()));
			if (intr.size() == 0) {
				maximal = false;
				break;
			}
		}

		if (maximal) {
			ret.insert(s);
		}
	}

	return ret;
}

void foreach_kplex_pivot_rec(int offset, int left, vector<NodeId> &pivot_candidates, vector<NodeId> &stack,
			     function<void(NodeSet &)> callback) {
	if (stack.size() > 0) {
		NodeSet pivot(stack.begin(), stack.end());
		callback(pivot);
	}
	if (left > 0) {
		for (int i = offset; i <= (int)(pivot_candidates.size()) - left; i++) {
			stack.push_back(pivot_candidates[i]);
			foreach_kplex_pivot_rec(i + 1, left - 1, pivot_candidates, stack, callback);
			stack.pop_back();
		}
	}
}

void foreach_kplex_pivot(int k, NodeSet &pivot_candidates, function<void(NodeSet &)> callback) {
	vector<NodeId> pivot_candidate_vec(pivot_candidates.begin(), pivot_candidates.end());
	vector<NodeId> stack;

	/* Empty set is subset of every set */
	NodeSet empty;
	callback(empty);

	foreach_kplex_pivot_rec(0, k, pivot_candidate_vec, stack, callback);
}