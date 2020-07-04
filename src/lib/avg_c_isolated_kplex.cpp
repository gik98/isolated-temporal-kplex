#include <isolation_splexes.hpp>

#include <algorithm>
#include <spdlog/spdlog.h>
#include <unordered_map>

#include <conf.hpp>

using std::max;
using std::min;
using std::set_difference;
using std::set_intersection;
using std::set_union;
using std::sort;
using std::unordered_map;

NodeSetSet avg_isolated_subsets(SGraph &g, int k, int c, NodeSet &candidate, int max_del) {
	NodeSetSet avg_subsets;
	vector<NodeId> deg_sorted(candidate.begin(), candidate.end());
	sort(deg_sorted.begin(), deg_sorted.end(), [&](NodeId a, NodeId b) { return g.degree(a) < g.degree(b); });

	NodeSetSet deletions, deletions_prime;
	deletions_prime.insert(NodeSet());

	int mindegree = g.mindegree(candidate);
	while ((deletions = deletions_prime).size() > 0) {
		deletions_prime.clear();
		for (const NodeSet &del : deletions) {
			NodeSet candidate_del;
			set_difference(candidate.begin(), candidate.end(), del.begin(), del.end(),
				       std::inserter(candidate_del, candidate_del.begin()));

			if (g.outdegree_sum(candidate_del) >= (int)(candidate_del.size()) * c) {
				int max_deletion = min((int)(max_del) - (int)(del.size()),
						       (mindegree - c - 3 * k) - (int)(del.size()));
				if (max_deletion > 0) {
					/* Create new deletions sets, adding candidates */
					for (int idx = max(mindegree - c - 3 * k, 0); idx < (int)(deg_sorted.size());
					     idx++) {
						if (del.find(deg_sorted[idx]) == del.end()) {
							NodeSet delcpy = del;
							delcpy.insert(deg_sorted[idx]);

							deletions_prime.insert(delcpy);
						}
					}
				}
				/*
				 * else this cannot be a avg-c-isolated subset
				 */
			} else {
				/* This is a valid avg-c-isolated k-plex. Store it */
				avg_subsets.insert(candidate_del);
			}
		}
	}

	return avg_subsets;
}

NodeSetSet avg_c_isolated_kplex(SGraph &g, int c, int k) {
	NodeSet unrestricted = g.getNodes();
	return avg_c_isolated_kplex_restricted(g, c, k, unrestricted);
}

NodeSetSet avg_c_isolated_kplex_restricted(SGraph &g, int c, int k, const NodeSet &restriction) {
	NodeSetSet sol, ret;

#pragma omp parallel if (parallelism)
	{
#pragma omp single
		for (NodeId pivot_node : restriction) {
			spdlog::trace("Pivot node {}", pivot_node);
#pragma omp task if (parallelism)
			{
				SGraph compl_graph, search_graph;
				NodeSet pivot_candidate, pivot_neigh, node_set, node_set_restricted;

				/* Candidate set */
				NodeSet candidate = pivot_neigh = g.neighbourhood(pivot_node);
				pivot_neigh.insert(pivot_node);

				int pivot_node_deg = candidate.size();

				/* Trimming stage */
				int max_del = c - 1;
				bool fixpoint = false;
				while (!fixpoint) {
					fixpoint = true;
					const NodeSet candidate_iter = candidate;

					for (NodeId u : candidate_iter) {
						int neigh_u_size = g.degree(u, candidate);

						if (restriction.find(u) == restriction.end() ||
						    g.outdegree(u, candidate) >= (k + (int)candidate.size()) * c ||
						    neigh_u_size <= pivot_node_deg - c - k) {
							fixpoint = false;
							candidate.erase(u);
							max_del--;
						}

						if (max_del < 0) {
							goto next_pivot;
						}
					}
				}

				/* Enumeration stage */

				/* We are interested in reachable nodes only */
				node_set = g.getReachableNodes(pivot_node);
				set_intersection(node_set.begin(), node_set.end(), restriction.begin(),
						 restriction.end(),
						 std::inserter(node_set_restricted, node_set_restricted.begin()));

				set_difference(node_set_restricted.begin(), node_set_restricted.end(),
					       pivot_neigh.begin(), pivot_neigh.end(),
					       std::inserter(pivot_candidate, pivot_candidate.begin()));

				foreach_kplex_pivot(k - 1, pivot_candidate, [&](NodeSet &pivot_set) {
					NodeSet candidate_plex;

					pivot_set.insert(pivot_node);

					set_union(pivot_set.begin(), pivot_set.end(), candidate.begin(),
						  candidate.end(),
						  std::inserter(candidate_plex, candidate_plex.begin()));

					/*
					 * Compute meaningful k-plexes: a k-plex shall have at least k + 2 vertices
					 * Moreover, we are interested in connected k-plexes only
					 */
					NodeSetSet screening_candidates;

					int bdd_max_del = std::min(max_del, (int)(candidate_plex.size()) - k - 2);

					if (bdd_max_del < 0) {
						goto next_kplex;
					} else if (bdd_max_del == 0) {
						NodeSet plex = candidate_plex;

						if (g.isKplex(plex, k) &&
						    g.outdegree_sum(plex) < c * (int)(plex.size())) {
							screening_candidates.insert(plex);
						}

					} else {
						search_graph = g.buildComplement(candidate_plex);

						NodeSetSet bdd_sets =
						    min_bdd_d_set(search_graph, bdd_max_del, k - 1, candidate);

						NodeSet plex;
						for (const NodeSet &bdd_set : bdd_sets) {
							plex = candidate_plex;

							for (NodeId deletion : bdd_set) {
								plex.erase(deletion);
							}

							if (!g.isKplex(plex, k)) {
								spdlog::error("{} is not a {}-plex! bdd-set: {}",
									      nodeset_to_string(plex), k,
									      nodeset_to_string(bdd_set));
							}

							/*
							 * Forward to screening avg-isolated subsets only
							 */

							NodeSetSet isolated_subsets = avg_isolated_subsets(
							    g, k, c, plex, bdd_max_del - (int)(bdd_set.size()));

							for (const NodeSet &plex_avg : isolated_subsets) {
								screening_candidates.insert(plex_avg);
							}
						}
					}

					for (const NodeSet &plex : screening_candidates) {
						if (!g.isKplex(plex, k)) {
							spdlog::error("After isolation screening, {} is not a {}-plex!",
								      nodeset_to_string(plex), k);
						}
						/* Screening #1: pivot vertex check */
						bool maximal = true;
						for (NodeId u : plex) {
							if (g.degree(u) < pivot_node_deg && g.outdegree(u, plex) < c) {
								/* Not maximal - drop */
								maximal = false;
								break;
							}
						}

#pragma omp critical(screening_candidates)
						if (maximal) {
							sol.insert(plex);
						} else {
							spdlog::trace(
							    "Ignoring k-plex - failed pivot vertex check (rule #1).");
						}
					}

				next_kplex:;
				});

			next_pivot:;
			}
		}
	}
	spdlog::trace("Enumeration stage returned {} {}-plexes", sol.size(), k);

/* Screening #2: maximality */
#pragma omp parallel if (parallelism)
	{
#pragma omp single

		{
			for (const NodeSet &s : sol) {
#pragma omp task if (parallelism)
				{
					bool is_maximal = true;
					for (const NodeSet &t : sol) {
						if (s == t) {
							continue;
						}
						if (t.size() < s.size()) {
							continue;
						}

						NodeSet isect;
						set_intersection(t.begin(), t.end(), s.begin(), s.end(),
								 std::inserter(isect, isect.begin()));

						if (isect == s) {
							is_maximal = false;
						}

						if (!is_maximal) {
							spdlog::trace(
							    "Ignoring k-plex - failed maximality check (rule #2).");
							break;
						}
					}

					if (is_maximal) {
#pragma omp critical(ret)
						ret.insert(s);
					}
				}
			}
		}
	}

	return ret;
}