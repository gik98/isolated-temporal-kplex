#include <isolation_splexes.hpp>

#include <algorithm>
#include <spdlog/spdlog.h>
#include <unordered_map>

#include <conf.hpp>

using std::set_difference;
using std::set_intersection;
using std::set_union;
using std::unordered_map;

NodeSetSet max_c_isolated_kplex(SGraph &g, int c, int k) {
	NodeSet unrestricted = g.getNodes();
	return max_c_isolated_kplex_restricted(g, c, k, unrestricted);
}

NodeSetSet max_c_isolated_kplex_restricted(SGraph &g, int c, int k, const NodeSet &restriction) {
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
						    neigh_u_size <= pivot_node_deg - c - k ||
						    g.outdegree(u, candidate) > c + k) {
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

#if 0
					/*
					 * Remove vertices which cannot be in a max-c-isolated plex
					 */
					{
						int max_del_cpy = max_del;
						const NodeSet candidate_plex_iter = candidate_plex;
						for (NodeId u : candidate_plex_iter) {
							if (g.outdegree(u, candidate_plex) >= c) {
								candidate_plex.erase(u);
								max_del_cpy--;
							}
							if (max_del_cpy < 0) {
								spdlog::trace("Dropping candidate k-plex - cannot be "
									      "max-{}-isolated. "
									      "Set: {}, max deletion: {}",
									      c, nodeset_to_string(candidate_plex_iter),
									      max_del);
								return;
							}
						}
					}
#endif

					/*
					 * Compute meaningful k-plexes: a k-plex shall have at least k + 2 vertices
					 * Moreover, we are interested in connected k-plexes only
					 */
					NodeSetSet screening_candidates;

					int bdd_max_del = std::min(max_del, (int)(candidate_plex.size()) - k - 2);

					if (bdd_max_del < 0) {
						goto next_kplex;
					} else {
						NodeSetSet bdd_sets;

						if (bdd_max_del == 0) {
							NodeSet empty;
							if (g.isKplex(candidate_plex, k)) {
								bdd_sets.insert(empty);
							}
						} else {
							search_graph = g.buildComplement(candidate_plex);

							bdd_sets =
							    min_bdd_d_set(search_graph, bdd_max_del, k - 1, candidate);
						}
#ifdef DEBUG_23MAY
						{
							if (candidate_plex.size() > DEBUG_23MAY) {
								spdlog::info("Found {} bdd-{} set with max size {}, "
									     "candidate set size {}",
									     bdd_sets.size(), k - 1, max_del,
									     candidate.size());
							}
							spdlog::trace("Found {} bdd-{} set with max size {}, candidate "
								      "set size {}",
								      bdd_sets.size(), k - 1, max_del,
								      candidate.size());
						}
#endif

						for (const NodeSet &bdd_set : bdd_sets) {
							NodeSet plex = candidate_plex;

							for (NodeId deletion : bdd_set) {
								plex.erase(deletion);
							}

							if (!g.isKplex(plex, k)) {
								spdlog::error("bdd enumerated a set which is not a "
									      "plex! {}; candidate is {}, bdd is {}",
									      nodeset_to_string(plex),
									      nodeset_to_string(candidate_plex),
									      nodeset_to_string(bdd_set));
							}

							int max_del_screening = bdd_max_del - (int)(bdd_set.size());
							if (max_del_screening < 0) {
								spdlog::error("max del screening < 0");
							}

							NodeSet plex_prime = plex;
							bool fixpoint = false;
							while (!fixpoint) {
								fixpoint = true;
								for (NodeId u : plex_prime) {
									if (g.outdegree(u, plex) >= c) {
										fixpoint = false;
										if (candidate.find(u) !=
										    candidate.end()) {
											/* This vertex is in the
											 * candidate set, try to remove
											 * it */
											plex.erase(u);
											max_del_screening--;

											if (max_del_screening < 0) {
												/* Too many vertices
												 * removed */
												goto break_outer_loop;
											}
										} else {
											/* This vertex is in the pivot
											 * set, drop the plex */
											goto break_outer_loop;
										}
									}
								}
								plex_prime = plex;
							}
						break_outer_loop:

							if (fixpoint) {
								screening_candidates.insert(plex);
							}
						}
					}

					for (const NodeSet &plex : screening_candidates) {

						if (!g.isKplex(plex, k)) {
							spdlog::error("After isolation screening, {} is not a {}-plex", nodeset_to_string(plex), k);
						}
						/* Screening #0: is max-c-isolated? */
						bool isolated = true;
						for (NodeId u : plex) {
							if (g.outdegree(u, plex) >= c) {
								isolated = false;
								break;
							}
						}

						/* Screening #1: pivot vertex check */
						if (isolated) {
							bool maximal = true;
							for (NodeId u : plex) {
								if (g.degree(u) < pivot_node_deg &&
								    g.outdegree(u, plex) < c) {
									/* Not maximal - drop */
									maximal = false;
									break;
								}
							}

							if (maximal) {
#pragma omp critical(sol)
								sol.insert(plex);
							} else {
								spdlog::trace("Ignoring k-plex - failed pivot vertex "
									      "check (rule #1).");
							}
						} else {
							spdlog::error(
							    "Ignoring k-plex - is not max-{}-isolated. Set: {}", c,
							    nodeset_to_string(plex));
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