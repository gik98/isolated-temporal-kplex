#include <isolation_splexes.hpp>

#include <algorithm>
#include <spdlog/spdlog.h>
#include <unordered_map>

#include <conf.hpp>

using std::set_difference;
using std::set_intersection;
using std::set_union;
using std::unordered_map;

NodeSetSet min_c_isolated_kplex(SGraph &g, int c, int k) {
	NodeSetSet sol, ret;

	g.forallNodes(
	    [&](NodeId pivot_node) {
		    spdlog::trace("Pivot node {}", pivot_node);
		    SGraph compl_graph, search_graph;
		    NodeSet pivot_candidate, pivot_neigh, node_set;

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

				    if (neigh_u_size <= pivot_node_deg - c - k) {
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
		    set_difference(node_set.begin(), node_set.end(), pivot_neigh.begin(), pivot_neigh.end(),
				   std::inserter(pivot_candidate, pivot_candidate.begin()));

		    foreach_kplex_pivot(k - 1, pivot_candidate, [&](NodeSet &pivot_set) {
			    NodeSet candidate_plex;

			    pivot_set.insert(pivot_node);

			    set_union(pivot_set.begin(), pivot_set.end(), candidate.begin(), candidate.end(),
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

				    if (g.isKplex(plex, k)) {
					    screening_candidates.insert(plex);
				    }

			    } else {
				    search_graph = g.buildComplement(candidate_plex);

				    NodeSetSet bdd_sets = min_bdd_d_set(search_graph, bdd_max_del, k - 1, candidate);
				    NodeSet plex;

				    for (const NodeSet &bdd_set : bdd_sets) {
					    plex = candidate_plex;

					    for (NodeId deletion : bdd_set) {
						    plex.erase(deletion);
					    }

					    screening_candidates.insert(plex);
				    }
			    }

			    for (const NodeSet &plex : screening_candidates) {
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
					    spdlog::trace("Ignoring k-plex - failed pivot vertex check (rule #1).");
				    }
			    }

		    next_kplex:;
		    });

	    next_pivot:;
	    },
	    parallelism);

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