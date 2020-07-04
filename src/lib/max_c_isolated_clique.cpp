#include <isolation_scliques.hpp>

#include <algorithm>
#include <unordered_map>

#include <conf.hpp>

using std::set_intersection;
using std::unordered_map;

NodeSetSet max_c_isolated_clique(SGraph &g, int c) {
	NodeSetSet sol;

	unordered_map<NodeSet, NodeId, boost::hash<NodeSet>> screening_candidates;
	/* Forall pivots, we enumerate the max_c_isolated cliques */
	g.forallNodes(
	    [&](NodeId pivot) {
		    /* Candidate set */
		    NodeSet candidate = g.neighbourhoodGreaterDeg(pivot);
		    candidate.insert(pivot);

		    SGraph search_graph;
		    NodeSetSet enumeration_delete_sets;

		    /* Trimming stage */
		    bool fixpoint;
		    int d = 0;
		    int max_vc_size;
		    do {
			    if (d >= c) {
				    /* We removed too many vertices, the pivot has c outgoing deg */
				    goto next_pivot;
			    }
			    fixpoint = true;
			    for (NodeId node : candidate) {
				    /* Invariant (a) */
				    if (g.degree(node) >= (int)(candidate.size()) + c - 1) {
					    candidate.erase(node);
					    d++;
					    fixpoint = false;
					    break;
				    }

				    /* Invariant (b) */
				    if (g.degree(node, candidate) < (int)(candidate.size()) - c) {
					    candidate.erase(node);
					    d++;
					    fixpoint = false;
					    break;
				    }
			    }
		    } while (!fixpoint);

		    /* Enumeration stage */
		    search_graph = g.buildComplement(candidate);
		    max_vc_size = c - d - 1;

		    enumeration_delete_sets = min_vertex_cover_bounded(search_graph, max_vc_size);

		    for (const NodeSet &dset : enumeration_delete_sets) {
			    int dd = d;
			    NodeSet diffset = candidate;
			    for (NodeId u : dset) {
				    diffset.erase(u);
				    dd++;
			    }

			    /* Check if the resulting clique respects the max-c-isolation condition */
			    bool fixpoint;
			    do {
				    fixpoint = true;
				    if (dd >= c) {
					    /* We removed too many nodes - this clique is not max-c */
					    goto next_candidate;
				    }

				    for (NodeId u : diffset) {
					    if (g.outdegree(u, diffset) >= c) {
						    /* This node violates max-c */
						    diffset.erase(u);
						    dd++;
						    fixpoint = false;
						    break;
					    }
				    }
			    } while (!fixpoint);

#pragma omp critical(screening_candidates)
			    { screening_candidates[diffset] = pivot; }

		    next_candidate:;
		    }
	    next_pivot:;
	    },
	    parallelism);

	/* Screening stage */
#pragma omp parallel if (parallelism)
	{
#pragma omp single
		{
			for (const auto &c1 : screening_candidates) {
				NodeSet vneigh = g.neighbourhoodSmallerDeg(c1.second);
				bool maximal = true;
				/* Lemma 5 */
				for (const auto &c2 : screening_candidates) {
					if (!maximal) {
						break;
					}
					if (c1.second != c2.second && vneigh.find(c2.second) != vneigh.end()) {
						NodeSet uneigh = g.neighbourhood(c2.second);
						NodeSet result;
						set_intersection(vneigh.begin(), vneigh.end(), uneigh.begin(),
								 uneigh.end(), std::inserter(result, result.begin()));

						if (result == c1.first) {
							maximal = false;
						}
					}
				}
#pragma omp critical(sol)
				if (maximal) {
					sol.insert(c1.first);
				}
			}
		}
	}
	return sol;
}