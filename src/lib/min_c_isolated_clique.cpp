#include <isolation_scliques.hpp>

#include <algorithm>
#include <unordered_map>

#include <conf.hpp>

using std::set_intersection;
using std::unordered_map;

NodeSetSet min_c_isolated_clique(SGraph &g, int c) {
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
				    /* Invariant (b) */
				    if (g.degree(node, candidate) < (int)(candidate.size()) - c) {
					    candidate.erase(node);
					    d++;
					    fixpoint = false;
					    break;
				    }
			    }
		    } while (!fixpoint);

		    /* The pivot must have less than c neighbours outside of the clique */
		    if (g.outdegree(pivot, candidate) >= c) {
			    goto next_pivot;
		    }

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
#pragma omp critical(screening_candidates)
			    { screening_candidates[diffset] = pivot; }
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
#pragma omp task if (parallelism)
				{
					NodeSet vneigh = g.neighbourhoodSmallerDeg(c1.second);
					bool maximal = true;
					/* If there is a vertex which is adjacent to all vertices of the clique, the
					 * clique is not maximal */
					for (const NodeId &u : vneigh) {
						NodeSet uneigh = g.neighbourhood(u);
						NodeSet result;
						set_intersection(c1.first.begin(), c1.first.end(), uneigh.begin(),
								 uneigh.end(), std::inserter(result, result.begin()));
						if (result.size() == c1.first.size()) {
							maximal = false;
							break;
						}
					}

					if (maximal) {
#pragma omp critical(sol)
						sol.insert(c1.first);
					}
				}
			}
		}
	}
	return sol;
}