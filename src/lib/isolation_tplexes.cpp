#include <isolation_tplexes.hpp>

#include <spdlog/spdlog.h>
#include <unordered_map>

#include <conf.hpp>
#include <isolation_splexes.hpp>

#include <iostream>
using namespace std;

using std::max;
using std::min;
using std::set_difference;
using std::sort;
using std::unordered_map;

NodeSetIntervalSet c_isolated_temporal_kplex(TGraph &g, int k, int c, TemporalIsolationType isolation) {
	NodeSetIntervalSet result;

	unordered_map<Interval, NodeSetSet, boost::hash<Interval>> interval_map;
	unordered_map<NodeSet, set<Interval>, boost::hash<NodeSet>> nodeset_map;

	spdlog::info("Starting c_isolated_temporal_kplex; is parallelism enabled? {}", parallelism);
	spdlog::info("TGraph lifetime: [{}, {}]", g.getLifetimeBegin(), g.getLifetimeEnd());

#pragma omp parallel for if (parallelism)
	for (NodeTime i = g.getLifetimeBegin(); i <= g.getLifetimeEnd(); i++) {
		SGraph gg = g.buildIntersectionGraph(g.getNodes(), i, i);

		NodeSetSet res;
		switch (isolation) {
		case ALLTIME_MAX:
		case USUALLY_MAX:
		case MAX_USUALLY:
			res = max_c_isolated_kplex(gg, c, k);
			break;
		default:
			res = avg_c_isolated_kplex(gg, c, k);
		}

		for (const NodeSet &s : res) {
#pragma omp critical(interval_map)
			{
				interval_map[Interval(i, i)].insert(s);
				if (!g.isKplex(s, k, i, i)) {
					spdlog::error("Instant {}, candidate {} is NOT a {}-plex", i,
						      nodeset_to_string(s), k);
				}
			}
		}
		if (res.size() == 0) {
			spdlog::debug("Instant {}, no candidates found.", i);
		}
	}

	spdlog::info("c_isolated_temporal_kplex: initialization done");

	for (NodeTime len = 2; len <= g.getLifetimeEnd() - g.getLifetimeBegin() + 1; len++) {
		NodeTime lifetime_end = g.getLifetimeEnd();
#pragma omp parallel for if (parallelism)
		for (NodeTime begin_w = 0; begin_w <= lifetime_end - len + 1; begin_w++) {
			NodeTime end_w = begin_w + len - 1;
			for (int i = 0; i < 2; i++) {
				NodeTime begin = begin_w - i + 1;
				NodeTime end = end_w - i;
				NodeTime crit = (i == 1 ? end_w : begin_w);

				if (interval_map.find(Interval(begin, end)) == interval_map.end()) {
					spdlog::trace("No candidate for interval [{}, {}]", begin, end);
				} else {
					for (const NodeSet &candidate : interval_map[Interval(begin, end)]) {
						SGraph g_star = g.buildAuxGraph(candidate, begin_w, end_w, crit);
						NodeSetSet candidate_k_set;
						switch (isolation) {
						case ALLTIME_MAX:
						case USUALLY_MAX:
						case MAX_USUALLY:
							candidate_k_set =
							    max_c_isolated_kplex_restricted(g_star, c, k, candidate);
							break;
						default:
							candidate_k_set =
							    avg_c_isolated_kplex_restricted(g_star, c, k, candidate);
							break;
						}

						for (const NodeSet &candidate_k : candidate_k_set) {
#pragma omp critical(interval_map)
							{
								interval_map[Interval(begin_w, end_w)].insert(
								    candidate_k);

								if (!g.isKplex(candidate_k, k, begin_w, end_w)) {
									spdlog::error("{} is not a {}-plex in [{}, {}]",
										      nodeset_to_string(candidate_k), k,
										      begin_w, end_w);
								}
							}
							NodeSetSet isolated_subsets;

							switch (isolation) {
							case ALLTIME_MAX:
								isolated_subsets = alltime_max_isolated_subset(
								    g, candidate_k, k, c, begin_w, end_w,
								    g_star.mindegree(candidate_k));
								break;
							case ALLTIME_AVG:
								isolated_subsets = alltime_avg_isolated_subset(
								    g, candidate_k, k, c, begin_w, end_w,
								    g_star.mindegree(candidate_k));
								break;
							case USUALLY_MAX:
								isolated_subsets = usually_max_isolated_subset(
								    g, candidate_k, k, c, begin_w, end_w,
								    g_star.mindegree(candidate_k));
								break;
							case USUALLY_AVG:
								isolated_subsets = usually_avg_isolated_subset(
								    g, candidate_k, k, c, begin_w, end_w,
								    g_star.mindegree(candidate_k));
								break;

							case AVG_ALLTIME:
								isolated_subsets = avg_alltime_isolated_subset(
								    g, candidate_k, k, c, begin_w, end_w,
								    g_star.mindegree(candidate_k));
								break;

							case MAX_USUALLY:
								isolated_subsets = max_usually_isolated_subset(
								    g, candidate_k, k, c, begin_w, end_w,
								    g_star.mindegree(candidate_k));
								break;
							}

#pragma omp critical(nodeset_map)
							{
								spdlog::debug("Found {} isolated subsets ({}).",
									      isolated_subsets.size(),
									      nodesetset_to_string(isolated_subsets));
								for (const NodeSet &isolated : isolated_subsets) {
									nodeset_map[isolated].insert(
									    Interval(begin_w, end_w));
									nodeset_map[isolated].erase(
									    Interval(begin, end));
								}
							}
						}
					}
				}
			}
		}

		spdlog::info("c_isolated_temporal_kplex enumeration: [{}/{}] interval length done", len,
			     g.getLifetimeEnd() - g.getLifetimeBegin() + 1);
	}

	spdlog::info("c_isolated_temporal_kplex: enumeration done");

	for (const auto &r : nodeset_map) {
		set<Interval> intervals = r.second;
		for (Interval i1 : r.second) {
			for (Interval i2 : r.second) {
				if (i1.first > i2.first || i1.second < i2.second) {
					intervals.erase(i1);
					break;
				}
			}
		}

		for (Interval i : intervals) {
			result.insert(NodeSetInterval(r.first, i));
		}
	}

	spdlog::info("c_isolated_temporal_kplex: maximality check done");

	return result;
}

NodeSetSet alltime_max_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start,
				       NodeTime t_stop, int delta) {
	NodeSet candidate, candidate_prime;
	candidate = candidate_prime = nodeset;
	NodeSetSet result;

	// int minsize = max(delta - c - k + 1, 0);
	int minsize = max(delta - c - k + 1, k + 2);
	bool fixpoint;

	do {
		fixpoint = true;

		for (NodeId v : candidate_prime) {
			if (g.outdegree_max(v, candidate, t_start, t_stop) >= c) {
				candidate.erase(v);
				fixpoint = false;
			}

			if ((int)(candidate.size()) < minsize) {
				fixpoint = false;
				goto break_outer;
			}
		}
		candidate_prime = candidate;
	} while (!fixpoint);
break_outer:;

	if (fixpoint && (int)candidate.size() >= minsize) {
		result.insert(candidate);
	}

	return result;
}

NodeSetSet max_usually_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start,
				       NodeTime t_stop, int delta) {
	NodeSetSet result;

	/* TODO */

	return result;
}

NodeSetSet avg_alltime_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start,
				       NodeTime t_stop, int delta) {
	NodeSetSet result;

	/* TODO */

	return result;
}

NodeSetSet usually_avg_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start,
				       NodeTime t_stop, int delta) {
	NodeSetSet result;
	vector<pair<NodeId, int>> degrees;

	int min_plex_size = max(delta - c - k, k + 2);

	if ((int)nodeset.size() < min_plex_size) {
		return result;
	}

	for (NodeId u : nodeset) {
		degrees.push_back(pair<NodeId, int>(u, g.degree_sum(u, t_start, t_stop)));
	}

	sort(degrees.begin(), degrees.end(), [&](const pair<NodeId, int> &rhs, const pair<NodeId, int> &lhs) {
		if (rhs.second == lhs.second) {
			return rhs.first > lhs.first;
		} else {
			return rhs.second > lhs.second;
		}
	});

	int del_indices = min(max(0, (int)(nodeset.size()) - delta + c + 3 * k), (int)(nodeset.size()));

	NodeSetSet deletion_sets, deletion_sets_prime;

	deletion_sets_prime.insert(NodeSet());
	while (!deletion_sets_prime.empty()) {
		deletion_sets = deletion_sets_prime;
		deletion_sets_prime.clear();

		for (const NodeSet &deletion : deletion_sets) {
			NodeSet candidate;
			set_difference(nodeset.begin(), nodeset.end(), deletion.begin(), deletion.end(),
				       std::inserter(candidate, candidate.begin()));

			int condition = 0;
			for (NodeId u : candidate) {
				condition += g.outdegree_time_sum(u, candidate, t_start, t_stop);
			}

			if (condition >= (int)(candidate.size()) * (c - 1) * (t_stop - t_start + 1)) {
				if ((int)candidate.size() > min_plex_size) {
					int j = 0;
					for (int i = 0; i < del_indices; i++) {
						if (deletion.find(degrees[i].first) != deletion.end()) {
							j = i;
						}
					}

					for (int i = j + 1; i < del_indices; i++) {
						NodeSet deletion_new = deletion;
						deletion_new.insert(degrees[i].first);
						deletion_sets_prime.insert(deletion_new);
					}
				}
			} else {
				/* It is isolated */
				result.insert(candidate);
			}
		}
	}

	return result;
}

NodeSetSet alltime_avg_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start,
				       NodeTime t_stop, int delta) {
	NodeSetSet result;
	int del_indices = min(max(0, (int)(nodeset.size()) - delta + c + 3 * k), (int)(nodeset.size()));
	int min_plex_size = max(delta - c - k, k + 2);

	if ((int)nodeset.size() < min_plex_size) {
		return result;
	}

	NodeSetSet deletions, deletions_prime;
	deletions_prime.insert(NodeSet());

	while (!deletions_prime.empty()) {
		deletions = deletions_prime;
		deletions_prime.clear();
		for (const NodeSet &deletion : deletions) {
			NodeSet candidate;
			set_difference(nodeset.begin(), nodeset.end(), deletion.begin(), deletion.end(),
				       std::inserter(candidate, candidate.begin()));

			for (NodeTime i = t_start; i <= t_stop; i++) {
				int condition = 0;
				for (NodeId v : candidate) {
					condition += g.outdegree(v, candidate, i);
				}

				if (condition >= (int)candidate.size() * c) {
					if ((int)candidate.size() > min_plex_size) {
						vector<pair<NodeId, int>> degrees;

						for (NodeId u : nodeset) {
							degrees.push_back(
							    pair<NodeId, int>(u, g.outdegree(u, candidate, i)));
						}

						sort(degrees.begin(), degrees.end(),
						     [&](const pair<NodeId, int> &rhs, const pair<NodeId, int> &lhs) {
							     if (rhs.second == lhs.second) {
								     return rhs.first > lhs.first;
							     } else {
								     return rhs.second > lhs.second;
							     }
						     });

						for (int i = 0; i < del_indices; i++) {
							if (deletion.find(degrees[i].first) == deletion.end()) {
								NodeSet deletion_new = deletion;
								deletion_new.insert(degrees[i].first);
								deletions_prime.insert(deletion_new);
							}
						}
					}

					goto next_deletion_set;
				}
			}

			result.insert(candidate);

		next_deletion_set:;
		}
	}

	return result;
}

NodeSetSet usually_max_isolated_subset(TGraph &g, const NodeSet &nodeset, int k, int c, NodeTime t_start,
				       NodeTime t_stop, int delta) {
	NodeSetSet result;
	int minsize = max(delta - c - k + 1, k + 2);
	if ((int)nodeset.size() < minsize) {
		return result;
	}

	vector<pair<NodeId, int>> degrees;

	for (NodeId u : nodeset) {
		degrees.push_back(pair<NodeId, int>(u, g.degree_sum(u, t_start, t_stop)));
	}

	sort(degrees.begin(), degrees.end(), [&](const pair<NodeId, int> &rhs, const pair<NodeId, int> &lhs) {
		if (rhs.second == lhs.second) {
			return rhs.first > lhs.first;
		} else {
			return rhs.second > lhs.second;
		}
	});
	NodeSetSet deletions, deletions_prime;

	deletions_prime.insert(NodeSet());

	while (!deletions_prime.empty()) {
		deletions = deletions_prime;
		deletions_prime.clear();

		for (const NodeSet &deletion : deletions) {
			NodeSet candidate;
			unordered_set<NodeId> max_vertices;

			set_difference(nodeset.begin(), nodeset.end(), deletion.begin(), deletion.end(),
				       std::inserter(candidate, candidate.begin()));

			int condition = 0;
			for (NodeTime t = t_start; t <= t_stop; t++) {
				int max_outdeg = -1;
				NodeId max_v = -1;
				for (NodeId v : candidate) {
					int outdeg = g.outdegree(v, candidate, t);
					if (outdeg > max_outdeg) {
						max_outdeg = outdeg;
						max_v = v;
					}
				}

				condition += max_outdeg;
				max_vertices.insert(max_v);
			}

			if (condition >= (t_stop - t_start + 1) * c) {
				if ((int)candidate.size() > minsize) {
					for (NodeId v : max_vertices) {
						if (deletion.find(v) == deletion.end()) {
							NodeSet deletion_new = deletion;
							deletion_new.insert(v);

							deletions_prime.insert(deletion_new);
						}
					}
				}
			} else {
				result.insert(candidate);
			}
		}
	}

	return result;
}
