#include <conf.hpp>
#include <utils.hpp>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::pair;
using std::unordered_map;
using std::unordered_set;
using std::vector;

bool parallelism = false;

SGraph load_sgraph(string path) {
	ifstream in(path);
	vector<SEdge> edge_list;
	SEdge e;

	while (in >> e.nodeFrom >> e.nodeTo) {
		edge_list.push_back(e);
	}

	return SGraph(edge_list);
}

TGraph load_tgraph(string path) {
#if 0
	ifstream in(path);
	vector<TEdge> edge_list;
	TEdge e;

	while (in >> e.nodeFrom >> e.nodeTo >> e.tStart >> e.tStop) {
		edge_list.push_back(e);
	}

	return TGraph(edge_list);
#endif
	return load_tgraph(path, false, 0, 1);
}

TGraph load_tgraph(string path, bool squash, NodeTime sliding_window, NodeTime downsample) {
	ifstream in(path);
	vector<TEdge> edge_list;
	unordered_map<pair<NodeId, NodeId>, vector<NodeTime>, boost::hash<pair<NodeId, NodeId>>> edge_map;
	set<NodeTime> timestamps;
	unordered_map<NodeTime, NodeTime> timestamps_map;

	NodeId nodeFrom, nodeTo;
	NodeTime timestamp;

	while (in >> timestamp >> nodeFrom >> nodeTo) {
		edge_map[pair<NodeId, NodeId>(nodeFrom, nodeTo)].push_back(timestamp);
		timestamps.insert(timestamp);
	}

	if (squash) {
		NodeTime next_timestamp = 0;
		for (NodeTime t : timestamps) {
			timestamps_map[t] = next_timestamp++;
		}

		for (auto &val : edge_map) {
			for (auto it = val.second.begin(); it < val.second.end(); it++) {
				*it = timestamps_map[*it];
			}
		}
	}

	for (auto &val : edge_map) {
		for (auto it = val.second.begin(); it < val.second.end(); it++) {
			*it /= downsample;
		}
		std::sort(val.second.begin(), val.second.end());

		NodeTime start_time = val.second[0], last_time = val.second[0];

		for (unsigned int i = 1; i < val.second.size(); i++) {
			if (last_time + sliding_window < val.second[i] - 1) {
				edge_list.push_back(TEdge(val.first.first, val.first.second, start_time, last_time));
				start_time = last_time = val.second[i];
			} else {
				last_time = val.second[i];
			}
		}

		edge_list.push_back(TEdge(val.first.first, val.first.second, start_time, last_time));
	}

	return TGraph(edge_list);
}