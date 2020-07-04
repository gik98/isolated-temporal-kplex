#include <ctype.h>
#include <fstream>
#include <iostream>
#include <spdlog/spdlog.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>

#include <Graph.hpp>
#include <conf.hpp>
#include <isolation_scliques.hpp>
#include <isolation_splexes.hpp>
#include <isolation_tplexes.hpp>
#include <utils.hpp>

using std::ofstream;
using std::string;

const string help_str =
    "Argument list:\n-d <path>:\t Load the specified dataset\n-c <isolation factor>:\t Set the isolation factor "
    "(int)\n-k <relaxation>:\t Set for maximal isolated k-plexes, where k is the relaxation parameter. Either set k or "
    "pass the -C flag.\n-C:\t Search for maximal isolated cliques. Either pass -C or set the relaxation parameter "
    "k.\n-m:\tSearch min-c-isolated communities\n-M:\tSearch max-c-isolated communities\n-a:\tSearch avg-c-isolated "
    "communities\n-p:\tEnable parallelism\n-v:\tVerbose logging\n-V:\tVery verbose logging\n-T <algorithm>: temporal "
    "graph analysis\n-s:\tSquash temporal dataset\n-D <n>: downsample temporal dataset\n-w <n>: sliding window for "
    "temporal dataset\n-o <path>: output file\n-X:\t Test output correctness";

int main(int argc, char **argv) {
	spdlog::set_level(spdlog::level::info);
	int opt;

	string dataset, output;
	int c, k;
	c = k = -1;

	bool min_c_isolation, max_c_isolation, avg_c_isolation, clique, dataset_set;
	min_c_isolation = max_c_isolation = avg_c_isolation = clique = dataset_set = false;

	bool temporal = false;
	bool squash = false;
	bool check = false;
	bool print_output = false;
	int downsample = 1;
	int sliding_window = 0;
	string temporal_algo;

	while ((opt = getopt(argc, argv, "d:c:k:mMaCpvVhT:D:sw:o:X")) != -1) {
		switch (opt) {
		case 'h':
			std::cout << help_str << std::endl;
			return 0;
		case 'a':
			avg_c_isolation = true;
			break;
		case 'm':
			min_c_isolation = true;
			break;
		case 'M':
			max_c_isolation = true;
			break;
		case 'c':
			c = atoi(optarg);
			break;
		case 'C':
			clique = true;
			break;
		case 'k':
			k = atoi(optarg);
			break;
		case 'p':
			parallelism = true;
			break;
		case 'd':
			dataset = string(optarg);
			dataset_set = true;
			break;
		case 'v':
			spdlog::set_level(spdlog::level::debug);
			break;
		case 'V':
			spdlog::set_level(spdlog::level::trace);
			break;
		case 'T':
			temporal_algo = string(optarg);
			temporal = true;
			break;
		case 's':
			squash = true;
			break;
		case 'D':
			downsample = atoi(optarg);
			break;
		case 'w':
			sliding_window = atoi(optarg);
			break;
		case 'o':
			output = string(optarg);
			print_output = true;
			break;
		case 'X':
			check = true;
			break;
		case '?':
			if (optopt == 'c' || optopt == 'd' || optopt == 'k') {
				spdlog::error("Option {} requires an argument!", optopt);
				return 1;
			} else if (isprint(optopt)) {
				spdlog::warn("Unknown option {}", optopt);
			} else {
				spdlog::warn("Unknown option character {}", optopt);
			}

			break;
		default:
			return 2;
		}
	}

	if (temporal) {
		NodeSetIntervalSet res;

		TGraph g = load_tgraph(dataset, squash, sliding_window, downsample);

#if 0
		auto adj = g.getAdjacencyList();
		for (auto p: adj) {
			std::cout << "Node " << p.first << ": ";
			for (TEdge edge: p.second) {
				std::cout << edge.nodeTo << " [" << edge.tStart << "," << edge.tStop << "];\t";
			}
			std::cout << std::endl;
		}
#endif

		spdlog::info("Loaded graph {}, {} nodes, {} edges, {} edges instants, lifetime: [{}, {}]", dataset,
			     g.getNodesCount(), g.getEdgesCount(), g.getEdgesInstantsCount(), g.getLifetimeBegin(),
			     g.getLifetimeEnd());

		TemporalIsolationType type;
		if (temporal_algo == "alltime-max") {
			type = ALLTIME_MAX;
		} else if (temporal_algo == "usually-avg") {
			type = USUALLY_AVG;
		} else if (temporal_algo == "alltime-avg") {
			type = ALLTIME_AVG;
		} else if (temporal_algo == "usually-max") {
			type = USUALLY_MAX;
		} else {
			spdlog::error("Isolation type {} not supported", temporal_algo);
			return 1;
		}
		auto begin = std::chrono::high_resolution_clock::now();
		res = c_isolated_temporal_kplex(g, k, c, type);
		auto end = std::chrono::high_resolution_clock::now();

		auto duration_us = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count() / 1000;

		spdlog::info("Graph {}, {}-{}-isolation returned {} {}-plexes. Took {} us", dataset,
			     temporal_algo, c, res.size(), k, duration_us);

		int i = 0;
		for (auto sol : res) {
			spdlog::info("k-plex #{}: {}", i++, nodesetinterval_to_string(sol));
		}

		bool check_errors = false;
		if (print_output) {
			ofstream out(output);
			// algo, c, D, k, w, N, E, E-instants, TIME
			out << temporal_algo << " " << c << " " << downsample << " " << k << " " << sliding_window
			    << " " << g.getNodesCount() << " " << g.getEdgesCount() << " " << g.getEdgesInstantsCount()
			    << " " << duration_us << std::endl;
			// #returned k-plex
			out << res.size() << std::endl;
			for (auto sol : res) {
				out << sol.second.first << " " << sol.second.second << " " << sol.first.size();
				for (NodeId u : sol.first) {
					out << " " << u;
				}
				out << std::endl;
			}

			out.close();
			spdlog::info("Solution has been written to output file {}.", output);
		}

		for (auto sol : res) {
			if (check) {
				if (!g.isKplex(sol.first, k, sol.second.first, sol.second.second)) {
					spdlog::error("A returned set does not satisfy the {}-plex condition: {}", k,
						      nodesetinterval_to_string(sol));
					check_errors = true;
				}
				if (!g.isIsolated(sol.first, sol.second.first, sol.second.second, type, c)) {
					spdlog::error("A returned k-plex does not satisfy the isolation condition: {}",
						      nodesetinterval_to_string(sol));
					check_errors = true;
				}
			}
		}

		if (check && check_errors) {
			spdlog::info("Solution is NOT correct");
		} else if (check && !check_errors) {
			spdlog::info("Solution is correct");
		}
	} else {
		if (!dataset_set) {
			spdlog::error("You must specify a dataset using the -d flag.");
			return 1;
		}

		if (c < 1) {
			spdlog::error("The isolation factor c is a mandatory flag and shall be a positive number");
			return 1;
		}

		if (k < 1 && !clique) {
			spdlog::error(
			    "Thek-plex paramter k is a mandatory flag if you did not set the clique flag -C.");
			return 1;
		}

		SGraph g = load_sgraph(dataset);
		spdlog::info("Input graph has {} nodes and {} edges", g.getNodesCount(), g.getEdgesCount());
		if (clique) {
			if (min_c_isolation) {
				spdlog::info("Running clique min-{}-isolation", c);
				auto res = min_c_isolated_clique(g, c);
				spdlog::info("{} maximal min-{}-isolated cliques found.", res.size(), c);

				int i = 1;
				for (const NodeSet &s : res) {
					spdlog::info("Clique #{}: {}", i++, nodeset_to_string(s));
				}
			}
			if (max_c_isolation) {
				spdlog::info("Running clique max-{}-isolation", c);
				auto res = max_c_isolated_clique(g, c);
				spdlog::info("{} maximal max-{}-isolated cliques found.", res.size(), c);

				int i = 1;
				for (const NodeSet &s : res) {
					spdlog::info("Clique #{}: {}", i++, nodeset_to_string(s));
				}
			}
			if (avg_c_isolation) {
				spdlog::error("clique avg-{}-isolation not supported yet.", c);
			}
		} else {
			if (min_c_isolation) {
				spdlog::info("Running {}-plex min-{}-isolation", k, c);
				auto res = min_c_isolated_kplex(g, c, k);
				spdlog::info("{} maximal min-{}-isolated {}-plexes found.", res.size(), c, k);

				int i = 1;
				for (const NodeSet &s : res) {
					spdlog::info("k-plex #{}: {}", i++, nodeset_to_string(s));
				}
			}
			if (max_c_isolation) {
				spdlog::info("Running {}-plex max-{}-isolation", k, c);
				auto res = max_c_isolated_kplex(g, c, k);
				spdlog::info("{} maximal max-{}-isolated {}-plex found.", res.size(), c, k);

				int i = 1;
				for (const NodeSet &s : res) {
					spdlog::info("k-plex #{}: {}", i++, nodeset_to_string(s));
				}
			}
			if (avg_c_isolation) {
				spdlog::info("Running {}-plex avg-{}-isolation", k, c);
				auto res = avg_c_isolated_kplex(g, c, k);
				spdlog::info("{} maximal avg-{}-isolated {}-plex found.", res.size(), c, k);

				int i = 1;
				for (const NodeSet &s : res) {
					spdlog::info("k-plex #{}: {}", i++, nodeset_to_string(s));
				}
			}
		}
	}

	return 0;
}
