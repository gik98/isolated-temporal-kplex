#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "Graph.hpp"
#include <fstream>
#include <string>

using std::ifstream;
using std::string;

SGraph load_sgraph(string path);

TGraph load_tgraph(string path);

TGraph load_tgraph(string path, bool squash, NodeTime sliding_window, NodeTime downsample);

#endif