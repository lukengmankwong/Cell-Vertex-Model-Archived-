#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <unordered_map>
#include <functional>
#include <string>
#include <fstream>

#include "libraries.h"
#include "globals.h"

#include "vertex.h"
#include "edge.h"
#include "cell.h"

void getInitialData(VD& vd, bool (*in)(const Point&));

void runSimulation(int time_steps);

void outputData(const std::unordered_map<int, Vertex>& vertex_map, const std::unordered_map<int, Edge>& edge_map, const std::unordered_map<int, Cell>& cell_map);

void WriteVTKFile(const std::unordered_map<int, Vertex>& vertex_map, const std::unordered_map<int, Edge>& edge_map, const std::unordered_map<int, Cell>& cell_map, const std::string& filename_graph, const std::string& filename_director);

#endif // FUNCTIONS_H
