#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>

#include "libraries.h"
#include "parameters.h"

#include "vertex.h"
#include "edge.h"
#include "cell.h"

void getInitialData(VD& vd, std::unordered_map<int, Vertex>& vertex_map, std::unordered_map<int, Edge>& edge_map, std::unordered_map<int, Cell>& cell_map, int& vertex_counter, int& edge_counter, int& cell_counter);

void runSimulation(std::unordered_map<int, Vertex>& vertex_map, std::unordered_map<int, Edge>& edge_map, std::unordered_map<int, Cell>& cell_map, int time_steps);

void outputData(const std::unordered_map<int, Vertex>& vertex_map, const std::unordered_map<int, Edge>& edge_map, const std::unordered_map<int, Cell>& cell_map);

void WriteVTKFile(const std::unordered_map<int, Vertex>& vertex_map, const std::unordered_map<int, Edge>& edge_map, const std::string& filename);

#endif // FUNCTIONS_H
