#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <iterator>
#include <string>
#include <fstream>

#include "libraries.h"
#include "parameters.h"

#include "vertex.h"
#include "edge.h"
#include "cell.h"


void getInitialData(VD& vd, std::vector<Cell>& cells, std::vector<Edge>& edges, std::vector<Vertex>& all_vertices);

void runSimulation(std::vector<Cell>& cells, std::vector<Edge>& edges, std::vector<Vertex>& vertices, int time_steps);


void outputData(const std::vector<Vertex>& vertices, const std::vector<Cell>& cells);

void WriteVTKFile(const std::vector<Vertex>& vertices, const std::vector<Edge>& edges, const std::string& filename);

#endif // FUNCTIONS_H


