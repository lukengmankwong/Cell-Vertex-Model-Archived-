#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <unordered_map>
#include <string>
#include <fstream>

#include "libraries.h"
#include "global.h"

#include "vertex.h"
#include "edge.h"
#include "cell.h"

void getInitialData(VD& vd, Global& global, bool (*in)(const Point&));
void runSimulation(Global& global, int time_steps);
void outputData(Global& global);
void WriteVTKFile(Global& global, const std::string& filename_graph, const std::string& filename_cell_defect, const std::string& filename_vertex_defect, const std::string& filename_director);

#endif // FUNCTIONS_H
