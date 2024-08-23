#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <vector>
#include <unordered_map>

#include <string>
#include <fstream>

#include "libraries.h"
#include "tissue.h"
#include "vertex.h"
#include "edge.h"
#include "cell.h"

void getInitialData(VD& vd, Tissue& tissue, bool (*in)(const Point&));
void outputData(const Tissue& Tissue);

void writeCellsFile(Tissue* tissue, const std::string& filename_cells);
void writeDirectorsFile(Tissue* tissue, const std::string& filename_directors);
void writeCellDefectsFile(Tissue* tissue, const std::string& filename_cell_defects);
void writeVertexDefectsFile(Tissue* tissue, const std::string& filename_vertex_defects);

#endif // FUNCTIONS_H
