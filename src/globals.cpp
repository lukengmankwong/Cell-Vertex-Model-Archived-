#include "globals.h"

std::unordered_map<int, Vertex> vertex_map; int vertex_counter = 0;
std::unordered_map<int, Edge> edge_map; int edge_counter = 0;
std::unordered_map<int, Cell> cell_map; int cell_counter = 0;


//parameters
const int cell_count = 2500;

const double dt = 0.001;
const int timesteps = 500;

const double A_0 = 0.001;

const double a = 1;
const double k = 1;
