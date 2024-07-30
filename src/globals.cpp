#include "globals.h"

std::unordered_map<int, Vertex> vertex_map; int vertex_counter = 0;
std::unordered_map<int, Edge> edge_map; int edge_counter = 0;
std::unordered_map<int, Cell> cell_map; int cell_counter = 0;


//parameters
const int cell_count = 1000;

const double dt = 1e-7;
const int timesteps = 1000;

const double A_0 = 1/cell_count;

const double k_A = 1;
const double k_L = 1;
const double T_l_0 = 0.25;

const double a = 1;
