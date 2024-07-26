#include "globals.h"

std::unordered_map<int, Vertex> vertex_map; int vertex_counter = 0;
std::unordered_map<int, Edge> edge_map; int edge_counter = 0;
std::unordered_map<int, Cell> cell_map; int cell_counter = 0;


//parameters
const int cell_count = 1000;

const double dt = 1e-10;
const int timesteps = 1;

const double A_0 = 0.001;

const double k_A = 0.01;
const double k_L = 0.5;
const double T_l_0 = 0.05;

const double a = 1;
