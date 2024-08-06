#include "global.h"

std::unordered_map<int, Vertex> vertex_map; int vertex_counter = 0;
std::unordered_map<int, Edge> edge_map; int edge_counter = 0;
std::unordered_map<int, Cell> cell_map; int cell_counter = 0;


//parameters
const int cell_count = 1000;

const double dt = 1e-6;
const int timesteps = 400;

const double A_0 = 1.0/cell_count;

const double k_A = 10000;
const double k_L = 1;
const double T_l_0 = 1;

const double a = 1;

Global::Global() : vertex_counter(0), edge_counter(0), cell_counter(0) {}
Global::~Global() {}


const int Global::vertexCounter() const { return cell_counter; }


void Global::createVertex(Point r)
{
	this->vertex_map.emplace(vertex_counter, Vertex(vertex_counter, r));
	vertex_counter++;
}

void Global::createEdge(int v1, int v2)
{
	this->edge_map.emplace(edge_counter, Edge(edge_counter, v1, v2));
	edge_counter++;
}

void Global::createCell(std::vector<int>& vertex_keys, std::vector<int>& edge_keys)
{
	this->cell_map.emplace(cell_counter, Cell(cell_counter, vertex_keys, edge_keys));
	cell_counter++;
}


void Global::destroyVertex(int v)
{
	this->vertex_map.erase(v);
}

void Global::destroyEdge(int e)
{
	this->edge_map.erase(e);
}

void Global::destroyCell(int c)
{
	this->cell_map.erase(c);
}

