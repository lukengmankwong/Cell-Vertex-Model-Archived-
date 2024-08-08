#include "global.h"

//parameters
const int cell_count = 1000;

const double dt = 1e-6;
const int timesteps = 400;

const double A_0 = 1.0/cell_count;

const double k_A = 1;
const double k_L = 1;
const double T_l_0 = 1;

const double a = 1;


//Global class

Global::Global() : vertex_counter(0), edge_counter(0), cell_counter(0) {}
Global::~Global() {}

std::unordered_map<int, Vertex>& Global::vertexMap() { return vertex_map; }
std::unordered_map<int, Edge>& Global::edgeMap() { return edge_map; }
std::unordered_map<int, Cell>& Global::cellMap() { return cell_map; }
const int Global::vertexCounter() const { return vertex_counter; }
const int Global::edgeCounter() const { return edge_counter; }
const int Global::cellCounter() const { return cell_counter; }


void Global::createVertex(Point r)
{
	vertex_map.emplace(vertex_counter, Vertex(this, vertex_counter, r));
	vertex_counter++;
}
const int Global::createEdge(int v1, int v2)
{
	edge_map.emplace(edge_counter, Edge(this, edge_counter, v1, v2));
	vertex_map.at(v1).addEdgeContact(edge_counter); //vertex v1 knows it's part of edge
	vertex_map.at(v2).addEdgeContact(edge_counter); //vertex v2 knows it's part of edge
	edge_counter++;
	return edge_counter-1; //return id of the created edge
}
const int Global::createCell(std::vector<int>& vertex_keys, std::vector<int>& edge_keys)
{
	cell_map.emplace(cell_counter, Cell(this, cell_counter, vertex_keys, edge_keys));
	for (int v : vertex_keys) { vertex_map.at(v).addCellContact(cell_counter); } //vertices know they are part of cell
	for (int e : edge_keys) { edge_map.at(e).addCellJunction(cell_counter); } //edges know they are part of cell
	cell_counter++;
	return cell_counter-1; //return id of created cell
}


void Global::cellNewVertex(int c, int v, int i) 
{ 
	cell_map.at(c).addVertex(v, i);  
	vertex_map.at(v).addCellContact(c);
}
void Global::cellNewEdge(int c, int e, int i)
{
	if (i == -1) { cell_map.at(c).addEdge(e, cell_map.at(c).getEdges().size()); } //add to back if no argument is given
	else { cell_map.at(c).addEdge(e, i); }
	edge_map.at(c).addCellJunction(c);
}


void Global::destroyVertex(int v)
{
	vertex_map.erase(v);
}
void Global::destroyEdge(int e)
{
	edge_map.erase(e);
}
void Global::destroyCell(int c)
{
	for (int v : cell_map.at(c).getVertices()) { vertex_map.at(v).removeCellContact(c); }
	for (int e : cell_map.at(c).getEdges()) { edge_map.at(e).removeCellJunction(c); }
	cell_map.erase(c);
}

