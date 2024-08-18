#include "global.h"

//constants

const double pi = 3.14159265358979323846264338;

//parameters
const int cell_count = 1000;

const double dt = 1e-8;
const int timesteps = 1200;

const double A_0 = 1.0/cell_count;

const double k_A = 10;
const double k_L = 2;
const double T_l_0 = 1;

const double a = 1;


//Global class

Global::Global() : v_c_(0), e_c_(0), c_c_(0), timestep(0) { }
Global::~Global() {}

void Global::nextStep() { timestep++; }
const int Global::Step() const { return timestep; }

void Global::addDefects()
{
	std::vector<int> cell_step_defects;
	for (const auto& cell : c_map) { if (std::fabs(cell.second.m()-0.5) < 1e-3 || std::fabs(cell.second.m()+0.5) < 1e-3 
		|| std::fabs(cell.second.m()-1) < 1e-3 || std::fabs(cell.second.m()+1) < 1e-3) cell_step_defects.push_back(cell.first); }
	cell_defects.push_back(cell_step_defects);
	
	std::vector<int> vertex_step_defects;
	for (const auto& vertex : v_map) { if (std::fabs(vertex.second.m()-0.5) < 1e-3 || std::fabs(vertex.second.m()+0.5) < 1e-3 
		|| std::fabs(vertex.second.m()-1) < 1e-3 || std::fabs(vertex.second.m()+1) < 1e-3) vertex_step_defects.push_back(vertex.first); }
	vertex_defects.push_back(vertex_step_defects);
}

const std::vector<int>& Global::cellStepDefects(int step) const { return cell_defects[step]; }

std::unordered_map<int, Vertex>& Global::vertexMap() { return v_map; }
std::unordered_map<int, Edge>& Global::edgeMap() { return e_map; }
std::unordered_map<int, Cell>& Global::cellMap() { return c_map; }

Vertex& Global::vert(int v) { return v_map.at(v); }
Edge& Global::edge(int e) { return e_map.at(e); }
Cell& Global::cell(int c) { return c_map.at(c); }

const int Global::v_c() const { return v_c_; }
const int Global::e_c() const { return e_c_; }
const int Global::c_c() const { return c_c_; }


const int Global::createVertex(Point r)
{
	v_map.emplace(v_c_, Vertex(this, r));
	return v_c_++; 																	//return id of created vertex
}
const int Global::createEdge(const int v1, const int v2)
{
	e_map.emplace(e_c_, Edge(this, v1, v2));
	v_map.at(v1).addEdgeContact(e_c_); 												//vertex v1 knows it's part of edge
	v_map.at(v2).addEdgeContact(e_c_); 												//vertex v2 knows it's part of edge
	return e_c_++; 																	//return id of created edge
}
const int Global::createCell(std::vector<int>& vertices, std::vector<std::pair<int,int>>& edges)
{
	c_map.emplace(c_c_, Cell(this, vertices, edges));
	for (int v : vertices) { v_map.at(v).addCellContact(c_c_); } 					//vertices know they are part of cell
	for (std::pair<int,int> e : edges) { e_map.at(e.first).addCellJunction(c_c_); } 	//edges know they are part of cell
	return c_c_++; 																	//return id of created cell
}

void Global::cellExchangeVertex(int c, int v_old, int v_new) { c_map.at(c).exchangeVertex(v_old, v_new); }

void Global::cellNewVertex(int c, int v, int i) 
{ 
	c_map.at(c).addVertex(v, i);  
	v_map.at(v).addCellContact(c);
}
void Global::cellNewEdge(int c, int e, int i)
{
	c_map.at(c).addEdge(e, i);
	e_map.at(e).addCellJunction(c);
}

void Global::cellRemoveVertex(int c, int v) { c_map.at(c).removeVertex(v); }
int Global::cellRemoveEdge(int c, int e) { return c_map.at(c).removeEdge(e); }


void Global::destroyVertex(int v) { v_map.erase(v);}
void Global::destroyEdge(int e) 
{ 
	for (int c : e_map.at(e).cellJunctions()) { c_map.at(c).removeEdge(e); }
	v_map.at(e_map.at(e).v1()).removeEdgeContact(e);
	v_map.at(e_map.at(e).v2()).removeEdgeContact(e);
	e_map.erase(e); 
}
void Global::destroyCell(int c)
{
	for (int v : c_map.at(c).Vertices()) { v_map.at(v).removeCellContact(c); }
	for (int e : c_map.at(c).Edges()) { e_map.at(e).removeCellJunction(c); }
	c_map.erase(c);
}

const bool Global::commonEdge(int c1, int c2) const
{
	for (int e : c_map.at(c1).Edges())
	{
		auto it = std::find(c_map.at(c2).Edges().begin(), c_map.at(c2).Edges().end(), e);
		if (it != c_map.at(c2).Edges().end()) { return true; }
	}
	return false;
}

void Global::extrusion()
{
	std::vector<int> small_cells;
	for (const auto& cell : c_map)
	{
		if (cell.second.A() < 0.2*A_0)
		{
			const std::vector<int>& vertices = c_map.at(cell.first).Vertices();
			int i = 0; bool contact = false;
			while (i < vertices.size() && contact == false)
			{
				for (int c : v_map.at(vertices[i]).cellContacts())
				{
					auto it = std::find(small_cells.begin(), small_cells.end(), c);
					if (it != small_cells.end()) { contact = true; }
				}
				i++;
			}
			if (!contact) small_cells.push_back(cell.first);
		}
	}
	for (int c : small_cells) { c_map.at(c).extrude(); }
}

void Global::division()
{
	std::vector<int> large_cells;
	for (const auto& cell : c_map)
	{
		if (cell.second.A() > 2*A_0)
		{
			const std::vector<int>& vertices = c_map.at(cell.first).Vertices();
			int i = 0; bool contact = false;
			while (i < vertices.size() && contact == false)
			{
				for (int c : v_map.at(vertices[i]).cellContacts())
				{
					auto it = std::find(large_cells.begin(), large_cells.end(), c);
					if (it != large_cells.end()) { contact = true; }
				}
				i++;
			}
			if (!contact) large_cells.push_back(cell.first);
		}
	}
	for (int c : large_cells) { c_map.at(c).divide(); }
}

void Global::T1()
{
	std::vector<int> fourfold_vertices;
	for (const auto& vertex : v_map)
	{
		if (vertex.second.edgeContacts().size() == 4) fourfold_vertices.push_back(vertex.first);
	}
	for (int v : fourfold_vertices) v_map.at(v).T1();
}

void Global::transitions()
{
	for (std::pair<const int,Cell>& cell : c_map) cell.second.calcA();
	extrusion();
	for (std::pair<const int,Cell>& cell : c_map) cell.second.calcA();
	division();
	//T1();
}

const double Global::D_angle(int c_i, int c_j) const
{
	double Z_i = c_map.at(c_i).Z(); double X_i = c_map.at(c_i).X(); double S_i = std::sqrt(X_i*X_i+Z_i*Z_i);
	double Z_j = c_map.at(c_j).Z(); double X_j = c_map.at(c_j).X(); double S_j = std::sqrt(X_j*X_j+Z_j*Z_j);
	return std::asin( (Z_i*X_j-X_i*Z_j) / (S_i*S_j) );
}


