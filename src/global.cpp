#include "global.h"

//constants

const double pi = 3.14159265358979323846264338;

//parameters
const int cell_count = 5000;

const double dt = 1e-7;
const int timesteps = 1000;

const double A_0 = 1.0/cell_count;

const double k_A = 10;
const double k_L = 2;
const double T_l_0 = 0;

const double a = 3;


//Global class

Global::Global() : v_c(0), e_c(0), c_c(0), timestep(0) {}
Global::~Global() {}

void Global::nextStep() { timestep++; }
const int Global::Step() const { return timestep; }

void Global::addDefects()
{
	std::vector<int> step_defects;
	for (auto c : c_map)
	{
		if (std::fabs(c.second.getm()-0.5) < 1e-3 || std::fabs(c.second.getm()+0.5) < 1e-3) step_defects.push_back(c.first);
	}
	defects.push_back(step_defects);
}

const std::vector<int>& Global::stepDefects(int step) const { return defects[step]; }

std::unordered_map<int, Vertex>& Global::vertexMap() { return v_map; }
std::unordered_map<int, Edge>& Global::edgeMap() { return e_map; }
std::unordered_map<int, Cell>& Global::cellMap() { return c_map; }

Vertex& Global::vert(int v) { return v_map.at(v); }
Edge& Global::edge(int e) { return e_map.at(e); }
Cell& Global::cell(int c) { return c_map.at(c); }

const int Global::vertexCounter() const { return v_c; }
const int Global::edgeCounter() const { return e_c; }
const int Global::cellCounter() const { return c_c; }


const int Global::createVertex(Point r)
{
	v_map.emplace(v_c, Vertex(this, r));
	return v_c++;
}
const int Global::createEdge(int v1, int v2)
{
	e_map.emplace(e_c, Edge(this, v1, v2));
	v_map.at(v1).addEdgeContact(e_c); //vertex v1 knows it's part of edge
	v_map.at(v2).addEdgeContact(e_c); //vertex v2 knows it's part of edge
	return e_c++; //return id of the created edge
}
const int Global::createCell(std::vector<int>& vertices, std::vector<std::pair<int,int>>& edges)
{
	c_map.emplace(c_c, Cell(this, vertices, edges));
	for (int v : vertices) { v_map.at(v).addCellContact(c_c); } //vertices know they are part of cell
	for (std::pair<int,int> e : edges) { e_map.at(e.first).addCellJunction(c_c); } //edges know they are part of cell
	return c_c++; //return id of created cell
}


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

void Global::cellRemoveEdge(int c, int e)
{
	c_map.at(c).removeEdge(e);
}


void Global::destroyVertex(int v) { v_map.erase(v);}
void Global::destroyEdge(int e) 
{ 
	for (int c : e_map.at(e).cellJunctions()) { c_map.at(c).removeEdge(e); }
	v_map.at(e_map.at(e).E().first).removeEdgeContact(e);
	v_map.at(e_map.at(e).E().second).removeEdgeContact(e);
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
		if (cell.second.getA() < 0.2*A_0)
		{
			const std::vector<int>& vertices = c_map.at(cell.first).Vertices();
			int i = 0; bool contact = false;
			while (i < vertices.size() && contact == false)
			{
				for (int c : v_map.at(vertices[i]).cellContacts())
				{
					auto it = std::find(small_cells.begin(), small_cells.end(), c);
					if (it != small_cells.end()) { contact = true; break; }
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
		if (cell.second.getA() > 2*A_0)
		{
			const std::vector<int>& vertices = c_map.at(cell.first).Vertices();
			int i = 0; bool contact = false;
			while (i < vertices.size() && contact == false)
			{
				for (int c : v_map.at(vertices[i]).cellContacts())
				{
					auto it = std::find(large_cells.begin(), large_cells.end(), c);
					if (it != large_cells.end()) { contact = true; break; }
				}
				i++;
			}
			if (!contact) large_cells.push_back(cell.first);
		}
	}
	for (int c : large_cells) { c_map.at(c).divide(); }
}

void Global::run()
{
	while (timestep < timesteps)
	{
	}
}


