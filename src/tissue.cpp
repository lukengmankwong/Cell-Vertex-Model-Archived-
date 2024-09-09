#include "tissue.h"


Tissue::Tissue() : v_c_(0), e_c_(0), c_c_(0), timestep(0)
{
	cell_defects.reserve(1000);
	vertex_defects.reserve(1000);
}
Tissue::~Tissue() {}

const std::unordered_map<int, Vertex>& Tissue::vertexMap() const { return v_map; }
const std::unordered_map<int, Edge>& Tissue::edgeMap() const { return e_map; }
const std::unordered_map<int, Cell>& Tissue::cellMap() const { return c_map; }
const int Tissue::v_c() const { return v_c_; }
const int Tissue::e_c() const { return e_c_; }
const int Tissue::c_c() const { return c_c_; }

const std::vector<int>& Tissue::cellStepDefects() const { return cell_defects; }
const std::vector<int>& Tissue::vertexStepDefects() const { return vertex_defects; }

Vertex& Tissue::vert(int v) { return v_map.at(v); }
Edge& Tissue::edge(int e) { return e_map.at(e); }
Cell& Tissue::cell(int c) { return c_map.at(c); }

void Tissue::addDefects()
{
	std::vector<int> cell_step_defects;
	for (const auto& cell : c_map) { if (std::fabs(cell.second.m()-0.5) < 1e-3 || std::fabs(cell.second.m()+0.5) < 1e-3
		|| std::fabs(cell.second.m()-1) < 1e-3 || std::fabs(cell.second.m()+1) < 1e-3) cell_step_defects.push_back(cell.first); }
	cell_defects = cell_step_defects;
	
	std::vector<int> vertex_step_defects;
	for (const auto& vertex : v_map) { if (std::fabs(vertex.second.m()-0.5) < 1e-3 || std::fabs(vertex.second.m()+0.5) < 1e-3
		|| std::fabs(vertex.second.m()-1) < 1e-3 || std::fabs(vertex.second.m()+1) < 1e-3) vertex_step_defects.push_back(vertex.first); }
	vertex_defects = vertex_step_defects;
}

const int Tissue::createVertex(Point r)
{
	v_map.emplace(v_c_, Vertex(this, r));
	return v_c_++; 																	//return id of created vertex
}
const int Tissue::createEdge(const int v1, const int v2)
{
	e_map.emplace(e_c_, Edge(this, v1, v2));
	v_map.at(v1).addEdgeContact(e_c_); 												//vertex v1 knows it's part of edge
	v_map.at(v2).addEdgeContact(e_c_); 												//vertex v2 knows it's part of edge
	return e_c_++; 																	//return id of created edge
}
const int Tissue::createCell(std::vector<int>& vertices, std::vector<int>& edges)
{
	c_map.emplace(c_c_, Cell(this, vertices, edges));
	for (int v : vertices) { v_map.at(v).addCellContact(c_c_); } 					//vertices know they are part of cell
	for (int e : edges) { e_map.at(e).addCellJunction(c_c_); } 	//edges know they are part of cell
	return c_c_++; 																	//return id of created cell
}

void Tissue::destroyVertex(int v) { v_map.erase(v);}
void Tissue::destroyEdge(int e) 
{ 
	for (int c : e_map.at(e).cellJunctions()) { c_map.at(c).removeEdge(e); }
	v_map.at(e_map.at(e).v1()).removeEdgeContact(e);
	v_map.at(e_map.at(e).v2()).removeEdgeContact(e);
	e_map.erase(e); 
}
void Tissue::destroyCell(int c)
{
	for (int v : c_map.at(c).Vertices()) { v_map.at(v).removeCellContact(c); }
	for (int e : c_map.at(c).Edges()) { e_map.at(e).removeCellJunction(c); }
	c_map.erase(c);
}


void Tissue::cellNewVertex(int c, int v, int i) 
{ 
	c_map.at(c).addVertex(v, i);  
	v_map.at(v).addCellContact(c);
}
void Tissue::cellRemoveVertex(int c, int v) { c_map.at(c).removeVertex(v); }
void Tissue::cellExchangeVertex(int c, int v_old, int v_new) { c_map.at(c).exchangeVertex(v_old, v_new); }

void Tissue::cellNewEdge(int c, int e, int i)
{
	c_map.at(c).addEdge(e, i);
	e_map.at(e).addCellJunction(c);
}
int Tissue::cellRemoveEdge(int c, int e) { return c_map.at(c).removeEdge(e); }


void Tissue::cellsFindNeighbours() { for (auto& cell : c_map) cell.second.findNeighbours(); }
void Tissue::verticesFindNeighbours() { for (auto& vertex : v_map) vertex.second.orderCellContacts(); }


const bool Tissue::commonEdge(int c1, int c2) const
{
	for (int e : c_map.at(c1).Edges())
	{
		auto it = std::find(c_map.at(c2).Edges().begin(), c_map.at(c2).Edges().end(), e);
		if (it != c_map.at(c2).Edges().end()) { return true; }
	}
	return false;
}

const double Tissue::D_angle(int c_i, int c_j) const
{
	double Z_i = c_map.at(c_i).Z(); double X_i = c_map.at(c_i).X(); double S_i = std::sqrt(X_i*X_i+Z_i*Z_i);
	double Z_j = c_map.at(c_j).Z(); double X_j = c_map.at(c_j).X(); double S_j = std::sqrt(X_j*X_j+Z_j*Z_j);
	return std::asin( (Z_i*X_j-X_i*Z_j) / (S_i*S_j) );
}


void Tissue::extrusion()
{
	std::vector<int> small_cells;
	for (const auto& cell : c_map)
	{
		if (cell.second.A() < param::A_min)
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

void Tissue::division()
{
	std::vector<int> large_cells;
	for (const auto& cell : c_map)
	{
		if (cell.second.A() > param::A_max)
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

void Tissue::transitions()
{
	for (std::pair<const int,Cell>& cell : c_map) cell.second.calcA();
	extrusion();
	for (std::pair<const int,Cell>& cell : c_map) cell.second.calcA();
	division();
}

void Tissue::T1()
{
	std::vector<int> short_edges;
	for (const auto& edge : e_map) if ( edge.second.l() < param::l_min ) { short_edges.push_back(edge.first); }
	for (int e : short_edges) { e_map.at(e).T1(); break; }
	std::vector<int> fourfold_vertices;
	for (const auto& vertex : v_map) if (vertex.second.edgeContacts().size() == 4) fourfold_vertices.push_back(vertex.first);
	for (int v : fourfold_vertices){ v_map.at(v).T1split(); }
}

void Tissue::run(int max_timestep)
{
	for (auto& vertex : v_map) vertex.second.onBoundaryCell();
	while (timestep < max_timestep)
	{
		transitions();
		
		for (auto& edge : e_map) 	{ edge.second.calcLength(); }
		for (auto& cell : c_map) 
		{ 
			cell.second.calcL();
			cell.second.calcA(); 
			cell.second.calcT_A();
			cell.second.calcG();
		}
		
		for (auto& edge : e_map)	{ edge.second.calcT_l(); }		
		for (auto& vertex : v_map) 	{ vertex.second.calcForce(); }		
        for (auto& vertex : v_map) 	{ vertex.second.applyForce(); }
        
		if (timestep % 20 == 0)
		{
			for (auto& cell : c_map) 	{ cell.second.calcm(); }
			for (auto& vertex : v_map) 	{ vertex.second.calcm(); }
			addDefects();
			writeCellsFile(this, "cells" + std::to_string(timestep) + ".vtk");
			writeDirectorsFile(this, "directors" + std::to_string(timestep) + ".vtk");
			writeCellDefectsFile(this, "cell defects" + std::to_string(timestep) + ".vtk");
			writeVertexDefectsFile(this, "vertex defects" + std::to_string(timestep) + ".vtk");
		}
		T1();
        timestep++;	
	}
}
