#include "tissue.h"

void removeDuplicates(std::vector<int>& vec) {
    std::unordered_set<int> seen;   // To track seen elements
    auto it = vec.begin();

    while (it != vec.end()) {
        // If the element is seen for the first time, keep it
        if (seen.insert(*it).second) {
            ++it;
        } else {
            // Otherwise, erase the duplicate
            it = vec.erase(it);
        }
    }
}

Tissue::Tissue(VD& vd, bool (*in)(const Point&)) : v_c_(0), e_c_(0), c_c_(0), timestep(0)
{
	v_in = {false};
	e_in = {false};
	c_in = {false};
	
	std::cout << "COLLECTING INITIAL DATA\n";
	for (VD::Vertex_iterator vit = vd.vertices_begin(); vit != vd.vertices_end(); vit++) createVertex(vit->point());
    
    for (VD::Face_iterator fi = vd.faces_begin(); fi != vd.faces_end(); fi++) 
    {
        std::vector<int> cell_vertices;
        VD::Ccb_halfedge_circulator ec_start = fi->ccb();
        VD::Ccb_halfedge_circulator ec = ec_start;
        
        do { 
			if (!ec->is_unbounded()) 
			{
				for (int v = 0; v < v_c_; v++)
				{
					if (v_in[v])
					{
						if (ec->source()->point() == v_arr[v].r())
						{
							cell_vertices.push_back(v);
							break;
						}
					}
				}
			}
			++ec;
        } while (ec != ec_start);
		
		std::vector<int> cell_edges; 
		for (int i = 0; i < cell_vertices.size(); i++)
		{
			int v1 = cell_vertices[i]; int v2 = cell_vertices[(i+1)%cell_vertices.size()];
			bool found = false;
			for (int e = 0; e < e_c_; e++)
			{
				if (e_in[e])
				{
					if ( (e_arr[e].v1() == v1 && e_arr[e].v2() == v2) || (e_arr[e].v1() == v2 && e_arr[e].v2() == v1) )
					{
						cell_edges.push_back(e);
						found = true;
						break;
					}
				}
			}
			if (!found) cell_edges.push_back(createEdge(v1, v2));
		}
		
		removeDuplicates(cell_edges); //temporary fix
        createCell(cell_vertices, cell_edges);
    } 
    
	std::unordered_set<int> cells_to_remove;
	for (int v = 0; v < v_c_; v++)
	{
		if (!in(v_arr[v].r())) cells_to_remove.insert(v_arr[v].cellContacts().begin(), v_arr[v].cellContacts().end());
	}
	for (int c = 0; c < c_c_; c++)
	{
		if (c_arr[c].Edges().size() < 3) cells_to_remove.insert(c);
	}
	for (int c : cells_to_remove) destroyCell(c);
	
	int V = vertices().size(); int E = edges().size(); int C = cells().size();
	int Euler = V-E+C;
    std::cout << "V=" << V << "\nE=" << E << "\nC=" << C << "\nV-E+C=" << Euler << '\n';	
	cellsFindNeighbours();
	verticesFindNeighbours();
}
Tissue::~Tissue() {}

const bool Tissue::v_alive(int v) const { return v_in[v]; }


std::vector<int> Tissue::vertices()
{
	std::vector<int> vertices;
	for (int v = 0; v < v_c_; v++) if (v_in[v]) vertices.push_back(v);
	return vertices;
}
std::vector<int> Tissue::edges()
{
	std::vector<int> edges;
	for (int e = 0; e < e_c_; e++) if (e_in[e]) edges.push_back(e);
	return edges;
}
std::vector<int> Tissue::cells()
{
	std::vector<int> cells;
	for (int c = 0; c < c_c_; c++) if (c_in[c]) cells.push_back(c);
	return cells;
}


const int Tissue::v_c() const { return v_c_; }
const int Tissue::e_c() const { return e_c_; }
const int Tissue::c_c() const { return c_c_; }

const std::vector<int>& Tissue::cellStepDefects() const { return cell_defects; }
const std::vector<int>& Tissue::vertexStepDefects() const { return vertex_defects; }

Vertex& Tissue::vert(int v) 
{ 
	if (v_in[v]) return v_arr[v]; 
}
Edge& Tissue::edge(int e) 
{ 
	if (e_in[e]) return e_arr[e]; 
}
Cell& Tissue::cell(int c) 
{ 
	if (c_in[c]) return c_arr[c]; 
}

void Tissue::addDefects()
{
	std::vector<int> cell_step_defects;
	for (int c = 0; c < c_c_; c++) {
		if (c_in[c]) {
			if (std::fabs(c_arr[c].m()-0.5) < 1e-3 || std::fabs(c_arr[c].m()+0.5) < 1e-3
				|| std::fabs(c_arr[c].m()-1) < 1e-3 || std::fabs(c_arr[c].m()+1) < 1e-3) cell_step_defects.push_back(c);
		}
	} cell_defects = cell_step_defects;
	
	std::vector<int> vertex_step_defects;
	for (int v = 0; v < v_c_; v++) {
		if (v_in[v]) {
			if (std::fabs(v_arr[v].m()-0.5) < 1e-3 || std::fabs(v_arr[v].m()+0.5) < 1e-3
				|| std::fabs(v_arr[v].m()-1) < 1e-3 || std::fabs(v_arr[v].m()+1) < 1e-3) vertex_step_defects.push_back(v);
		} 
	} vertex_defects = vertex_step_defects;
}

const int Tissue::createVertex(Point r)
{
	v_arr[v_c_] = Vertex(this, r); v_in[v_c_] = true;
	return v_c_++; 																	//return id of created vertex
}
const int Tissue::createEdge(const int v1, const int v2)
{	
	e_arr[e_c_] = Edge(this, v1, v2); e_in[e_c_] = true;
	v_arr[v1].addEdgeContact(e_c_);													//vertex v1 knows it's part of edge
	v_arr[v2].addEdgeContact(e_c_);													//vertex v1 knows it's part of edge
	return e_c_++; 																	//return id of created edge
}
const int Tissue::createCell(std::vector<int>& vertices, std::vector<int>& edges)
{
	for (int v : vertices) v_arr[v].addCellContact(c_c_);							//vertices know they are part of cell
	for (int e : edges) e_arr[e].addCellJunction(c_c_);								//edges know they are part of cell	
	c_arr[c_c_] = Cell(this, vertices, edges); c_in[c_c_] = true;	
	return c_c_++; 																	//return id of created cell
}

void Tissue::destroyVertex(int v) { v_in[v] = false; }
void Tissue::destroyEdge(int e) 
{ 
	for (int c : e_arr[e].cellJunctions()) c_arr[c].removeEdge(e);
	v_arr[e_arr[e].v1()].removeEdgeContact(e);
	v_arr[e_arr[e].v2()].removeEdgeContact(e);
	e_in[e] = false;
	
}
void Tissue::destroyCell(int c)
{
	for (int v : c_arr[c].Vertices()) v_arr[v].removeCellContact(c);
	for (int e : c_arr[c].Edges()) e_arr[e].removeCellJunction(c);
	c_in[c] = false;
}


void Tissue::cellNewVertex(int c, int v, int i) 
{ 
	c_arr[c].addVertex(v, i);
	v_arr[v].addCellContact(c);	
}
void Tissue::cellRemoveVertex(int c, int v) { c_arr[c].removeVertex(v); }
void Tissue::cellExchangeVertex(int c, int v_old, int v_new) { c_arr[c].exchangeVertex(v_old, v_new); }

void Tissue::cellNewEdge(int c, int e, int i)
{
	c_arr[c].addEdge(e, i);
	e_arr[e].addCellJunction(c);	
}
int Tissue::cellRemoveEdge(int c, int e) { return c_arr[c].removeEdge(e); }


void Tissue::cellsFindNeighbours() { for (int c = 0; c < c_c_; c++) if (c_in[c]) c_arr[c].findNeighbours(); }
void Tissue::verticesFindNeighbours() { for (int v = 0; v < v_c_; v++) if (v_in[v]) v_arr[v].orderCellContacts(); }

const double Tissue::D_angle(int c_i, int c_j) const
{	
	double Z_i = c_arr[c_i].Z(); double X_i = c_arr[c_i].X(); double S_i = std::sqrt(X_i*X_i+Z_i*Z_i);
	double Z_j = c_arr[c_j].Z(); double X_j = c_arr[c_j].X(); double S_j = std::sqrt(X_j*X_j+Z_j*Z_j);
	return std::asin( (Z_i*X_j-X_i*Z_j) / (S_i*S_j) );
}


void Tissue::extrusion()
{	
	std::vector<int> small_cells;
	for (int i = 0; i < c_c_; i++)
	{
		if (c_in[i])
		{
			if (c_arr[i].A() < param::A_min)
			{
				const std::vector<int>& vertices = c_arr[i].Vertices();
				int j = 0; bool contact = false;
				while (j < vertices.size() && !contact)
				{
					for (int c : v_arr[vertices[j]].cellContacts())
					{
						std::vector<int>::const_iterator it = std::find(small_cells.begin(), small_cells.end(), c);
						if (it != small_cells.end()) contact = true;
					}
					j++;
				}
				if (!contact) small_cells.push_back(i);
			}
		}
	}
	for (int c : small_cells) c_arr[c].extrude();
}

void Tissue::division()
{	
	std::vector<int> large_cells;
	for (int i = 0; i < c_c_; i++)
	{
		if (c_in[i])
		{
			if (c_arr[i].A() > param::A_max)
			{
				const std::vector<int>& vertices = c_arr[i].Vertices();
				int j = 0; bool contact = false;
				while (j < vertices.size() && !contact)
				{
					for (int c : v_arr[vertices[j]].cellContacts())
					{
						std::vector<int>::const_iterator it = std::find(large_cells.begin(), large_cells.end(), c);
						if (it != large_cells.end()) contact = true;
					}
					j++;
				}
				if (!contact) large_cells.push_back(i);
			}
		}
	}
	for (int c : large_cells) c_arr[c].divide();
}

void Tissue::transitions()
{	
	for (int i = 0; i < c_c_; i++) if (c_in[i]) c_arr[i].calcA();
	extrusion();
	for (int i = 0; i < c_c_; i++) if (c_in[i]) c_arr[i].calcA();
	division();

}

void Tissue::T1()
{
	std::vector<int> short_edges;
	for (int e = 0; e < e_c_; e++) if (e_in[e]) if (e_arr[e].l() < param::l_min) short_edges.push_back(e);
	for (int e : short_edges) { e_arr[e].T1(); break; }
	std::vector<int> fourfold_vertices;
	for (int v = 0; v < v_c_; v++) if (v_in[v]) if (v_arr[v].edgeContacts().size() == 4) fourfold_vertices.push_back(v);
	for (int v : fourfold_vertices) v_arr[v].T1split();
}

void Tissue::run(int max_timestep)
{
	for (int v = 0; v < v_c_; v++) if (v_in[v]) v_arr[v].onBoundaryCell();
	while (timestep < max_timestep)
	{
		transitions();
		
        for (int e = 0; e < e_c_; e++) if (e_in[e]) e_arr[e].calcLength();
        for (int c = 0; c < c_c_; c++)
        {
			if (c_in[c])
			{
				c_arr[c].calcL();
				c_arr[c].calcA();
				c_arr[c].calcT_A();
				c_arr[c].calcG();
			}
		}
		
		for (int e = 0; e < e_c_; e++) if (e_in[e]) e_arr[e].calcT_l();
		for (int v = 0; v < v_c_; v++) if (v_in[v]) v_arr[v].calcForce();
		for (int v = 0; v < v_c_; v++) if (v_in[v]) v_arr[v].applyForce();
        
		if (timestep % 10 == 0)
		{
			for (int c = 0; c < c_c_; c++) if (c_in[c]) c_arr[c].calcm();
			for (int v = 0; v < v_c_; v++) if (v_in[v]) v_arr[v].calcm();
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
