#include "tissue.h"

void removeDuplicates(std::vector<Edge*>& vec) {
    std::unordered_set<Edge*> seen;   // To track seen elements
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

Tissue::Tissue(VD& vd, bool (*in)(const Point&)) : v_0_(v_arr.data()), e_0_(e_arr.data()), c_0_(c_arr.data()), timestep(0)
{
	v_in = {false};
	e_in = {false};
	c_in = {false};
	v_c_ = v_0_;
	e_c_ = e_0_;
	c_c_ = c_0_;
	
	def_PLUSHALF_c = {0};
	def_PLUSONE_c = {0};
	def_MINUSHALF_c = {0};
	def_MINUSONE_c = {0};
	
	std::cout << "COLLECTING INITIAL DATA\n";
	for (VD::Vertex_iterator vit = vd.vertices_begin(); vit != vd.vertices_end(); vit++) createVertex(vit->point());
    
    for (VD::Face_iterator fi = vd.faces_begin(); fi != vd.faces_end(); fi++) 
    {
        std::vector<Vertex*> cell_vertices;
        VD::Ccb_halfedge_circulator ec_start = fi->ccb();
        VD::Ccb_halfedge_circulator ec = ec_start;
        
        do { 
			if (!ec->is_unbounded()) 
			{
				for (Vertex* v = v_0_; v < v_c_; v++)
				{
					if (v_in[v-v_0_])
					{
						if (ec->source()->point() == v->r())
						{
							cell_vertices.push_back(v);
							break;
						}
					}
				}
			}
			++ec;
        } while (ec != ec_start);
		
		std::vector<Edge*> cell_edges; 
		size_t n = cell_vertices.size();
		for (int i = 0; i < n; i++)
		{
			Vertex* v_1 = cell_vertices[i]; 
			Vertex* v_2 = cell_vertices[(i+1)%n];
			bool found = false;
			for (Edge* e = e_0_; e < e_c_; e++)
			{
				if (e_in[e-e_0_])
				{
					if ( (e->v1() == v_1 && e->v2() == v_2) || (e->v1() == v_2 && e->v2() == v_1) )
					{
						cell_edges.push_back(e);
						found = true;
						break;
					}
				}
			}
			if (!found) cell_edges.push_back(createEdge(v_1, v_2));
		}
		
		removeDuplicates(cell_edges); //temporary fix
        createCell(cell_vertices, cell_edges);
    } 
    
	std::unordered_set<Cell*> cells_to_remove;
	for (Vertex* v = v_0_; v < v_c_; v++)
	{
		if (!in(v->r())) cells_to_remove.insert(v->cellContacts().begin(), v->cellContacts().end());
	}
	for (Cell* c = c_0_; c < c_c_; c++)
	{
		if (c->edges().size() < 3) cells_to_remove.insert(c);
	}
	for (Cell* c : cells_to_remove) destroyCell(c);
	    
    for (Cell* c = c_0_; c < c_c_; c++) if (c_in[c-c_0_]) c->findNeighbours(); 			//cells find neighbours
	for (Vertex* v = v_0_; v < v_c_; v++) if (v_in[v-v_0_]) v->orderCellContacts();		//vertices order cell contacts
	
	//sanity check using Euler characteristic: we expect Euler = 1
	int V = vertices().size(); int E = edges().size(); int C = cells().size();
	int Euler = V-E+C;
    std::cout << "V=" << V << "\nE=" << E << "\nC=" << C << "\nV-E+C=" << Euler << '\n';
}
Tissue::~Tissue() {}

const bool Tissue::v_alive(Vertex* v) const { return v_in[v-v_0_]; }

std::vector<Vertex*> Tissue::vertices()
{
	std::vector<Vertex*> vertices;
	for (Vertex* v = v_0_; v < v_c_; v++) if (v_in[v-v_0_]) vertices.push_back(v);
	return vertices;
}
std::vector<Edge*> Tissue::edges()
{
	std::vector<Edge*> edges;
	for (Edge* e = e_0_; e < e_c_; e++) if (e_in[e-e_0_]) edges.push_back(e);
	return edges;
}
std::vector<Cell*> Tissue::cells()
{
	std::vector<Cell*> cells;
	for (Cell* c = c_0_; c < c_c_; c++) if (c_in[c-c_0_]) cells.push_back(c);
	return cells;
}

Vertex* const Tissue::v_0() const { return v_0_; }
Edge* const Tissue::e_0() const { return e_0_; }
Cell* const Tissue::c_0() const { return c_0_; }
Vertex* const Tissue::v_c() const { return v_c_; }
Edge* const Tissue::e_c() const { return e_c_; }
Cell* const Tissue::c_c() const { return c_c_; }

const std::vector<Cell*>& Tissue::c_def_PLUSHALF() const { return c_def_PLUSHALF_; }
const std::vector<Cell*>& Tissue::c_def_PLUSONE() const { return c_def_PLUSONE_; }
const std::vector<Cell*>& Tissue::c_def_MINUSHALF() const { return c_def_MINUSHALF_; }
const std::vector<Cell*>& Tissue::c_def_MINUSONE() const { return c_def_MINUSONE_; }
const std::vector<Vertex*>& Tissue::v_def_PLUSHALF() const { return v_def_PLUSHALF_; }
const std::vector<Vertex*>& Tissue::v_def_PLUSONE() const { return v_def_PLUSONE_; }
const std::vector<Vertex*>& Tissue::v_def_MINUSHALF() const { return v_def_MINUSHALF_; }
const std::vector<Vertex*>& Tissue::v_def_MINUSONE() const { return v_def_MINUSONE_; }

void Tissue::findDefects()
{
	c_def_PLUSHALF_ = {};
	c_def_PLUSONE_ = {};
	c_def_MINUSHALF_ = {};
	c_def_MINUSONE_ = {};
	for (Cell* c = c_0_; c < c_c_; c++)
	{
		if (c_in[c-c_0_])
		{
			double m = c->m();
			if 		(std::fabs(m - 0.5) < 1e-3) c_def_PLUSHALF_.push_back(c); 
			else if (std::fabs(m - 1) < 1e-3) c_def_PLUSONE_.push_back(c); 
			else if (std::fabs(m + 0.5) < 1e-3) c_def_MINUSHALF_.push_back(c); 
			else if (std::fabs(m + 1) < 1e-3) c_def_MINUSONE_.push_back(c); 
		}
	}
	
	v_def_PLUSHALF_ = {};
	v_def_PLUSONE_ = {};
	v_def_MINUSHALF_ = {};
	v_def_MINUSONE_ = {};
	for (Vertex* v = v_0_; v < v_c_; v++)
	{
		if (v_in[v-v_0_])
		{
			double m = v->m();
			if 		(std::fabs(m - 0.5) < 1e-3) v_def_PLUSHALF_.push_back(v); 
			else if (std::fabs(m - 1) < 1e-3) v_def_PLUSONE_.push_back(v); 
			else if (std::fabs(m + 0.5) < 1e-3) v_def_MINUSHALF_.push_back(v); 
			else if (std::fabs(m + 1) < 1e-3) v_def_MINUSONE_.push_back(v); 
		}
	}
}

void Tissue::countDefects()
{
	for (Cell* c = c_0_; c < c_c_; c++)
	{
		if (c_in[c-c_0_])
		{
			double m = c->m();
			if 		(std::fabs(m - 0.5) < 1e-3) def_PLUSHALF_c[timestep]++; 
			else if (std::fabs(m - 1) < 1e-3) def_PLUSONE_c[timestep]++; 
			else if (std::fabs(m + 0.5) < 1e-3) def_MINUSHALF_c[timestep]++; 
			else if (std::fabs(m + 1) < 1e-3) def_MINUSONE_c[timestep]++; 
		}
	}
	for (Vertex* v = v_0_; v < v_c_; v++)
	{
		if (v_in[v-v_0_])
		{
			double m = v->m();
			if 		(std::fabs(m - 0.5) < 1e-3) def_PLUSHALF_c[timestep]++; 
			else if (std::fabs(m - 1) < 1e-3) def_PLUSONE_c[timestep]++; 
			else if (std::fabs(m + 0.5) < 1e-3) def_MINUSHALF_c[timestep]++; 
			else if (std::fabs(m + 1) < 1e-3) def_MINUSONE_c[timestep]++; 
		}
	}
}

Vertex* const Tissue::createVertex(Point r)
{
	*v_c_ = Vertex(this, r); v_in[v_c_-v_0_] = true;
	return v_c_++; 																	//return id of created vertex
}
Edge* const Tissue::createEdge(Vertex* v1, Vertex* v2)
{	
	v1->addEdgeContact(e_c_);														//vertex v1 knows it's part of edge
	v2->addEdgeContact(e_c_);														//vertex v1 knows it's part of edge
	*e_c_ = Edge(this, v1, v2); e_in[e_c_-e_arr.data()] = true;
	return e_c_++; 																	//return id of created edge
}
Cell* const Tissue::createCell(std::vector<Vertex*>& vertices, std::vector<Edge*>& edges)
{
	for (Vertex* v : vertices) v->addCellContact(c_c_);								//vertices know they are part of cell
	for (Edge* e : edges) e->addCellJunction(c_c_);									//edges know they are part of cell	
	*c_c_ = Cell(this, vertices, edges); c_in[c_c_-c_0_] = true;	
	return c_c_++; 																	//return id of created cell
}

void Tissue::destroyVertex(Vertex* v) { v_in[v-v_0_] = false; }
void Tissue::destroyEdge(Edge* e) 
{ 
	for (Cell* c : e->cellJunctions()) c->removeEdge(e);
	e->v1()->removeEdgeContact(e);
	e->v2()->removeEdgeContact(e);
	e_in[e-e_0_] = false;
}
void Tissue::destroyCell(Cell* c)
{
	for (Vertex* v : c->vertices()) v->removeCellContact(c);
	for (Edge* e : c->edges()) e->removeCellJunction(c);
	c_in[c-c_0_] = false;
}


void Tissue::cellNewVertex(Cell* c, Vertex* v, int i) 
{ 
	c->addVertex(v, i);
	v->addCellContact(c);
}
void Tissue::cellRemoveVertex(Cell* c, Vertex* v) { c->removeVertex(v); }
void Tissue::cellExchangeVertex(Cell* c, Vertex* v_old, Vertex* v_new) { c->exchangeVertex(v_old, v_new); }

void Tissue::cellNewEdge(Cell* c, Edge* e, int i)
{
	c->addEdge(e, i);
	e->addCellJunction(c);	
}
int Tissue::cellRemoveEdge(Cell* c, Edge* e) { return c->removeEdge(e); }

const double Tissue::D_angle(Cell* c_i, Cell* c_j) const
{	
	double Z_i = c_i->Z(); double X_i = c_i->X(); double S_i = std::sqrt(X_i*X_i+Z_i*Z_i);
	double Z_j = c_j->Z(); double X_j = c_j->X(); double S_j = std::sqrt(X_j*X_j+Z_j*Z_j);
	return std::asin( (Z_i*X_j-X_i*Z_j) / (S_i*S_j) );
}


void Tissue::extrusion()
{	
	std::vector<Cell*> small_cells;
	for (Cell* c = c_0_; c < c_c_; c++)
	{
		if (c_in[c-c_0_])
		{
			if (c->A() < param::A_min)
			{
				const std::vector<Vertex*>& vertices = c->vertices();
				int j = 0; bool contact = false;
				while (j < vertices.size() && !contact)
				{
					for (Cell* v_cell : vertices[j]->cellContacts())
					{
						std::vector<Cell*>::const_iterator it = std::find(small_cells.begin(), small_cells.end(), v_cell);
						if (it != small_cells.end()) contact = true;
					}
					j++;
				}
				if (!contact) small_cells.push_back(c);
			}
		}
	}
	for (Cell* c : small_cells) c->extrude();
}

void Tissue::division()
{	
	std::vector<Cell*> large_cells;
	for (Cell* c = c_0_; c < c_c_; c++)
	{
		if (c_in[c-c_0_])
		{
			if (c->A() > param::A_max)
			{
				const std::vector<Vertex*>& vertices = c->vertices();
				int j = 0; bool contact = false;
				while (j < vertices.size() && !contact)
				{
					for (Cell* v_cell : vertices[j]->cellContacts())
					{
						std::vector<Cell*>::const_iterator it = std::find(large_cells.begin(), large_cells.end(), v_cell);
						if (it != large_cells.end()) contact = true;
					}
					j++;
				}
				if (!contact) large_cells.push_back(c);
			}
		}
	}
	for (Cell* c : large_cells) c->divide();
}

void Tissue::transitions()
{	
	for (int c = 0; c < c_c_-c_0_; c++) if (c_in[c]) c_arr[c].calcA();
	extrusion();
	for (int c = 0; c < c_c_-c_0_; c++) if (c_in[c]) c_arr[c].calcA();
	division();
}

void Tissue::T1()
{
	std::vector<Edge*> short_edges;
	for (Edge* e = e_0_; e < e_c_; e++)
	{
		if (e_in[e-e_0_])
		{
			if (e->l() < param::l_min)
			{
				const std::unordered_set<Cell*>& cells = e->cellJunctions();
				bool contact = false;
				for (Cell* c : cells)
				{
					for (Edge* c_edge : c->edges())
					{
						std::vector<Edge*>::const_iterator it = std::find(short_edges.begin(), short_edges.end(), c_edge);
						if (it != short_edges.end()) { contact = true; break; }
					}
				}
				if (!contact) short_edges.push_back(e);
			}
		}
	}
	for (Edge* e : short_edges) { e->T1(); }
	std::vector<int> fourfold_vertices;
	for (int v = 0; v < v_c_-v_0_; v++) if (v_in[v]) if (v_arr[v].edgeContacts().size() == 4) fourfold_vertices.push_back(v);
	for (int v : fourfold_vertices) v_arr[v].T1split();
}

void Tissue::run(int max_timestep, std::string title)
{
	for (int v = 0; v < v_c_-v_0_; v++) if (v_in[v]) v_arr[v].onBoundaryCell();
	while (timestep < max_timestep)
	{
		if (timestep % 1000 == 0) std::cout << timestep << '\n';
		transitions();
		
        for (int e = 0; e < e_c_-e_0_; e++) if (e_in[e]) e_arr[e].calcLength();
        for (int c = 0; c < c_c_-c_0_; c++)
        {
			if (c_in[c])
			{
				c_arr[c].calcL();
				c_arr[c].calcA();
				c_arr[c].calcT_A();
				c_arr[c].calcG();
			}
		}
		
		for (int e = 0; e < e_c_-e_0_; e++) if (e_in[e]) e_arr[e].calcT_l();
		for (int v = 0; v < v_c_-v_0_; v++) if (v_in[v]) v_arr[v].calcForce();
		for (int v = 0; v < v_c_-v_0_; v++) if (v_in[v]) v_arr[v].applyForce();
		
		for (int c = 0; c < c_c_-c_0_; c++) if (c_in[c]) c_arr[c].calcm();
		for (int v = 0; v < v_c_-v_0_; v++) if (v_in[v]) v_arr[v].calcm();
		countDefects();
        
		/*if (timestep % 20 == 0)
		{
			findDefects();
			writeCellsFile(this, "cells" + std::to_string(timestep) + ".vtk");
			writeDirectorsFile(this, "directors" + std::to_string(timestep) + ".vtk");
			
			writeCellDefectsFile(this, c_def_PLUSHALF_, "cell defects PLUSHALF" + std::to_string(timestep) + ".vtk");
			writeCellDefectsFile(this, c_def_PLUSONE_, "cell defects PLUSONE" + std::to_string(timestep) + ".vtk");
			writeCellDefectsFile(this, c_def_MINUSHALF_, "cell defects MINUSHALF" + std::to_string(timestep) + ".vtk");
			writeCellDefectsFile(this, c_def_MINUSONE_, "cell defects MINUSONE" + std::to_string(timestep) + ".vtk");
			writeVertexDefectsFile(this, v_def_PLUSHALF_, "vertex defects PLUSHALF" + std::to_string(timestep) + ".vtk");
			writeVertexDefectsFile(this, v_def_PLUSONE_, "vertex defects PLUSONE" + std::to_string(timestep) + ".vtk");
			writeVertexDefectsFile(this, v_def_MINUSHALF_, "vertex defects MINUSHALF" + std::to_string(timestep) + ".vtk");
			writeVertexDefectsFile(this, v_def_MINUSONE_, "vertex defects MINUSONE" + std::to_string(timestep) + ".vtk");
		}*/
		T1();
        timestep++;	
	}
	std::ofstream plushalf(title + "PLUSHALF.txt");
	for (int i = 0; i < max_timestep; i++) plushalf << def_PLUSHALF_c[i] << "\n";
	plushalf.close();
	
	std::ofstream plusone(title + "PLUSONE.txt");
	for (int i = 0; i < max_timestep; i++) plusone << def_PLUSONE_c[i] << "\n";
	plusone.close();
	
	std::ofstream minushalf(title + "MINUSHALF.txt");
	for (int i = 0; i < max_timestep; i++) minushalf << def_MINUSHALF_c[i] << "\n";
	minushalf.close();
	
	std::ofstream minusone(title + "MINUSONE.txt");
	for (int i = 0; i < max_timestep; i++) minusone << def_MINUSONE_c[i] << "\n";
	minusone.close();
	
}
