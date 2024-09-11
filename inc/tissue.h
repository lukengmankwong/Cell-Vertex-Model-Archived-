#ifndef TISSUE_H
#define TISSUE_H

#include <unordered_map>
#include <array>

#include "libraries.h"
#include "vertex.h"
#include "edge.h"
#include "cell.h"

#include "parameters.h"
#include "functions.h"

#define V_ARR_SIZE 6000
#define E_ARR_SIZE 9000
#define C_ARR_SIZE 3000

class Tissue
{
private:

	std::array<Vertex, V_ARR_SIZE> v_arr;
	std::array<Edge, E_ARR_SIZE> e_arr;
	std::array<Cell, C_ARR_SIZE> c_arr;
	std::array<bool, V_ARR_SIZE> v_in;
	std::array<bool, E_ARR_SIZE> e_in;
	std::array<bool, C_ARR_SIZE> c_in;
	
	Vertex* v_0_; Vertex* v_c_;
	Edge* e_0_; Edge* e_c_;
	Cell* c_0_; Cell* c_c_;
	
	std::vector<Cell*> c_def_PLUSHALF_;
	std::vector<Cell*> c_def_PLUSONE_;
	std::vector<Cell*> c_def_MINUSHALF_;
	std::vector<Cell*> c_def_MINUSONE_;
	std::vector<Vertex*> v_def_PLUSHALF_;
	std::vector<Vertex*> v_def_PLUSONE_;
	std::vector<Vertex*> v_def_MINUSHALF_;
	std::vector<Vertex*> v_def_MINUSONE_;
	
	int timestep;
	
	void extrusion();
	void division();
	void transitions();
	void T1();
	void findDefects();
	
public:

	Tissue(VD& vd, bool (*in)(const Point&));
	~Tissue();
	
	const bool v_alive(Vertex* v) const;
    
    std::vector<Vertex*> vertices();
	std::vector<Edge*> edges();
    std::vector<Cell*> cells();	
    
    Vertex* const v_c() const; Vertex* const v_0() const;
    Edge* const e_c() const; Edge* const e_0() const;
    Cell* const c_c() const; Cell* const c_0() const;

	const std::vector<Cell*>& c_def_PLUSHALF() const;
	const std::vector<Cell*>& c_def_PLUSONE() const;
	const std::vector<Cell*>& c_def_MINUSHALF() const;
	const std::vector<Cell*>& c_def_MINUSONE() const;
	const std::vector<Vertex*>& v_def_PLUSHALF() const;
	const std::vector<Vertex*>& v_def_PLUSONE() const;
	const std::vector<Vertex*>& v_def_MINUSHALF() const;
	const std::vector<Vertex*>& v_def_MINUSONE() const;
    
	Vertex* const createVertex(Point r);
	Edge* const createEdge(Vertex* v_1, Vertex* v_2);
	Cell* const createCell(std::vector<Vertex*>& vertices, std::vector<Edge*>& edges);
	
	void destroyVertex(Vertex* v);
	void destroyEdge(Edge* e);
	void destroyCell(Cell* c);
		
	void cellNewVertex(Cell* c, Vertex* v, int i); 	//add new vertex to cell vertices at index i
	void cellRemoveVertex(Cell* c, Vertex* v);
	void cellExchangeVertex(Cell* c, Vertex* v_old, Vertex* v_new);
	
	void cellNewEdge(Cell* c, Edge* e, int i); 		//add new edge to cell edges at index i
	int cellRemoveEdge(Cell* c, Edge* e);
	
	const double D_angle(Cell* c_i, Cell* c_j) const; 
	
	void run(int max_timestep);
	
};

#endif // TISSUE_H
