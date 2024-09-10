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

	std::array<Vertex, V_ARR_SIZE> v_arr; 	int v_c_;
	std::array<Edge, E_ARR_SIZE> e_arr; 	int e_c_;
	std::array<Cell, C_ARR_SIZE> c_arr;		int c_c_;
	std::array<bool, V_ARR_SIZE> v_in;
	std::array<bool, E_ARR_SIZE> e_in;
	std::array<bool, C_ARR_SIZE> c_in;
	
	int timestep;
	std::vector<int> cell_defects;
	std::vector<int> vertex_defects;
	
	void extrusion();
	void division();
	void transitions();
	void T1();
	void addDefects();
	
public:

	Tissue(VD& vd, bool (*in)(const Point&));
	~Tissue();
	
	const bool v_alive(int v) const;
    
    std::vector<int> cells();
	std::vector<int> edges();
	std::vector<int> vertices();
    
    const int v_c() const; 
    const int e_c() const;
    const int c_c() const;
    
    const std::vector<int>& cellStepDefects() const;
    const std::vector<int>& vertexStepDefects() const;
    
    Vertex& vert(int v);
    Edge& edge(int e);
    Cell& cell(int c);
    
	const int createVertex(Point r);
	const int createEdge(int v1, int v2);
	const int createCell(std::vector<int>& vertices, std::vector<int>& edges);
	
	void destroyVertex(int v);
	void destroyEdge(int e);
	void destroyCell(int c);
		
	void cellNewVertex(int c, int v, int i); 	//add new vertex to cell vertices at index i
	void cellRemoveVertex(int c, int v);
	void cellExchangeVertex(int c, int v_old, int v_new);
	
	void cellNewEdge(int c, int e, int i); 		//add new edge to cell edges at index i
	int cellRemoveEdge(int c, int e);
	
	void cellsFindNeighbours();
	void verticesFindNeighbours();
	
	const double D_angle(int c_i, int c_j) const; 
	
	void run(int max_timestep);
	
};

#endif // TISSUE_H
