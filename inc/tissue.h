#ifndef TISSUE_H
#define TISSUE_H

#include <unordered_map>

#include "libraries.h"
#include "vertex.h"
#include "edge.h"
#include "cell.h"
class Vertex;
class Edge;
class Cell;

#include "parameters.h"
#include "functions.h"

class Tissue
{
private:

	std::unordered_map<int, Vertex> v_map; int v_c_;
	std::unordered_map<int, Edge> e_map; int e_c_;
	std::unordered_map<int, Cell> c_map; int c_c_;
	
	int timestep;
	std::vector<int> cell_defects;
	std::vector<int> vertex_defects;
	
	void extrusion();
	void division();
	void transitions();
	void T1();
	void addDefects();

	
public:
	Tissue();
	~Tissue();
    
    const std::vector<int>& cellStepDefects() const;
    const std::vector<int>& vertexStepDefects() const;
    
    const std::unordered_map<int, Vertex>& vertexMap() const;
    const std::unordered_map<int, Edge>& edgeMap() const ;
    const std::unordered_map<int, Cell>& cellMap() const;
    
    Vertex& vert(int v);
    Edge& edge(int e);
    Cell& cell(int c);
    
    const int v_c() const; 
    const int e_c() const;
    const int c_c() const;

	const int createVertex(Point r);
	const int createEdge(int v1, int v2);
	const int createCell(std::vector<int>& vertices, std::vector<int>& edges);
	
	void cellNewVertex(int c, int v, int i); //add new vertex to cell vertices at index i
	void cellNewEdge(int c, int e, int i); //add new edge to cell edges at index i
	
	void cellRemoveVertex(int c, int v);
	int cellRemoveEdge(int c, int e);
	
	void cellExchangeVertex(int c, int v_old, int v_new);
	
	void destroyVertex(int v);
	void destroyEdge(int e);
	void destroyCell(int c);
	
	const bool commonEdge(int c1, int c2) const;
	const double D_angle(int c_i, int c_j) const; 
	
	void cellFindNeighbours();
	void run(int max_timestep);
	
};

#endif // TISSUE_H
