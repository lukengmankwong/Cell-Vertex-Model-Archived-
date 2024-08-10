#ifndef GLOBALS_H
#define GLOBALS_H

#include <utility>
#include <unordered_map>

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel                   K;
typedef CGAL::Vector_2<K>                                                   Vec;
typedef CGAL::Point_2<K>												  Point;

#include "vertex.h"
#include "edge.h"
#include "cell.h"

class Vertex;
class Edge;
class Cell;

//constants

extern const double pi;

//parameters
extern const int cell_count;

extern const double dt;
extern const int timesteps;

extern const double A_0; //characteristic cell area

extern const double k_A; //spring constant for surface tension
extern const double k_L; //spring constant for edge tension
extern const double T_l_0; //independent edge tension

extern const double a;

class Global //singleton
{
private:
	
	Global();
	~Global();

	std::unordered_map<int, Vertex> vertex_map; int vertex_counter;
	std::unordered_map<int, Edge> edge_map; int edge_counter;
	std::unordered_map<int, Cell> cell_map; int cell_counter;
	
public:

	Global(const Global&) = delete;
    Global& operator=(const Global&) = delete;
    static Global& get() 
    {
        static Global instance;
        return instance;
    }
    
    std::unordered_map<int, Vertex>& vertexMap();
    std::unordered_map<int, Edge>& edgeMap();
    std::unordered_map<int, Cell>& cellMap();
    
    const int vertexCounter() const; 
    const int edgeCounter() const;
    const int cellCounter() const;

	const int createVertex(Point r);
	const int createEdge(int v1, int v2);
	const int createCell(std::vector<int>& vertex_keys, std::vector<int>& edge_keys);
	
	void cellNewVertex(int c, int v, int i); //add new vertex to cell vertex_keys at index i
	void cellNewEdge(int c, int e, int i=-1); //add new edge to cell edge_keys at index i
	
	void destroyVertex(int v);
	void destroyEdge(int e);
	void destroyCell(int c);
	
	const bool commonEdge(int c1, int c2) const;
	
	void extrusion();
	
};

#endif // GLOBALS_H
