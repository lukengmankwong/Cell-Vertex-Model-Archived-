#ifndef GLOBALS_H
#define GLOBALS_H

#include <utility>
#include <unordered_map>

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

extern std::unordered_map<int, Vertex> vertex_map; extern int vertex_counter;
extern std::unordered_map<int, Edge> edge_map; extern int edge_counter;
extern std::unordered_map<int, Cell> cell_map; extern int cell_counter;


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
    
    const int vertexCounter() const; 

	void createVertex(Point r);
	void createEdge(int v1, int v2);
	void createCell(std::vector<int>& vertex_keys, std::vector<int>& edge_keys);
	
	void destroyVertex(int v);
	void destroyEdge(int e);
	void destroyCell(int c);
	
};

#endif // GLOBALS_H
