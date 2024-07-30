#ifndef CELL_H
#define CELL_H

#include <cmath>
#include <vector>
#include <unordered_set>

#include <CGAL/Origin.h>
#include "globals.h"


class Cell
{
private:

	const int id;
	
    std::vector<int> vertex_keys;
    std::unordered_set<int> edge_keys;

    Point centroid;
    double A; double dA;
    double L; //cell perimeter
    double T_A; //surface tension
    
    double G[3]; //gyration tensor symmetric so only need 3 values (a b, b c)
    double lambda; //gyration tensor largest eigenvalue
    Vec director;
    
public:

    Cell(int id, std::vector<int>& vertex_keys, std::vector<int>& edge_keys);
     
    const int getID() const;
    const Point& getCentroid() const;
    const Vec& getDirector() const;
    
    const double getA() const;
    const double getdA() const;
    const double getL() const;
    const double getT_A() const;
    
    const std::vector<int>& getVertices() const;
    const std::unordered_set<int>& getEdges() const;
   
	bool removeEdge(int edge_id);
	bool removeVertex(int vertex_id);
        
    void removeEdges();
    void removeVertices();
    
    void extrude();
    
    void calcCentroid();
    void calcG();
    void calcA();
    void calcL();
    void calcT_A();

};

#endif // CELL_H
