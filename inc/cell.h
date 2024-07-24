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
	
    std::unordered_set<int> vertex_keys;
    std::unordered_set<int> edge_keys;

    Point centroid;
    double A; double dA;
    double T; //surface tension
    
    double G[3]; //gyration tensor symmetric so only need 3 values (a b, b c)
    double lambda; //gyration tensor largest eigenvalue
    Vec director;
    
public:

    Cell(int id, std::vector<int>& vertex_keys, std::vector<int>& edge_keys);
     
    const int getID() const;
    const Point& getCentroid() const;
    double getA() const;
    double getdA() const;
    double getT() const;
    const std::unordered_set<int>& getVertices() const;
    const std::unordered_set<int>& getEdges() const;
   
	bool removeEdge(int edge_id);
	bool removeVertex(int vertex_id);
        
    void removeEdges();
    void removeVertices();
    
    void extrude();
    
    void calcCentroid();
    void calcG();
    void calcArea();
    void calcT();

};

#endif // CELL_H
