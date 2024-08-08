#ifndef CELL_H
#define CELL_H

#include <cmath>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <array>

#include <CGAL/Origin.h>
#include "global.h"
class Global;


class Cell
{
private:

	Global* g;
	const int id;
    std::vector<int> vertex_keys;
    std::vector<int> edge_keys;

    Point centroid;
    double A; double S; //cell area, area sign (+1 or -1)
    double L; //cell perimeter
    double T_A; //surface tension
    
    double G[3]; //gyration tensor symmetric so only need 3 values (a b, b c)
    double TL[3]; //traceless matrix of gyration tensor
    double lambda; //gyration tensor largest eigenvalue
    Vec director;
    double m; //winding number around cell nearest neighbors
    
    const int longestEdge() const;
    
public:
    Cell(Global* g, int id, std::vector<int>& vertex_keys, std::vector<int>& edge_keys);
     
    const int getID() const;
    const Point& getCentroid() const;
    const Vec& getDirector() const;
    
    const double getA() const;
    const double getS() const;
    const double getL() const;
    const double getT_A() const;
    
    const std::vector<int>& getVertices() const;
    const std::vector<int>& getEdges() const;
    
    
    void addVertex(int v, int i);
    void addEdge(int e, int i);
 
    void removeEdges();
    void removeVertices();
    
    void extrude();
	void divide();
    
    void calcCentroid();
    void calcG();
    void calcA();
    void calcL();
    void calcT_A();
    void calcm();

};

#endif // CELL_H
