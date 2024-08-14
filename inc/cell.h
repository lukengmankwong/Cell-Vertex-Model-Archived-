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
    std::vector<int> vertices;
    std::vector<int> edges;
    std::vector<int> nearest_neighbours;
    std::vector<bool> edge_directions;

    Point centroid;
    double A; double S; //cell area, area sign (+1 or -1)
    double L; //cell perimeter
    double T_A; //surface tension
    
    double G[3]; //gyration tensor symmetric so only need 3 values (a b, b c)
    double TL[3]; //traceless matrix of gyration tensor
    double lambda; //gyration tensor largest eigenvalue
    Vec director; //normalised
    double m; //winding number around cell nearest neighbors
    
    const int longestEdge_i() const;
    std::vector<int> nearestNeighbours();
    
public:
    Cell(Global* g, std::vector<int>& vertices, std::vector<std::pair<int,int>>& edges);
    
    const std::vector<int>& Vertices() const;
    const std::vector<int>& Edges() const;
    
    const double getA() const;
    const double getS() const;
    const double getL() const;
    const double getT_A() const;
    const Point& getCentroid() const;
    const Vec& getDirector() const;
    const double getm() const;
    
    const bool hasEdge(int e) const;
    
    void addVertex(int v, int i);
    void removeVertex(int v);
    void addEdge(int e, int i);
    void removeEdge(int e);
 
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
