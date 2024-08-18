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

    Point r_0;
    double A; double S; //cell area, area sign (+1 or -1)
    double L; //cell perimeter
    double T_A; //surface tension
    
    double G[3]; //gyration tensor symmetric so only need 3 values (a b, b c)
    double lambda; //gyration tensor largest eigenvalue
    double Z, X; //for defect analysis
    Vec n; //normalised director
    double m; //winding number around cell nearest neighbors
    
    const int longestEdge_i() const;
    std::vector<int> nearestNeighbours();
    
public:
    Cell(Global* g, std::vector<int>& vertices, std::vector<std::pair<int,int>>& edges);
    
    void outputVertices() const;
    void outputEdges() const;
    void outputEdgeVertices() const;
    
    const std::vector<int>& Vertices() const;
    const std::vector<int>& Edges() const;
    
    const double getA() const;
    const double getS() const;
    const double getL() const;
    const double getT_A() const;
    const Point& R_0() const;
    const Vec& N() const;
    const double getZ() const; const double getX() const;
    const double getm() const;
    
    const bool hasEdge(int e) const;
    const bool onBoundary() const;
    
    void addVertex(int v, int i);
    void removeVertex(int v);
    void addEdge(int e, int i);
    int removeEdge(int e);
 
    void exchangeVertex(int v_old, int v_new);
    void rotateVertices();
    
    void extrude();
	void divide();
    
    void calcR_0();
    void calcG();
    void calcA();
    void calcL();
    void calcT_A();
    void calcm();
    
    bool valid();

};

#endif // CELL_H
