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

    Point r_0_;
    double A_; double S_; //cell area, area sign (+1 or -1)
    double L_; //cell perimeter
    double T_A_; //surface tension
    
    double G[3]; //gyration tensor symmetric so only need 3 values (a b, b c)
    double lambda; //gyration tensor largest eigenvalue
    double Z_, X_; //for defect analysis
    Vec n_; //normalised director
    double m_; //winding number around cell nearest neighbors
    
    const int longestEdge_i() const;
    std::vector<int> nearestNeighbours();
    
public:
    Cell(Global* g, std::vector<int>& vertices, std::vector<int>& edges);
    
    void outputVertices() const;
    void outputEdges() const;
    void outputEdgeVertices() const;
    
    const std::vector<int>& Vertices() const;
    const std::vector<int>& Edges() const;
    const Point& r_0() const;
    const double A() const; 
    const double S() const;
    const double L() const;
    const double T_A() const;
    const Vec& n() const;
    const double Z() const; 
    const double X() const;
    const double m() const;
    
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
    void calcA();
    void calcL();
    void calcT_A();
    void calcG();
    void calcm();
    
    bool valid();

};

#endif // CELL_H
