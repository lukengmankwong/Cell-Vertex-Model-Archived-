#ifndef VERTEX_H
#define VERTEX_H

#include <cmath>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <array>
#include <iostream>

#include "parameters.h"
#include "libraries.h"

class Tissue;

class Edge;
class Cell;


class Vertex
{
private:

	Tissue* T;
    int id;
    Point r_;
    Vec force_;
    double m_;
    int not_boundary_cell;
    
    std::unordered_set<int> cell_contacts;
    std::unordered_set<int> edge_contacts;
    
    std::unordered_set<Edge*> edge_contacts_;
    std::unordered_set<Cell*> cell_contacts_;
    
	std::vector<std::pair<int, double>> cell_contacts_ordered;
    
    Vec calcSurfaceForce();
    Vec calcLineForce();

public:
    
    Vertex(Tissue* T, Point r);
    Vertex();
    bool operator==(const Vertex& other) const;
    bool onBoundaryCell();
    
    const Point& r() const;
    const double m() const;
    const std::unordered_set<int>& cellContacts() const;
    const std::unordered_set<int>& edgeContacts() const;

    void addCellContact(int cell_id);
    void removeCellContact(int cell_id);
    
    void addEdgeContact(int edge_id);
    void removeEdgeContact(int edge_id);

    void calcForce();
    void applyForce();
    void shearForce();

	void orderCellContacts();
	void T1split();
	void calcm();
	
};

#endif // VERTEX_H
