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
    Point r_;
    Vec force_;
    double m_;
    int not_boundary_cell;
    
    std::unordered_set<Edge*> edge_contacts_;
    std::unordered_set<Cell*> cell_contacts_;
    
	std::vector<std::pair<Cell*, double>> cell_contacts_ordered;
    
    Vec calcSurfaceForce();
    Vec calcLineForce();

public:
    
    Vertex(Tissue* T, Point r);
    Vertex();
    bool operator==(const Vertex& other) const;
    bool onBoundaryCell();
    
    const Point& r() const;
    const double m() const;
    const std::unordered_set<Cell*>& cellContacts() const;
    const std::unordered_set<Edge*>& edgeContacts() const;

    void addCellContact(Cell* c);
    void removeCellContact(Cell* c);
    
    void addEdgeContact(Edge* e);
    void removeEdgeContact(Edge* e);

    void calcForce();
    void applyForce();
    void shearForce();

	void orderCellContacts();
	void T1split();
	void calcm();
	
};

#endif // VERTEX_H
