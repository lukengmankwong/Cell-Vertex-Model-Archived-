#ifndef VERTEX_H
#define VERTEX_H

#include <unordered_set>
#include <vector>

#include "tissue.h"
class Tissue;

#include "parameters.h"

class Vertex
{
private:

	Tissue* T;
    const int id;
    Point r_;
    Vec force_;
    double m_;
    
    std::unordered_set<int> cell_contacts;
    std::unordered_set<int> edge_contacts;
    
    Vec calcSurfaceForce();
    Vec calcLineForce();
    
    int not_boundary_cell;
    
	void orderCellContacts();
    std::vector<std::pair<int, double>> cell_contacts_ordered;

public:
    
    Vertex(Tissue* T, Point r);
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

	void T1split();
	
	void calcm();
	
};

#endif // VERTEX_H
