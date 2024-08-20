#ifndef VERTEX_H
#define VERTEX_H

#include <utility>
#include <unordered_set>

#include "global.h"
class Global;

#include "parameters.h"

class Vertex
{
private:

	Global* g;
    const int id;
    Point r_;
    Vec force_;
    double m_;
    
    std::unordered_set<int> cell_contacts;
    std::unordered_set<int> edge_contacts;
    
    Vec calcSurfaceForce();
    Vec calcLineForce();
    
    std::vector<int> orderCellContacts();

public:
    Vertex(Global* g, Point r);
    bool operator==(const Vertex& other) const;
    
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

	void T1split();
	
	void calcm();
	
};

#endif // VERTEX_H
