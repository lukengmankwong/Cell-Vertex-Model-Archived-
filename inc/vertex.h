#ifndef VERTEX_H
#define VERTEX_H

#include <unordered_set>

#include "libraries.h"
#include "parameters.h"
#include "edge.h"

class Vertex
{
private:
    std::unordered_map<int, Vertex>* vertex_map_ptr; int* vertex_counter_ptr;
    std::unordered_map<int, Edge>* edge_map_ptr; 

    const int id;
    Point r;
    Vec force;
    
    std::unordered_set<int> cell_contacts;
    std::unordered_set<int> edge_contacts;

public:
    Vertex(int id, int* vertex_counter_ptr, std::unordered_map<int, Vertex>* vertex_map_ptr, std::unordered_map<int, Edge>* edge_map_ptr, Point r);
    bool operator==(const Vertex& other) const;

    void addCellContact(int cell_id);
    void removeCellContact(int cell_id);
    const std::unordered_set<int>& getCellContacts() const;
    
    void addEdgeContact(int edge_id);
    void removeEdgeContact(int edge_id);

    const Point& getR() const;
    const int getID() const;
    
    bool inEdge(int edge_index) const;

    void calcForce(Point centroid);
    void applyForce();

};

#endif // VERTEX_H
