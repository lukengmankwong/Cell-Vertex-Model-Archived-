#ifndef VERTEX_H
#define VERTEX_H

#include <set>

#include "libraries.h"
#include "parameters.h"

class Vertex
{
private:
    std::unordered_map<int, Vertex>* vertex_map;

    int id;
    Point r;
    std::set<int> cell_contacts;

    Vec force;

public:
    Vertex(std::unordered_map<int, Vertex>* vertex_map, int id, Point r);
    bool operator==(const Vertex& other) const;

    void addCellContact(int cell_id);
    void removeCellContact(int cell_id);
    const std::set<int>& getCellContacts() const;

    const Point& getR() const;

    void calcForce(Point centroid);
    void applyForce();

};

#endif // VERTEX_H
