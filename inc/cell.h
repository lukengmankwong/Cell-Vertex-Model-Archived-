#ifndef CELL_H
#define CELL_H

#include <unordered_map>
#include <vector>
#include <unordered_map>

#include <CGAL/Origin.h>

#include "vertex.h"
#include "edge.h"

class Cell
{
private:
    int id;

    std::unordered_map<int, Vertex>* vertex_map;
    std::unordered_map<int, Edge>* edge_map;

    std::vector<int> vertex_indices; //indices to cell vertices in vertex list
    std::vector<int> edge_indices; //indices to cell edges in edge list

    Point centroid;
    double area;

public:

    Cell(int id, std::unordered_map<int, Vertex>* vertex_map, std::vector<int>& vertex_indices, std::unordered_map<int, Edge>* edge_map, std::vector<int>& edge_indices);
    ~Cell();
    void kill();
    void extrude();

    void calcCentroid();
    void calcArea();

    const int getID() const;
    const Point& getCentroid() const;
    double getArea() const;
    const std::vector<int>& getVertices() const;
    const std::vector<int>& getEdges() const;

};
#endif // CELL_H
