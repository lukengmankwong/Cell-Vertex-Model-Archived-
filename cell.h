#ifndef CELL_H
#define CELL_H

#include <vector>

#include <CGAL/Origin.h>

#include "vertex.h"

class Cell
{
private:
    std::vector<int> vertices; //indices to cell vertices in vertex list
    std::vector<int> edges; //indices to cell edges in edge list

    Point centroid;

public: 

    Cell(std::vector<int>& vertices, std::vector<int>& edges);
    void calcCentroid(std::vector<Vertex>& all_vertices);

    Point getCentroid() const;
    std::vector<int>& getVertices();
    
};
#endif // CELL_H
