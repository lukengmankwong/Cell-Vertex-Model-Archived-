#ifndef EDGE_H
#define EDGE_H

#include <vector>
#include <cmath>
#include <utility>

#include "vertex.h"

class Edge
{
private:
    std::pair<int, int> e;
    double length;

public:
    Edge(int v1, int v2);
    bool operator==(const Edge& other) const;

    std::pair<int, int> getE() const;
    
    void calcLength(std::vector<Vertex>& vertices);

};

#endif // EDGE_H