#include "edge.h"

Edge::Edge(int v1, int v2)
{
    (v1 < v2) ? e = std::make_pair(v1, v2) : e = std::make_pair(v2, v1);
}

bool Edge::operator==(const Edge& other) const { return ((e.first == other.e.first) && (e.second == other.e.second)); }


std::pair<int, int> Edge::getE() const { return e; }


void Edge::calcLength(std::vector<Vertex>& vertices) {
    length = std::sqrt((vertices[e.first].getR()-vertices[e.second].getR()).squared_length());
}


