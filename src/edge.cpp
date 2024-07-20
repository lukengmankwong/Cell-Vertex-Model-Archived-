#include "edge.h"

Edge::Edge(std::unordered_map<int, Edge>* edge_map, int id, int v1, int v2) :
    edge_map(edge_map), id(id)
{
    (v1 < v2) ? e = std::make_pair(v1, v2) : e = std::make_pair(v2, v1);
    cell_junction_count = 0;
}

bool Edge::operator==(const Edge& other) const { return ((e.first == other.e.first) && (e.second == other.e.second)); }


void Edge::addCellJunction() { cell_junction_count++; }
void Edge::removeCellJunction() {
    cell_junction_count--;
    if (cell_junction_count == 0) {
        (*edge_map).erase(id);
    }
}
int Edge::getCellJunctions() const { return cell_junction_count; }


const std::pair<int, int>& Edge::getE() const { return e; }


/*void Edge::calcLength(std::vector<Vertex>& vertices) {
    length = std::sqrt((vertices[e.first].getR()-vertices[e.second].getR()).squared_length());
}*/


