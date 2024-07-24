#include "edge.h"


Edge::Edge(int id, int v1, int v2) : id(id)
{
    (v1 < v2) ? e = std::make_pair(v1, v2) : e = std::make_pair(v2, v1);
    cell_junction_count = 0;
}

bool Edge::operator==(const Edge& other) const { return ((e.first == other.e.first) && (e.second == other.e.second)); }

const int Edge::getID() const { return id; }
const std::pair<int, int>& Edge::getE() const { return e; }
int Edge::getCellJunctions() const { return cell_junction_count; }

void Edge::addCellJunction() { cell_junction_count++; }

void Edge::removeCellJunction() {
    cell_junction_count--;
    if (cell_junction_count == 0) {
        edge_map.erase(id);
    }
}


bool Edge::swapVertex(int v_old, int v_new)
{
	if (v_old == e.first) {
		e.first = v_new;
		return true;
	}
	else if (v_old == e.second) {
		e.second = v_new;
		return true;
	}
	else { return false;}
}


void Edge::calcLength()
{

}
