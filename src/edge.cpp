#include "edge.h"
#include <iostream>


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
        edge_map->erase(id);
    }
}
int Edge::getCellJunctions() const { return cell_junction_count; }


const int Edge::getID() const { return id; }
const std::pair<int, int>& Edge::getE() const { return e; }

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

