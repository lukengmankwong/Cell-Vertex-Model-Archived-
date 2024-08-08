#include "edge.h"


Edge::Edge(Global* g, int id, int v1, int v2) : g(g), id(id)
{
    (v1 < v2) ? e = std::make_pair(v1, v2) : e = std::make_pair(v2, v1);
}

bool Edge::operator==(const Edge& other) const { return ((e.first == other.e.first) && (e.second == other.e.second)); }

const int Edge::getID() const { return id; }
const std::pair<int, int>& Edge::getE() const { return e; }
const std::unordered_set<int>& Edge::getCellJunctions() const { return cell_junctions; }
const double Edge::getl() const { return l; }
const double Edge::getT_l() const { return T_l; }


void Edge::addCellJunction(int cell_id) { cell_junctions.insert(cell_id); }

void Edge::removeCellJunction(int cell_id) 
{
    cell_junctions.erase(cell_id);
    if (cell_junctions.size() == 0) { 
		g->vertexMap().at(e.first).removeEdgeContact(id);
		g->vertexMap().at(e.second).removeEdgeContact(id);
		g->edgeMap().erase(id); 
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


void Edge::calcLength() { l = std::sqrt((g->vertexMap().at(e.first).getR()-g->vertexMap().at(e.second).getR()).squared_length()); }

void Edge::calcT_l()
{
	T_l = T_l_0;
	for (int c: cell_junctions) {
		T_l += k_L*g->cellMap().at(c).getL();
	}
}

