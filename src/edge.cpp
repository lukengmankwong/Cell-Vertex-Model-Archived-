#include "edge.h"


Edge::Edge(Global* g, int v1, int v2) : g(g), id(g->edgeCounter())
{
	e = std::make_pair(v1, v2);
}

bool Edge::operator==(const Edge& other) const { return ((e.first == other.e.first) && (e.second == other.e.second)) || ((e.first == other.e.second) && (e.second == other.e.first)); }

const std::pair<int, int>& Edge::E() const { return e; }
const std::unordered_set<int>& Edge::cellJunctions() const { return cell_junctions; }
const double Edge::getl() const { return l; }
const double Edge::getT_l() const { return T_l; }


void Edge::addCellJunction(int cell_id) { cell_junctions.insert(cell_id); }

void Edge::removeCellJunction(int cell_id) 
{
    cell_junctions.erase(cell_id);
    if (cell_junctions.size() == 0) 
    { 
		g->vert(e.first).removeEdgeContact(id);
		g->vert(e.second).removeEdgeContact(id);
		g->edgeMap().erase(id); 
	}
}

const bool Edge::hasVertex(int v) const { return (v == e.first || v == e.second); }

bool Edge::swapVertex(int v_old, int v_new)
{
	if (v_old == e.first) {
		e.first = v_new;
		for (int c : cell_junctions) { g->cell(c).removeVertex(v_old); }
		g->vert(v_old).removeEdgeContact(id);
		return true;
	}
	else if (v_old == e.second) {
		e.second = v_new;
		for (int c : cell_junctions) { g->cell(c).removeVertex(v_old); }
		g->vert(v_old).removeEdgeContact(id);
		return true;
	}
	else { return false;}
}


void Edge::calcLength() { l = std::sqrt((g->vert(e.first).R()-g->vert(e.second).R()).squared_length()); }

void Edge::calcT_l()
{
	T_l = T_l_0;
	for (int c: cell_junctions) {
		T_l += k_L*g->cell(c).getL();
	}
}

