#include "edge.h"


Edge::Edge(Global* g, int v1, int v2) : g(g), id(g->e_c()), v_1(v1), v_2(v2) {}

bool Edge::operator==(const Edge& other) const { return ((v_1 == other.v_1) && (v_2 == other.v_2)) || ((v_1== other.v_2) && (v_2 == other.v_1)); }

const int Edge::v1() const { return v_1; }
const int Edge::v2() const { return v_2; }
const std::unordered_set<int>& Edge::cellJunctions() const { return cell_junctions; }
const double Edge::l() const { return l_; }
const double Edge::T_l() const { return T_l_; }


void Edge::addCellJunction(int cell_id) { cell_junctions.insert(cell_id); }

void Edge::removeCellJunction(int cell_id) 
{
    cell_junctions.erase(cell_id);
    if (cell_junctions.size() == 0) 
    { 
		g->vert(v_1).removeEdgeContact(id);
		g->vert(v_2).removeEdgeContact(id);
		g->edgeMap().erase(id); 
	}
}

const bool Edge::hasVertex(int v) const { return (v == v_1 || v == v_2); }

bool Edge::swapVertex_rep(int v_old, int v_new)
{
	if (v_old == v_1) 
	{
		v_1 = v_new;
		for (int c : cell_junctions) g->cellExchangeVertex(c, v_old, v_new);
		g->vert(v_new).addEdgeContact(id);
		g->vert(v_old).removeEdgeContact(id);
		return true;
	}
	else if (v_old == v_2) 
	{
		v_2 = v_new;
		for (int c : cell_junctions) g->cellExchangeVertex(c, v_old, v_new);
		g->vert(v_new).addEdgeContact(id);
		g->vert(v_old).removeEdgeContact(id);
		return true;
	}
	else { return false;}
}


void Edge::calcLength() { l_ = std::sqrt((g->vert(v_1).r()-g->vert(v_2).r()).squared_length()); }

void Edge::calcT_l()
{
	T_l_ = T_l_0;
	for (int c: cell_junctions) T_l_ += k_L*g->cell(c).L();
}


void Edge::T1merge()
{
	Point p = CGAL::midpoint(g->vert(v_1).r(), g->vert(v_2).r());
	const int v = g->createVertex(p);

	
	
}




