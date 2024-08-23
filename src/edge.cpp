#include "edge.h"


Edge::Edge(Tissue* T, int v1, int v2) : T(T), id(T->e_c()), v_1(v1), v_2(v2) {}

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
    if (cell_junctions.size() == 0) T->destroyEdge(id);
}

const bool Edge::hasVertex(int v) const { return (v == v_1 || v == v_2); }

bool Edge::swapVertex(int v_old, int v_new)
{
	if (v_old == v_1) 
	{
		v_1 = v_new;
		for (int c : cell_junctions) T->cellExchangeVertex(c, v_old, v_new);
		T->vert(v_new).addEdgeContact(id);
		T->vert(v_old).removeEdgeContact(id);
		return true;
	}
	else if (v_old == v_2) 
	{
		v_2 = v_new;
		for (int c : cell_junctions) T->cellExchangeVertex(c, v_old, v_new);
		T->vert(v_new).addEdgeContact(id);
		T->vert(v_old).removeEdgeContact(id);
		return true;
	}
	else { return false;}
}

bool Edge ::swapVertex_noedit(int v_old, int v_new)
{
	if (v_old == v_1) 
	{
		v_1 = v_new;
		T->vert(v_new).addEdgeContact(id);
		T->vert(v_old).removeEdgeContact(id);
		return true;
	}
	else if (v_old == v_2) 
	{
		v_2 = v_new;
		T->vert(v_new).addEdgeContact(id);
		T->vert(v_old).removeEdgeContact(id);
		return true;
	}
	else { return false;}
}


void Edge::calcLength() { l_ = std::sqrt((T->vert(v_1).r()-T->vert(v_2).r()).squared_length()); }

void Edge::calcT_l()
{
	T_l_ = param::LAMBDA;
	for (int c: cell_junctions) T_l_ += param::GAMMA*T->cell(c).L();
}


void Edge::T1merge()
{
	if (T->vert(v_1).edgeContacts().size() != 3 || T->vert(v_2).edgeContacts().size() != 3) { return; }
	if (cell_junctions.size() < 2) { return; }
	
	/*std::cout << "Edge: " << id << '\t' << v_1 << ' ' << v_2 << '\n';
	for (int c : T->vert(v_1).cellContacts()) T->cell(c).outputVertices();
	for (int c : T->vert(v_2).cellContacts()) T->cell(c).outputVertices();
	std::cout << '\n';*/
	Point p = CGAL::midpoint(T->vert(v_1).r(), T->vert(v_2).r());
	const int v_new = T->createVertex(p);
	
	auto side = [this, v_new](int v)
	{
		std::unordered_set<int>::const_iterator e_it = T->vert(v).edgeContacts().begin();
		std::array<int, 2> corner; int i = 0;
		while (i < 2)
		{
			if (*e_it != id) { corner[i] = *e_it; i++; }
			e_it = std::next(e_it);
		}
 		T->edge(corner[0]).swapVertex_noedit(v, v_new);
		T->edge(corner[1]).swapVertex_noedit(v, v_new);		
		
		for (int c : T->vert(v).cellContacts()) if (T->cell(c).hasEdge(corner[0]) && T->cell(c).hasEdge(corner[1])) T->cellExchangeVertex(c, v, v_new);
		
	};
	side(v_1); side(v_2);
	
	for (int c : cell_junctions) { T->cellExchangeVertex(c, v_1, v_new); T->cellRemoveVertex(c, v_2); }
	std::unordered_set<int>::const_iterator c_a_it = cell_junctions.begin();
	std::unordered_set<int>::const_iterator c_b_it = std::next(c_a_it);
	T->cell(*c_a_it).removeEdge(id); T->cell(*c_b_it).removeEdge(id);
	
	//for (int c : cell_junctions) { T->cell(c).outputEdgeVertices(); T->cell(c).outputVertices(); }
	//for (int c : T->vert(v_1).cellContacts()) { T->cell(c).outputEdgeVertices(); T->cell(c).outputVertices(); }
	//for (int c : T->vert(v_2).cellContacts()) { T->cell(c).outputEdgeVertices(); T->cell(c).outputVertices(); }

}

