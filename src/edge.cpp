#include "edge.h"
#include "tissue.h"


Edge::Edge(Tissue* T, Vertex* v_1, Vertex* v_2) : 
	T(T)
{
	this->v_1 = v_1;
	this->v_2 = v_2;
}
Edge::Edge() = default;

bool Edge::operator==(const Edge& other) const { return ((v_1 == other.v_1) && (v_2 == other.v_2)) || ((v_1== other.v_2) && (v_2 == other.v_1)); }

Vertex* const Edge::v1()		const { return v_1; }
Vertex* const Edge::v2() 		const { return v_2; }
const double Edge::l() 		const { return l_; }
const double Edge::T_l()	const { return T_l_; }
const std::unordered_set<Cell*>& Edge::cellJunctions()	const { return cell_junctions_; }


void Edge::addCellJunction(Cell* c) { cell_junctions_.insert(c); }
void Edge::removeCellJunction(Cell* c) 
{
    cell_junctions_.erase(c);
    if (cell_junctions_.size() == 0) T->destroyEdge(this);
}

const bool Edge::hasVertex(Vertex* v) const { return (v == v_1 || v == v_2); }

bool Edge::swapVertex(Vertex* v_old, Vertex* v_new)
{
	if (v_old == v_1) 
	{
		v_1 = v_new;
		for (Cell* c : cell_junctions_) T->cellExchangeVertex(c, v_old, v_new);
		v_new->addEdgeContact(this);
		v_old->removeEdgeContact(this);
		return true;
	}
	else if (v_old == v_2) 
	{
		v_2 = v_new;
		for (Cell* c : cell_junctions_) T->cellExchangeVertex(c, v_old, v_new);
		v_new->addEdgeContact(this);
		v_old->removeEdgeContact(this);
		return true;
	}
	else { return false;}
}

bool Edge ::swapVertex_noedit(Vertex* v_old, Vertex* v_new)
{
	if (v_old == v_1) 
	{
		v_1 = v_new;
		v_new->addEdgeContact(this);
		v_old->removeEdgeContact(this);
		return true;
	}
	else if (v_old == v_2) 
	{
		v_2 = v_new;
		v_new->addEdgeContact(this);
		v_old->removeEdgeContact(this);
		return true;
	}
	else { return false;}
}


void Edge::calcLength() { l_ = std::sqrt((v_1->r()-v_2->r()).squared_length()); }

void Edge::calcT_l()
{
	T_l_ = param::LAMBDA;
	for (const Cell* c: cell_junctions_) T_l_ += param::GAMMA*c->L();
}


void Edge::T1() //problem with order of v_p and v_q;
{
	if (cell_junctions_.size() != 2) return;
	if (v_1->cellContacts().size() != 3 || v_2->cellContacts().size() != 3) return;
	
	//cells either side of edge
	Cell* const c_a = *(cell_junctions_.begin());
	Cell* const c_b = *std::next(cell_junctions_.begin());
	const std::unordered_set<Cell*> cellsAB = {c_a, c_b};
	
	//copy of edge contacts
	const std::unordered_set<Edge*> v_1_edges = v_1->edgeContacts();
	const std::unordered_set<Edge*> v_2_edges = v_2->edgeContacts();
	
	auto other_cell = [this, cellsAB](Vertex* const v) 
	{ 
		for (Cell* c : v->cellContacts()) if (cellsAB.find(c) == cellsAB.end()) return c; 
	};
	Cell* const c_p = other_cell(v_1);
	Cell* const c_q = other_cell(v_2);
	
	/*std::cout << "Cell P:\n"; T->cell(c_p).outputVertices(); T->cell(c_p).outputEdgeVertices();
	std::cout << "Cell Q:\n"; T->cell(c_q).outputVertices(); T->cell(c_q).outputEdgeVertices();*/
	
	//create new vertices
	Point cen = CGAL::midpoint(v_1->r(), v_2->r());
	Vec u = v_2->r() - v_1->r(); Vec s(-u.y(), u.x()); //s is u rotated 90 anticlockwise
	s *= (param::l_new/s.squared_length());
	Point a = cen + param::l_new*s; Point b = cen - param::l_new*s;
	c_a->calcR_0(); Point r_0 = c_a->r_0();
	if ( CGAL::squared_distance(a, r_0) > CGAL::squared_distance(b, r_0) ) std::swap(a,b);
	Vertex* const v_a = T->createVertex(a); Vertex* const v_b = T->createVertex(b);
	Edge* const e_new = T->createEdge(v_a, v_b);
	
	//replace edge in cells a and b with vertices a and b respectively
	auto edgeToVertex = [this](Cell* const c_x, Vertex* const v_x)
	{
		const std::vector<Edge*>& edges = c_x->edges();
		for (Edge* e : edges)
		{
			if (e != this) 
			{
				e->swapVertex(v_1, v_x);
				e->swapVertex(v_2, v_x);
			}
		}
		T->cellRemoveEdge(c_x, this);
	};
	edgeToVertex(c_a, v_a); edgeToVertex(c_b, v_b);
	/*std::cout << "Cell P:\n"; T->cell(c_p).outputVertices(); T->cell(c_p).outputEdgeVertices();
	std::cout << "Cell Q:\n"; T->cell(c_q).outputVertices(); T->cell(c_q).outputEdgeVertices();
	std::cout << "Cell A:\n"; T->cell(c_a).outputVertices(); T->cell(c_a).outputEdgeVertices();
	std::cout << "Cell B:\n"; T->cell(c_b).outputVertices(); T->cell(c_b).outputEdgeVertices();*/
	
	auto VertexToEdge = [this, c_a, c_b, v_a, v_b, e_new](Cell* const c_x, const std::unordered_set<Edge*>& v_edges)
	{
		const std::vector<Vertex*>& c_x_vertices = c_x->vertices();
		std::vector<Vertex*>::const_iterator it_v_a = std::find(c_x_vertices.begin(), c_x_vertices.end(), v_a);
		int i_v_a = std::distance(c_x_vertices.begin(), it_v_a);
		
		//std::cout << i_v_a << ' ' << vertices.size()-1 << '\n';
		
		Edge* e_a = nullptr; Edge* e_b = nullptr;
		for (Edge* e : v_edges)
		{
			if (c_a->hasEdge(e)) e_a = e;
			else if (c_b->hasEdge(e)) e_b = e;
		}
		
		const std::vector<Edge*>& c_x_edges = c_x->edges(); int n = c_x_edges.size();
		for (int i = 0; i < n; i++)
		{
			if (c_x_edges[i] == e_a && c_x_edges[(i+1)%n] == e_b)
			{
				T->cellNewEdge(c_x, e_new, (i+1)%n);
				T->cellNewVertex(c_x, v_b, (i+2)%n);
				break;
			}
			else if (c_x_edges[i] == e_b && c_x_edges[(i+1)%n] == e_a)
			{
				T->cellNewEdge(c_x, e_new, (i+1)%n);
				T->cellNewVertex(c_x, v_b, (i+1)%n);
				break;
			}
		}
	};
	VertexToEdge(c_p, v_1_edges); VertexToEdge(c_q, v_2_edges);

	/*std::cout << "Cell P:\n"; T->cell(c_p).outputVertices(); T->cell(c_p).outputEdgeVertices();
	std::cout << "Cell Q:\n"; T->cell(c_q).outputVertices(); T->cell(c_q).outputEdgeVertices();
	if (!(T->cell(c_p).valid())) { std::cout << "P invalid\n"; std::cin.get(); }
	if (!(T->cell(c_q).valid())) { std::cout << "Q invalid\n"; std::cin.get(); }*/
	while (!(c_p->valid())) { c_p->rotateVertices(); }
	while (!(c_q->valid())) { c_q->rotateVertices(); }
	/*if (!(T->cell(c_p).valid())) { std::cout << "P invalid\n"; std::cin.get(); }
	if (!(T->cell(c_q).valid())) { std::cout << "Q invalid\n"; std::cin.get(); }
	std::cout << "Cell P:\n"; T->cell(c_p).outputVertices(); T->cell(c_p).outputEdgeVertices();
	std::cout << "Cell Q:\n"; T->cell(c_q).outputVertices(); T->cell(c_q).outputEdgeVertices();*/
	v_a->orderCellContacts(); v_b->orderCellContacts();
	std::cout << "T1\n";
}

