#include "edge.h"
#include "tissue.h"


Edge::Edge(Tissue* T, int v1, int v2) : T(T), id(T->e_c()), v_1(v1), v_2(v2) {}
//Edge::Edge(Tissue* T, Vertex* v_i_, Vertex* v_j_) : T(T), id(T->e_c()), v_i_(v_i_), v_j_(v_j_) {}
Edge::Edge() = default;

bool Edge::operator==(const Edge& other) const { return ((v_1 == other.v_1) && (v_2 == other.v_2)) || ((v_1== other.v_2) && (v_2 == other.v_1)); }

const int Edge::v1()		const { return v_1; }
const int Edge::v2() 		const { return v_2; }
const double Edge::l() 		const { return l_; }
const double Edge::T_l()	const { return T_l_; }
const std::unordered_set<int>& Edge::cellJunctions()	const { return cell_junctions; }


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


void Edge::T1() //problem with order of v_p and v_q;
{
	if (cell_junctions.size() != 2) return;
	if (T->vert(v_1).cellContacts().size() != 3 || T->vert(v_2).cellContacts().size() != 3) return;
	
	//cells either side of edge
	const int c_a = *(cell_junctions.begin());
	const int c_b = *std::next(cell_junctions.begin());
	const std::unordered_set<int> cellsAB = {c_a, c_b};
	
	//copy of edge contacts
	const std::unordered_set<int> v_1_edges = T->vert(v_1).edgeContacts();
	const std::unordered_set<int> v_2_edges = T->vert(v_2).edgeContacts();
	
	auto other_cell = [this, cellsAB](const int v) 
	{ 
		for (int c : T->vert(v).cellContacts()) if (cellsAB.find(c) == cellsAB.end()) return c; 
		return -1;
	};
	const int c_p = other_cell(v_1);
	const int c_q = other_cell(v_2);
	
	/*std::cout << "Cell P:\n"; T->cell(c_p).outputVertices(); T->cell(c_p).outputEdgeVertices();
	std::cout << "Cell Q:\n"; T->cell(c_q).outputVertices(); T->cell(c_q).outputEdgeVertices();*/
	
	//create new vertices
	Point cen = CGAL::midpoint(T->vert(v_1).r(), T->vert(v_2).r());
	Vec u = T->vert(v_2).r() - T->vert(v_1).r(); Vec s(-u.y(), u.x()); //s is u rotated 90 anticlockwise
	s *= (param::l_new/s.squared_length());
	Point a = cen + param::l_new*s; Point b = cen - param::l_new*s;
	T->cell(c_a).calcR_0(); Point r_0 = T->cell(c_a).r_0();
	if ( CGAL::squared_distance(a, r_0) > CGAL::squared_distance(b, r_0) ) std::swap(a,b);
	const int v_a = T->createVertex(a); const int v_b = T->createVertex(b);
	const int e_new = T->createEdge(v_a, v_b);
	
	//replace edge in cells a and b with vertices a and b respectively
	auto edgeToVertex = [this](const int c_x, const int v_x)
	{
		const std::vector<int>& edges = T->cell(c_x).Edges();
		for (int e : edges)
		{
			if (e != id) 
			{
				T->edge(e).swapVertex(v_1, v_x);
				T->edge(e).swapVertex(v_2, v_x);
			}
		}
		T->cellRemoveEdge(c_x, id);
	};
	edgeToVertex(c_a, v_a); edgeToVertex(c_b, v_b);
	/*std::cout << "Cell P:\n"; T->cell(c_p).outputVertices(); T->cell(c_p).outputEdgeVertices();
	std::cout << "Cell Q:\n"; T->cell(c_q).outputVertices(); T->cell(c_q).outputEdgeVertices();
	std::cout << "Cell A:\n"; T->cell(c_a).outputVertices(); T->cell(c_a).outputEdgeVertices();
	std::cout << "Cell B:\n"; T->cell(c_b).outputVertices(); T->cell(c_b).outputEdgeVertices();*/
	
	auto VertexToEdge = [this, c_a, c_b, v_a, v_b, e_new](const int c_x, const std::unordered_set<int>& v_edges)
	{
		const std::vector<int>& vertices = T->cell(c_x).Vertices();
		std::vector<int>::const_iterator it_v_a = std::find(vertices.begin(), vertices.end(), v_a);
		int i_v_a = std::distance(vertices.begin(), it_v_a);
		
		//std::cout << i_v_a << ' ' << vertices.size()-1 << '\n';
		
		int e_a = -1; int e_b = -1;
		for (int e : v_edges)
		{
			if (T->cell(c_a).hasEdge(e)) e_a = e;
			else if (T->cell(c_b).hasEdge(e)) e_b = e;
		}
		
		const std::vector<int>& edges = T->cell(c_x).Edges(); int n = edges.size();
		for (int i = 0; i < n; i++)
		{
			if (edges[i] == e_a && edges[(i+1)%n] == e_b)
			{ //correct
				T->cellNewEdge(c_x, e_new, (i+1)%n);
				T->cellNewVertex(c_x, v_b, (i+2)%n);
				break;
			}
			else if (edges[i] == e_b && edges[(i+1)%n] == e_a)
			{ //incorrect
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
	while (!(T->cell(c_p).valid())) { T->cell(c_p).rotateVertices(); }
	while (!(T->cell(c_q).valid())) { T->cell(c_q).rotateVertices(); }
	/*if (!(T->cell(c_p).valid())) { std::cout << "P invalid\n"; std::cin.get(); }
	if (!(T->cell(c_q).valid())) { std::cout << "Q invalid\n"; std::cin.get(); }
	std::cout << "Cell P:\n"; T->cell(c_p).outputVertices(); T->cell(c_p).outputEdgeVertices();
	std::cout << "Cell Q:\n"; T->cell(c_q).outputVertices(); T->cell(c_q).outputEdgeVertices();*/
	T->vert(v_a).orderCellContacts(); T->vert(v_b).orderCellContacts();
	std::cout << "T1\n";
}

