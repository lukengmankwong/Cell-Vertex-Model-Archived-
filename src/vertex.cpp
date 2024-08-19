#include "vertex.h"


Vertex::Vertex(Global* g, Point r) : g(g), id(g->v_c()), r_(r), force_(Vec(0,0)) {}

bool Vertex::operator==(const Vertex& other) const { return r_ == other.r_; }

const Point& Vertex::r() const { return r_; }
const double Vertex::m() const { return m_; }
const std::unordered_set<int>& Vertex::cellContacts() const { return cell_contacts; }
const std::unordered_set<int>& Vertex::edgeContacts() const { return edge_contacts; }

void Vertex::addCellContact(int cell_id) { cell_contacts.insert(cell_id); }
void Vertex::removeCellContact(int cell_id) { cell_contacts.erase(cell_id); }

void Vertex::addEdgeContact(int edge_id) { edge_contacts.insert(edge_id); }
void Vertex::removeEdgeContact(int edge_id) 
{
	edge_contacts.erase(edge_id);
	if (edge_contacts.size() == 0) { g->vertexMap().erase(id); }
}


Vec Vertex::calcSurfaceForce()
{
	Vec f_A(0,0);
	for (int c : cell_contacts) 
	{
		const std::vector<int>& c_vertices = g->cell(c).Vertices();
		int n = c_vertices.size(); int S = g->cell(c).S();
		auto it = std::find(c_vertices.begin(), c_vertices.end(), id);
		int j = std::distance(c_vertices.begin(), it);
		
		double dAdx = S*0.5*(g->vert(c_vertices[(j+1)%n]).r().y() - g->vert(c_vertices[(j-1+n)%n]).r().y());
		double dAdy = S*0.5*(g->vert(c_vertices[(j-1+n)%n]).r().x() - g->vert(c_vertices[(j+1)%n]).r().x());
		f_A -= g->cell(c).T_A()*Vec(dAdx, dAdy);
	}
	return f_A;
}
Vec Vertex::calcLineForce()
{
	Vec f_L(0,0);
	for (int e : edge_contacts)
	{
		int v_; //other vertex in edge
		(id == g->edge(e).v1()) ? v_ = g->edge(e).v2() : v_ = g->edge(e).v1();
		
		double x_diff = r_.x() - g->vert(v_).r().x();
		double y_diff = r_.y() - g->vert(v_).r().y();
		double dldx = x_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		double dldy = y_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		f_L -= g->edge(e).T_l()*Vec(dldx, dldy);
	}
	return f_L;
}

void Vertex::calcForce() { force_ += calcSurfaceForce()+calcLineForce(); }
void Vertex::applyForce() { r_ += a*dt*force_; }


std::vector<int> Vertex::orderCellContacts()
{
	std::vector<int> contact_order;
	for (int c : cell_contacts) { if (g->cell(c).onBoundary()) return contact_order; }
	contact_order.push_back(*(cell_contacts.begin()));
	std::unordered_set<int> edge_set = edge_contacts; int i = 0;
	while (contact_order.size() < cell_contacts.size())
	{
		for (int e : edge_set)
		{
			const std::unordered_set<int>& cell_junctions = g->edge(e).cellJunctions();
			std::unordered_set<int>::const_iterator c_it = std::find(cell_junctions.begin(), cell_junctions.end(), contact_order[i]);
			if (c_it != cell_junctions.end())
			{
				if (c_it == cell_junctions.begin()) contact_order.push_back(*(std::next(c_it)));
				else contact_order.push_back(*(cell_junctions.begin()));
				i++; edge_set.erase(e); break;
			}
		}
	}
	//ensure anticlockwise orientation
	for (int c : cell_contacts) g->cell(c).calcR_0();
	double A = 0;
	for (int i = 0; i < contact_order.size(); i++) 
	{	A += g->cell(contact_order[i]).r_0().x()*g->cell(contact_order[(i+1)%contact_order.size()]).r_0().y()
			- g->cell(contact_order[(i+1)%contact_order.size()]).r_0().x()*g->cell(contact_order[i]).r_0().y(); } 
	if (A < 0) std::reverse(contact_order.begin(), contact_order.end());
	return contact_order;
	
}


void Vertex::T1split()
{
	if (cell_contacts.size() != 4 || edge_contacts.size() != 4) return;
	for (int c : cell_contacts) if (g->cell(c).onBoundary()) return; 
	
	//std::cout << "vertex id: " << id << '\n';
	//affected cells, a,b change vertex p,q gets new edge
	std::vector<int> nn_cells = orderCellContacts();
	const int c_a = nn_cells[0]; const int c_b = nn_cells[2];
	const int c_p = nn_cells[1]; const int c_q = nn_cells[3];
	/*std::cout << "BEFORE\n";
	g->cell(c_p).outputVertices(); g->cell(c_p).outputEdgeVertices();	
	g->cell(c_q).outputVertices(); g->cell(c_q).outputEdgeVertices();*/

	const std::vector<int>& vertices_p = g->cell(c_p).Vertices();
	auto it_id = std::find(vertices_p.begin(), vertices_p.end(), id);
	int i = std::distance(vertices_p.begin(), it_id);
	/*std::cout << "MIDDLE\n";
	g->cell(c_p).outputVertices(); g->cell(c_p).outputEdgeVertices();
	g->cell(c_q).outputVertices(); g->cell(c_q).outputEdgeVertices();*/

	int v_a, v_b; //new vertices
	auto updateAB = [this](const int c_x, int& v_x)
	{
		//find position of and create vertex
		g->cell(c_x).calcR_0();
		Vec vec_x = g->cell(c_x).r_0() - r_;
		Point x = r_ + 0.1*vec_x;
		v_x = g->createVertex(x);

		//attatch relevent edges to vertex
		for (int e : g->cell(c_x).Edges()) g->edge(e).swapVertex_rep(id, v_x);
		if (!(g->cell(c_x).valid())) g->cell(c_x).rotateVertices();
	};
	updateAB(c_a, v_a); updateAB(c_b, v_b);

	//edge that vertex is split into
	const int e_new = g->createEdge(v_a, v_b);

	//update edges for cells p, q
	auto updateEdges = [this, c_a, c_b, v_a, v_b, e_new](const int c_x)
	{
		const std::vector<int>& vertices_c_x = g->cell(c_x).Vertices();
		const std::vector<int>& edges_c_x = g->cell(c_x).Edges();
		for (int i = 0; i < edges_c_x.size(); i++)
		{
			if (g->edge(edges_c_x[i]).hasVertex(v_a) && g->edge(edges_c_x[(i+1)%edges_c_x.size()]).hasVertex(v_b))
			{
				g->cellNewEdge(c_x,e_new, i+1);
				break;
			}
			else if (g->edge(edges_c_x[i]).hasVertex(v_b) && g->edge(edges_c_x[(i+1)%edges_c_x.size()]).hasVertex(v_a))
			{
				g->cellNewEdge(c_x,e_new, i+1);
				break;
			}
		}		
	};
	updateEdges(c_p); updateEdges(c_q);
	
	auto updateVertices = [this, c_a, c_b, v_a, v_b, e_new](const int c_x, bool before)
	{
		const std::vector<int>& vertices = g->cell(c_x).Vertices();
		std::vector<int>::const_iterator it_va = std::find(vertices.begin(), vertices.end(), v_a);
		
		if (it_va != vertices.end())
		{
			int i = std::distance(vertices.begin(), it_va);
			if (before) g->cellNewVertex(c_x, v_b, i);
			else g->cellNewVertex(c_x, v_b, i+1);
			while (!(g->cell(c_x).valid())) g->cell(c_x).rotateVertices();
		}
	};
	updateVertices(c_p, true); updateVertices(c_q, false);
	/*if (!(g->cell(c_a).valid())) { std::cout << "a invalid\n"; std::cin.get(); }
	if (!(g->cell(c_b).valid())) { std::cout << "b invalid\n"; std::cin.get(); }
	if (!(g->cell(c_p).valid())) { std::cout << "p invalid\n"; std::cin.get(); }
	if (!(g->cell(c_q).valid())) { std::cout << "q invalid\n"; std::cin.get(); }
	
	std::cout << "END\n";
	g->cell(c_p).outputVertices(); g->cell(c_p).outputEdgeVertices(); 
	g->cell(c_q).outputVertices(); g->cell(c_q).outputEdgeVertices();*/
	g->destroyVertex(id);
}


void Vertex::calcm()
{
	std::vector<int> nn_cells = orderCellContacts(); int n = nn_cells.size();
	double w = 0;
	for (int i = 0; i < nn_cells.size(); i++) w += g->D_angle(nn_cells[i], nn_cells[(i+1)%n]);
	m_ = 0.5*w/pi;
}


