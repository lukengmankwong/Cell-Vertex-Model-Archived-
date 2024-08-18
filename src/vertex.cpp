#include "vertex.h"


Vertex::Vertex(Global* g, Point r) : g(g), id(g->vertexCounter()), r(r), force(Vec(0,0)) {}

bool Vertex::operator==(const Vertex& other) const { return r == other.r; }

const Point& Vertex::R() const { return r; }
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
		int n = c_vertices.size(); int S = g->cell(c).getS();
		auto it = std::find(c_vertices.begin(), c_vertices.end(), id);
		int j = std::distance(c_vertices.begin(), it);
		
		double dAdx = S*0.5*(g->vert(c_vertices[(j+1)%n]).R().y() - g->vert(c_vertices[(j-1+n)%n]).R().y());
		double dAdy = S*0.5*(g->vert(c_vertices[(j-1+n)%n]).R().x() - g->vert(c_vertices[(j+1)%n]).R().x());
		f_A -= g->cell(c).getT_A()*Vec(dAdx, dAdy);
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
		
		double x_diff = r.x() - g->vert(v_).R().x();
		double y_diff = r.y() - g->vert(v_).R().y();
		double dldx = x_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		double dldy = y_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		f_L -= g->edge(e).getT_l()*Vec(dldx, dldy);
	}
	return f_L;
}

void Vertex::calcForce() { force += calcSurfaceForce()+calcLineForce(); }
void Vertex::applyForce() { r += a*dt*force; }


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
	return contact_order;
}

void Vertex::T1()
{
	//no transition unless vertex is fourfold and not on the boundary
	
	if (cell_contacts.size() !=4) { return; }
	for (int c : cell_contacts) { if (g->cell(c).onBoundary()) return; }

	std::cout << "vertex: " << id << '\n';
	std::vector<int> contact_order = orderCellContacts();
	int c_a = contact_order[0]; int c_b = contact_order[2];
	int c_p = contact_order[1]; int c_q = contact_order[3];
	
	std::cout << "A, B, BEFORE:\n";
	g->cell(c_a).outputVertices(); g->cell(c_b).outputVertices();
	std::cout << '\n';
	std::cout << "P, Q, BEFORE:\n";
	g->cell(c_p).outputVertices(); g->cell(c_q).outputVertices();
	std::cout << '\n';

	//it_x points to this vertex, it_x1 the vertex before and it_x2 the vertex after
	std::vector<int>::const_iterator it_a, it_a1, it_a2, it_b, it_b1, it_b2;	
	auto vertexIterators = [this](const int c_x, std::vector<int>::const_iterator& it_x, std::vector<int>::const_iterator& it_x1, std::vector<int>::const_iterator& it_x2)
	{
		const std::vector<int>& vertices = g->cell(c_x).Vertices();
		it_x = std::find(vertices.begin(), vertices.end(), id);
		if (it_x == vertices.begin()) it_x1 = vertices.end()-1;
		else it_x1 = it_x-1;
		it_x2 = it_x + ((vertices.size()+1) % vertices.size());
	};	
	vertexIterators(c_a, it_a, it_a1, it_a2); vertexIterators(c_b, it_b, it_b1, it_b2);
 	
	//midpoints of two edges connected to this vertex for cells a and b
	Point a = CGAL::midpoint(g->vert(*it_a1).R(), g->vert(*it_a2).R());
	Point b = CGAL::midpoint(g->vert(*it_b1).R(), g->vert(*it_b2).R());
	int v_a = g->createVertex(a); int v_b = g->createVertex(b);
	std::cout << "v_a: " << v_a << " v_b: " << v_b << '\n';
	
	auto updateEdgeVertices = [this](const int c_x, const int v_x, std::vector<int>::const_iterator& it_x, std::vector<int>::const_iterator& it_x1)
	{
		const std::vector<int>& vertices = g->cell(c_x).Vertices();
		int i = std::distance(vertices.begin(), it_x);
		int i1 = std::distance(vertices.begin(), it_x1);
		
		const std::vector<int>& edges = g->cell(c_x).Edges(); 	//need to check that edges and vertices have same starting point after cell division
		g->edge(edges[i]).swapVertex(id, v_x);
		g->edge(edges[i1]).swapVertex(id, v_x);				
	};
	addEdgeContact(-1); //temporary to prevent this vertex from being destroyed, there is no edge -1
	updateEdgeVertices(c_a, v_a, it_a, it_a1); updateEdgeVertices(c_b, v_b, it_b, it_b1);

	int e_new = g->createEdge(v_a, v_b);
	int j_p, j_q;
	//update vertices and edges for cells p and q
	auto updateVerticesEdges = [this, v_a, v_b, e_new](const int c_x, int& j_x)
	{
		const std::vector<int>& vertices = g->cell(c_x).Vertices();
		const std::vector<int>& edges = g->cell(c_x).Edges(); 
		//find edge between cell A and cell P etc
		int i_a, i_b;
		for (int i = 0; i < edges.size(); i++)
		{
			if (g->edge(edges[i]).hasVertex(v_a)) i_a = i;
			if (g->edge(edges[i]).hasVertex(v_b)) i_b = i;
		}
		if (i_a < i_b && i_b != edges.size()-1)  		j_x = i_b;
		else if (i_a < i_b && (i_b == edges.size()-1) )	j_x = 0;
		else if (i_b < i_a && (i_a != edges.size()-1) )	j_x = i_a;
		else if (i_b < i_a && (i_a == edges.size()-1) )	j_x = 0;
		g->cellNewEdge(c_x, e_new, j_x);	
		
		
		std::vector<int>::const_iterator it_va = std::find(vertices.begin(), vertices.end(), v_a);
		std::vector<int>::const_iterator it_vb = std::find(vertices.begin(), vertices.end(), v_b);
		int i;
		if (it_va == vertices.end()) //v_a missing from vertices
		{
			if (g->edge(edges[(j_x+1)%edges.size()]).hasVertex(v_a)) i = (j_x+1)%edges.size();
			else i = (j_x+edges.size()-1)%edges.size();
			g->cellNewVertex(c_x, v_a, i);
		}
		else if (it_vb == vertices.end()) //v_b missing from vertices
		{
			if (g->edge(edges[(j_x+1)%edges.size()]).hasVertex(v_b)) i = (j_x+1)%edges.size();
			else i = (j_x+edges.size()-1)%edges.size();
			g->cellNewVertex(c_x, v_b, i);
		}			
	};	
	updateVerticesEdges(c_p, j_p); updateVerticesEdges(c_q, j_q);
	
	std::cout << "A, B, AFTER:\n";
	g->cell(c_a).outputVertices(); g->cell(c_b).outputVertices();
	std::cout << '\n';
	std::cout << "P, Q, AFTER:\n";
	g->cell(c_p).outputVertices(); g->cell(c_q).outputVertices();
	g->cell(c_p).outputEdgeVertices(); g->cell(c_q).outputEdgeVertices();
	std::cout << '\n';

	removeEdgeContact(-1);
}



