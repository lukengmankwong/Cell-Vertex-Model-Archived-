#include "vertex.h"


Vertex::Vertex(Global* g, Point r) : g(g), id(g->vertexCounter()), r(r), force(Vec(0,0)) {}

bool Vertex::operator==(const Vertex& other) const { return this->r == other.r; }

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
		(id == g->edge(e).E().first) ? v_ = g->edge(e).E().second : v_ = g->edge(e).E().first;
		
		double x_diff = r.x() - g->vert(v_).R().x();
		double y_diff = r.y() - g->vert(v_).R().y();
		double dldx = x_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		double dldy = y_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		f_L -= g->edge(e).getT_l()*Vec(dldx, dldy);
	}
	return f_L;
}


void Vertex::calcForce()
{
	force += calcSurfaceForce();
	force += calcLineForce();
}

void Vertex::applyForce() 
{ 
	dr = a*dt*force;
    r += dr;
}

void Vertex::T1()
{
	if (edge_contacts.size()<4) { return; }
	
	
	
	
}



