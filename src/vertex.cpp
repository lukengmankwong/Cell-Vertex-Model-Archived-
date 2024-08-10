#include "vertex.h"


Vertex::Vertex(Global* g, Point r) : g(g), id(g->vertexCounter()), r(r), force(Vec(0,0)) {}

bool Vertex::operator==(const Vertex& other) const { return this->r == other.r; }

const Point& Vertex::R() const { return r; }
const int Vertex::ID() const { return id; }
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
		const std::vector<int>& c_vertex_keys = g->cellMap().at(c).getVertices();
		int n = c_vertex_keys.size(); int S = g->cellMap().at(c).getS();
		auto it = std::find(c_vertex_keys.begin(), c_vertex_keys.end(), id);
		int j = std::distance(c_vertex_keys.begin(), it);
		
		double dAdx = S*0.5*(g->vertexMap().at(c_vertex_keys[(j+1)%n]).R().y() - g->vertexMap().at(c_vertex_keys[(j-1+n)%n]).R().y());
		double dAdy = S*0.5*(g->vertexMap().at(c_vertex_keys[(j-1+n)%n]).R().x() - g->vertexMap().at(c_vertex_keys[(j+1)%n]).R().x());
		f_A -= g->cellMap().at(c).getT_A()*Vec(dAdx, dAdy);
	}
	return f_A;
}

Vec Vertex::calcLineForce()
{
	Vec f_L(0,0);
	for (int e : edge_contacts)
	{
		int v_; //other vertex in edge
		(id == g->edgeMap().at(e).getE().first) ? v_ = g->edgeMap().at(e).getE().second : v_ = g->edgeMap().at(e).getE().first;
		
		double x_diff = r.x() - g->vertexMap().at(v_).R().x();
		double y_diff = r.y() - g->vertexMap().at(v_).R().y();
		double dldx = x_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		double dldy = y_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		f_L -= g->edgeMap().at(e).getT_l()*Vec(dldx, dldy);
	}
	return f_L;
}


void Vertex::calcForce()
{
	force += calcSurfaceForce();
	force += calcLineForce();
}

void Vertex::addNoise()
{
	force += Vec(static_cast<double>(std::rand())/RAND_MAX-0.5, static_cast<double>(std::rand())/RAND_MAX-0.5);
}

void Vertex::applyForce() 
{ 
	dr = a*dt*force;
    r += dr;
}
