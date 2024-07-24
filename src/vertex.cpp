#include "vertex.h"


Vertex::Vertex(int id, Point r) : id(id), r(r), force(Vec(0,0)) 
{ 
	dr = Vec(0.01, 0.01);
}

bool Vertex::operator==(const Vertex& other) const { return this->r == other.r; }

const Point& Vertex::getR() const { return r; }
const int Vertex::getID() const { return id; }


void Vertex::addCellContact(int cell_id) { cell_contacts.insert(cell_id); }

void Vertex::removeCellContact(int cell_id) {
    cell_contacts.erase(cell_id);
    if (cell_contacts.size() == 0) { vertex_map.erase(id); } //remove vertex from vertex_map if all cells with that vertex are destroyed
}

const std::unordered_set<int>& Vertex::getCellContacts() const { return cell_contacts; }


void Vertex::addEdgeContact(int edge_id) { edge_contacts.insert(edge_id); }

void Vertex::removeEdgeContact(int edge_id) { //possibly unfinished definition
	edge_contacts.erase(edge_id);
}

bool Vertex::inEdge(int edge_index) const { return edge_map.at(edge_index).getE().first == id || edge_map.at(edge_index).getE().second == id; }


void Vertex::calcForce(Point centroid)
{
    for (int c : cell_contacts) {
		double T = cell_map.at(c).getT();
		force += Vec(-T*cell_map.at(c).getdA()/dr.x(), -T*cell_map.at(c).getdA()/dr.y());
	}
}

void Vertex::applyForce() 
{ 
	dr = a*dt*force;
    r += dr;
    force = Vec(0,0);
}
