#include "vertex.h"


Vertex::Vertex(int id, int* vertex_counter_ptr, std::unordered_map<int, Vertex>* vertex_map_ptr, std::unordered_map<int, Edge>* edge_map_ptr, Point r) :
    vertex_counter_ptr(vertex_counter_ptr), vertex_map_ptr(vertex_map_ptr), edge_map_ptr(edge_map_ptr), r(r), id(id), force(Vec(0,0)) 
    { 
		//(*vertex_counter_ptr)++; 
	}

bool Vertex::operator==(const Vertex& other) const { return this->r == other.r; }

void Vertex::addCellContact(int cell_id) { cell_contacts.insert(cell_id); }
void Vertex::removeCellContact(int cell_id) {
    cell_contacts.erase(cell_id);
    if (cell_contacts.size() == 0) { vertex_map_ptr->erase(id); } //remove vertex from vertex_map if all cells with that vertex are destroyed
}
const std::unordered_set<int>& Vertex::getCellContacts() const { return cell_contacts; }

void Vertex::addEdgeContact(int edge_id) { edge_contacts.insert(edge_id); }
void Vertex::removeEdgeContact(int edge_id) { 
	edge_contacts.erase(edge_id);
}

const Point& Vertex::getR() const { return r; }
const int Vertex::getID() const { return id; }

bool Vertex::inEdge(int edge_index) const { return edge_map_ptr->at(edge_index).getE().first == id || edge_map_ptr->at(edge_index).getE().second == id; }


void Vertex::calcForce(Point centroid)
{
    force += (static_cast<double>(std::rand())/RAND_MAX)*k*(centroid-r);
}

void Vertex::applyForce() 
{ 
    r += a*dt*force;
    force = Vec(0,0);
}


