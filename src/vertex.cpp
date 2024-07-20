#include "vertex.h"


Vertex::Vertex(std::unordered_map<int, Vertex>* vertex_map, int id, Point r) :
    id(id), vertex_map(vertex_map), r(r), force(Vec(0,0)) {}

bool Vertex::operator==(const Vertex& other) const { return this->r == other.r; }

void Vertex::addCellContact(int cell_id) {
	cell_contacts.insert(cell_id);
}
void Vertex::removeCellContact(int cell_id) {
    cell_contacts.erase(cell_id);
    if (cell_contacts.size() == 0) { (*vertex_map).erase(id); } //remove vertex from vertex_map if all cells with that vertex are destroyed

}
const std::set<int>& Vertex::getCellContacts() const { return cell_contacts; }

const Point& Vertex::getR() const { return r; }


void Vertex::calcForce(Point centroid)
{
    force += (static_cast<double>(std::rand())/RAND_MAX)*k*(centroid-r);
}

void Vertex::applyForce() 
{ 
    r += a*dt*force;
    force = Vec(0,0);
}


