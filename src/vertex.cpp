#include "vertex.h"


Vertex::Vertex(int id, Point r) : id(id), r(r), force(Vec(0,0)) { }

bool Vertex::operator==(const Vertex& other) const { return this->r == other.r; }

const Point& Vertex::getR() const { return r; }
const int Vertex::getID() const { return id; }

const std::unordered_set<int>& Vertex::getCellContacts() const { return cell_contacts; }


void Vertex::addCellContact(int cell_id) { cell_contacts.insert(cell_id); }

void Vertex::removeCellContact(int cell_id) {
    cell_contacts.erase(cell_id);
    if (cell_contacts.size() == 0) { 
		for (int e : edge_contacts) { edge_map.erase(e); }
		vertex_map.erase(id); 
	} //remove vertex from vertex_map and edges with that vertex if no cells have the vertex
}


void Vertex::addEdgeContact(int edge_id) { edge_contacts.insert(edge_id); }

void Vertex::removeEdgeContact(int edge_id) { edge_contacts.erase(edge_id); }


bool Vertex::inEdge(int edge_index) const { return edge_map.at(edge_index).getE().first == id || edge_map.at(edge_index).getE().second == id; }


Vec Vertex::calcSurfaceForce()
{
	Vec f_A(0,0);
	for (int c : cell_contacts) 
	{
		const std::vector<int>& c_vertex_keys = cell_map.at(c).getVertices();
		int n = c_vertex_keys.size();
		auto it = std::find(c_vertex_keys.begin(), c_vertex_keys.end(), id);
		int j = std::distance(c_vertex_keys.begin(), it);
		
		double dAdx = 0.5*(vertex_map.at(c_vertex_keys[(j+1)%n]).getR().y() - vertex_map.at(c_vertex_keys[(j-1+n)%n]).getR().y());
		double dAdy = 0.5*(vertex_map.at(c_vertex_keys[(j-1+n)%n]).getR().x() - vertex_map.at(c_vertex_keys[(j+1)%n]).getR().x());
		f_A -= cell_map.at(c).getT_A()*Vec(dAdx, dAdy);
	}
	return f_A;
}

Vec Vertex::calcLineForce()
{
	Vec f_L(0,0);
	for (int e : edge_contacts)
	{
		int v_;
		(id == edge_map.at(e).getE().first) ? v_ = edge_map.at(e).getE().second : v_ = edge_map.at(e).getE().first;
		
		double x_diff = r.x() - vertex_map.at(v_).getR().x();
		double y_diff = r.y() - vertex_map.at(v_).getR().y();
		double dldx = x_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		double dldy = y_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		f_L -= edge_map.at(e).getT_l()*Vec(dldx, dldy);
	}
	return f_L;
}


void Vertex::calcForce(Point centroid)
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




