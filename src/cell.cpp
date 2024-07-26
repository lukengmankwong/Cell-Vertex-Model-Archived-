#include "cell.h"


Cell::Cell(int id, std::vector<int>& vertex_keys, std::vector<int>& edge_keys) : id(id)
{
	this->vertex_keys = std::unordered_set<int>(vertex_keys.begin(), vertex_keys.end());
	this->edge_keys = std::unordered_set<int>(edge_keys.begin(), edge_keys.end());
	A = A_0;
}

const int Cell::getID() const { return id; }
const Point& Cell::getCentroid() const { return centroid; }
const double Cell::getA() const { return A; }
const double Cell::getdA() const { return dA; }
const double Cell::getL() const { return L; }
const double Cell::getT_A() const { return T_A; }
const std::unordered_set<int>& Cell::getVertices() const { return vertex_keys; }
const std::unordered_set<int>& Cell::getEdges() const { return edge_keys; }


/*bool Cell::removeEdge(int edge_id)
{
    auto it = edge_keys.find(edge_id);
    if (it != edge_keys.end()) {
		edge_keys.erase(edge_id);
		edge_map.at(edge_id).removeCellJunction(id);
		return true;
	} else { return false; }
}
bool Cell::removeVertex(int vertex_id)
{
    auto it = vertex_keys.find(vertex_id);
    if (it != vertex_keys.end()) {
		vertex_keys.erase(vertex_id);
		vertex_map.at(vertex_id).removeCellContact(id);
		return true;
	} else { return false; }
}*/


void Cell::removeVertices()
{
    for (int v : vertex_keys) {
        vertex_map.at(v).removeCellContact(id);
    }
}
void Cell::removeEdges()
{
    for (int e : edge_keys) {
        edge_map.at(e).removeCellJunction(id);
    }
}


void Cell::extrude()
{
}


void Cell::calcCentroid()
{
    K::FT x_sum = 0;
    K::FT y_sum = 0;

    for (int v : vertex_keys) {
        x_sum += vertex_map.at(v).getR().x();
        y_sum += vertex_map.at(v).getR().y();
    }

    centroid = Point(x_sum/vertex_keys.size(), y_sum/vertex_keys.size());
}

void Cell::calcG()
{
	G[0]=0; G[1]=0;
	G[2]=0;
	
	this->calcCentroid();
	double x_0 = centroid.x(); double y_0 = centroid.y();
	for (int v : vertex_keys)
	{
		double x_v = vertex_map.at(v).getR().x();
		double y_v = vertex_map.at(v).getR().y();
		G[0]+=(x_v-x_0)*(x_v-x_0); G[1]+=(x_v-x_0)*(y_v-y_0);
		G[2]+=(x_v-x_0)*(y_v-y_0);
	}
	
	double f = 1/vertex_keys.size();
	G[0]*=f; G[1]*=f;
	G[2]*=f;
	
	lambda = 0.5*(G[0]+G[2] + std::sqrt((G[0]+G[2])*(G[0]+G[2])+4*G[1]*G[1]));
	director = Vec(1, (lambda-G[0])/G[1]);
}

void Cell::calcA()
{
	this->calcCentroid();
	double x_0 = centroid.x(); double y_0 = centroid.y();
	double A_old = A;
	A = 0;
	for (int e : edge_keys) {
		double x_1 = vertex_map.at(edge_map.at(e).getE().first).getR().x();
		double y_1 = vertex_map.at(edge_map.at(e).getE().first).getR().y();
		double x_2 = vertex_map.at(edge_map.at(e).getE().second).getR().x();
		double y_2 = vertex_map.at(edge_map.at(e).getE().second).getR().y();
		A += std::fabs( x_0*(y_1-y_2) + x_1*(y_2-y_0) + x_2*(y_0-y_1) );
	} A *= 0.5;
	dA = A - A_old;
}

void Cell::calcL()
{
	L = 0;
	for (int e : edge_keys) { //make sure edge_lengths are already calculated
		L += edge_map.at(e).getl();
	}
}

void Cell::calcT_A() { T_A = k_A*(A-A_0); }



