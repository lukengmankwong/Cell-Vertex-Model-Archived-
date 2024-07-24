#include "cell.h"


Cell::Cell(int* vertex_counter_ptr, std::unordered_map<int, Vertex>* vertex_map_ptr, std::vector<int>& vertex_keys, 
	 int* edge_counter_ptr, std::unordered_map<int, Edge>* edge_map_ptr, std::vector<int>& edge_keys,
	 int* cell_counter_ptr, std::unordered_map<int, Cell>* cell_map_ptr) :
	 
vertex_counter_ptr(vertex_counter_ptr), vertex_map_ptr(vertex_map_ptr), edge_counter_ptr(edge_counter_ptr), edge_map_ptr(edge_map_ptr), cell_counter_ptr(cell_counter_ptr), cell_map_ptr(cell_map_ptr), id(*cell_counter_ptr)
{
	this->vertex_keys = std::unordered_set<int>(vertex_keys.begin(), vertex_keys.end());
	this->edge_keys = std::unordered_set<int>(edge_keys.begin(), edge_keys.end());
	//(*cell_counter_ptr)++;
}

Cell::~Cell() {}


bool Cell::removeEdge(int edge_id)
{
    auto it = edge_keys.find(edge_id);
    if (it != edge_keys.end()) {
		edge_keys.erase(edge_id);
		edge_map_ptr->at(edge_id).removeCellJunction();
		return true;
	} else { return false; }
}

bool Cell::removeVertex(int vertex_id)
{
    auto it = vertex_keys.find(vertex_id);
    if (it != vertex_keys.end()) {
		vertex_keys.erase(vertex_id);
		vertex_map_ptr->at(vertex_id).removeCellContact(id);
		return true;
	} else { return false; }
}


void Cell::removeVertices()
{

    for (int vertex_key : vertex_keys) {
        vertex_map_ptr->at(vertex_key).removeCellContact(id);
    }

}
void Cell::removeEdges()
{
    for (int edge_key : edge_keys) {
        edge_map_ptr->at(edge_key).removeCellJunction();
    }

}

void Cell::extrude()
{
	this->calcCentroid();
	int vertex_count = *vertex_counter_ptr;
	vertex_map_ptr->emplace(vertex_count, Vertex(vertex_counter_ptr, vertex_map_ptr, edge_map_ptr, centroid)); //create new vertex at centroid of this cell
	
	std::unordered_set<int> cell_neighbour_keys;
	for (int vertex_key : vertex_keys) { cell_neighbour_keys.insert(vertex_map_ptr->at(vertex_key).getCellContacts().begin(), vertex_map_ptr->at(vertex_key).getCellContacts().end()); }
	//find cells which share a vertex with this cell
	
	for (int cell_neighbour_key : cell_neighbour_keys) {
		for (int edge_key : edge_keys) { cell_map_ptr->at(cell_neighbour_key).removeEdge(edge_key); } //remove edges shared by cell to be extruded from other cells
	} this->removeEdges();
	
}



void Cell::calcCentroid()
{
    K::FT x_sum = 0;
    K::FT y_sum = 0;

    for (int vertex_index : vertex_keys) {
        x_sum += vertex_map_ptr->at(vertex_index).getR().x();
        y_sum += vertex_map_ptr->at(vertex_index).getR().y();
    }

    centroid = Point(x_sum/vertex_keys.size(), y_sum/vertex_keys.size());
}

void Cell::calcArea()
{
    area = 0;
}

void Cell::calcG()
{
	G[0]=0; G[1]=0;
	G[2]=0;
	
	this->calcCentroid();
	double x_0 = centroid.x(); double y_0 = centroid.y();
	for (int v : vertex_keys)
	{
		double x_v = vertex_map_ptr->at(v).getR().x();
		double y_v = vertex_map_ptr->at(v).getR().y();
		G[0]+=(x_v-x_0)*(x_v-x_0); G[1]+=(x_v-x_0)*(y_v-y_0);
		G[2]+=(x_v-x_0)*(y_v-y_0);
	}
	
	double f = 1/vertex_keys.size();
	G[0]*=f; G[1]*=f;
	G[2]*=f;
	
	lambda = 0.5*(G[0]+G[2] + std::sqrt((G[0]+G[2])*(G[0]+G[2])+4*G[1]*G[1]));
	director = Vec(1, (lambda-G[0])/G[1]);
}


const int Cell::getID() const { return id; }
const Point& Cell::getCentroid() const { return centroid; }
double Cell::getArea() const { return area; }
const std::unordered_set<int>& Cell::getVertices() const { return vertex_keys; }
const std::unordered_set<int>& Cell::getEdges() const {return edge_keys; }

