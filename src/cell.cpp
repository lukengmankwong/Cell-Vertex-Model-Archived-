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
	/*this->calcCentroid();
	int vertex_count = *vertex_counter_ptr;
	vertex_map_ptr->emplace(vertex_count, Vertex(vertex_counter_ptr, vertex_map_ptr, edge_map_ptr, centroid));*/
	
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
	//unfinished definition
    area = 0.5*std::abs(area);
}

const int Cell::getID() const { return id; }
const Point& Cell::getCentroid() const { return centroid; }
double Cell::getArea() const { return area; }
const std::unordered_set<int>& Cell::getVertices() const { return vertex_keys; }
const std::unordered_set<int>& Cell::getEdges() const {return edge_keys; }

