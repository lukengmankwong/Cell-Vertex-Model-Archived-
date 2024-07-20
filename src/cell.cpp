#include "cell.h"

Cell::Cell(int id, std::unordered_map<int, Vertex>* vertex_map, std::vector<int>& vertex_indices, std::unordered_map<int, Edge>* edge_map, std::vector<int>& edge_indices) :
id(id), vertex_map(vertex_map), vertex_indices(vertex_indices), edge_map(edge_map), edge_indices(edge_indices) {}

Cell::~Cell()
{

}


void Cell::kill()
{
    for (int edge_index : edge_indices) {
        (*edge_map).at(edge_index).removeCellJunction();
    }
    for (int edge_index : edge_indices) {
        if ((*edge_map).at(edge_index).getCellJunctions() == 0) {
            (*edge_map).erase(edge_index);
        }
    }

    for (int vertex_index : vertex_indices) {
        std::cout << vertex_index << ' ';
        (*vertex_map).at(vertex_index).removeCellContact(id);
    } std::cout << '\n';

}

void Cell::extrude()
{

}


void Cell::calcCentroid()
{
    K::FT x_sum = 0;
    K::FT y_sum = 0;

    for (int vertex_index : vertex_indices) {
        x_sum += (*vertex_map).at(vertex_index).getR().x();
        y_sum += (*vertex_map).at(vertex_index).getR().y();
    }

    centroid = Point(x_sum/vertex_indices.size(), y_sum/vertex_indices.size());
}

void Cell::calcArea()
{
    area = 0;
    for (int i = 0; i < vertex_indices.size(); i++) {
        int j1 = vertex_indices[i]; int j2 = vertex_indices[(i+1)%vertex_indices.size()];
        area += vertex_map->at(j1).getR().x()*vertex_map->at(j2).getR().y()
            - vertex_map->at(j2).getR().x()*vertex_map->at(j1).getR().y();
    }
    area = 0.5*std::abs(area);
}

const int Cell::getID() const { return id; }
const Point& Cell::getCentroid() const { return centroid; }
double Cell::getArea() const { return area; }
const std::vector<int>& Cell::getVertices() const { return vertex_indices; }
const std::vector<int>& Cell::getEdges() const {return edge_indices; }

