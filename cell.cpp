#include "cell.h"


Cell::Cell(std::vector<int>& vertices, std::vector<int>& edges) :
vertices(vertices), edges(edges) {}

void Cell::calcCentroid(std::vector<Vertex>& all_vertices)
{
    K::FT x_sum = 0;
    K::FT y_sum = 0;

    for (int vertex_id : vertices) {
        x_sum += all_vertices[vertex_id].getR().x();
        y_sum += all_vertices[vertex_id].getR().y();
    }

    centroid = Point(x_sum/vertices.size(), y_sum/vertices.size());
}

Point Cell::getCentroid() const { return centroid; }
std::vector<int>& Cell::getVertices() { return vertices; };
