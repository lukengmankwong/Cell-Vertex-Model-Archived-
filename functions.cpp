#include "functions.h"


void getInitialData(VD& vd,  std::vector<Cell>& cells, std::vector<Edge>& edges, std::vector<Vertex>& vertices)
{
	for (VD::Vertex_iterator vit = vd.vertices_begin(); vit != vd.vertices_end(); vit++) { //iterate through vertices in VD container and create vertex objects from points
        Vertex vertex(vit->point());
        vertices.push_back(vertex);
    }
    
    for (VD::Face_iterator fi = vd.faces_begin(); fi != vd.faces_end(); fi++) //iterate through faces in VD container
    {
        std::vector<int> cell_vertices;
        VD::Ccb_halfedge_circulator ec_start = fi->ccb();
        VD::Ccb_halfedge_circulator ec = ec_start;
        do { //iterate through vertices of face
                if (!ec->is_unbounded()) {
                    Vertex vertex(ec->source()->point()); //create Vertex object for the given vertex
                    auto it = std::find(vertices.begin(), vertices.end(), vertex); 
                    cell_vertices.push_back(std::distance(vertices.begin(), it));
                } //find index of vertex in vertices vector and push back to cell_vertices
                ++ec;
            } while (ec != ec_start);

        std::vector<int> cell_edges; 
        for (int i = 0; i < cell_vertices.size(); i++) {
            Edge edge(cell_vertices[i%cell_vertices.size()], cell_vertices[(i+1)%cell_vertices.size()]); 
            auto it = std::find(edges.begin(), edges.end(), edge);
            if (it == edges.end()) { 
                edges.push_back(edge); //if edge has not been added to edges yet we add it
                cell_edges.push_back(edges.size()-1); //add index of added edge to cell_edges
            }
            else { cell_edges.push_back(std::distance(edges.begin(), it)); } 
        }

        cells.push_back(Cell(cell_vertices, cell_edges));
    }

}


void runSimulation(std::vector<Cell> &cells, std::vector<Edge>& edges, std::vector<Vertex> &vertices, int time_steps)
{
    for (int step = 0; step < time_steps; step++)
    {
        for (Cell& cell : cells) { //calculate forces
            cell.calcCentroid(vertices);
            for (int vertex_id : cell.getVertices()) {
                vertices[vertex_id].calcForce(cell.getCentroid());
            }
        }
        for (Vertex& vertex : vertices) {
            vertex.applyForce();
        }
        std::string filename = "graph" + std::to_string(step) + ".vtk";
        WriteVTKFile(vertices, edges, filename);
    }
}



void outputData(const std::vector<Vertex>& vertices, const std::vector<Cell> &cells)
{
    for (Cell cell : cells) {
    std::cout << "Cell vertices: ";
    for (int vertex_id : cell.getVertices()) {
        std::cout << vertices[vertex_id].getR() << ", ";
    }
    std::cout << '\n';
    }
}

void WriteVTKFile(const std::vector<Vertex>& vertices, const std::vector<Edge>& edges, const std::string& filename) {

    std::ofstream vtkFile(filename);
    vtkFile << "# vtk DataFile Version 3.0\nGraph\nASCII\nDATASET POLYDATA\n";

    vtkFile << "POINTS " << vertices.size() << " float" << std::endl;
    for (const auto& vertex : vertices) {
        vtkFile << vertex.getR().x() << " " << vertex.getR().y() << " " << "0" << std::endl;
    }


    vtkFile << "LINES " << edges.size() << " " << edges.size() * 3 << std::endl;
    for (const auto& edge : edges) {
        vtkFile << "2 " << edge.getE().first << " " << edge.getE().second << std::endl;
    }

    vtkFile.close();
}