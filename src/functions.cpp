#include "functions.h"


void getInitialData(VD& vd,  std::unordered_map<int, Vertex>& vertex_map, std::unordered_map<int, Edge>& edge_map, std::unordered_map<int, Cell>& cell_map, int& vertex_counter, int& edge_counter, int& cell_counter)
{
    std::cout << "COLLECTING INITIAL DATA\n";
	for (VD::Vertex_iterator vit = vd.vertices_begin(); vit != vd.vertices_end(); vit++) { //iterate through vertices in VD container and create vertex objects from points
        vertex_map.emplace( vertex_counter, Vertex(&vertex_counter, &vertex_map, &edge_map, vit->point()) ); //save vertex to vertex_map
        vertex_counter++;

    }

    for (VD::Face_iterator fi = vd.faces_begin(); fi != vd.faces_end(); fi++) //iterate through faces in VD container
    {
        std::vector<int> cell_vertex_indices; //vertices of given face
        VD::Ccb_halfedge_circulator ec_start = fi->ccb();
        VD::Ccb_halfedge_circulator ec = ec_start;
        do { //iterate through vertices of face/cell
                if (!ec->is_unbounded()) {
                    for (auto& vertex : vertex_map) { //find vertex in vertex_map that is equal to current vertex
                        if (ec->source()->point() == vertex.second.getR()) {
							cell_vertex_indices.push_back(vertex.first); //tell cell that it has this vertex
                            vertex.second.addCellContact(cell_counter); //tell vertex it is part of this face/cell
                            break;
                        }
                    }
                }
                ++ec;
        } while (ec != ec_start);

        std::vector<int> cell_edge_indices; //edges of given cell/face
        for (int i = 0; i < cell_vertex_indices.size(); i++)  //iterate through cell vertices
        {
            Edge edge(&edge_map, edge_counter, cell_vertex_indices[i], cell_vertex_indices[(i+1)%cell_vertex_indices.size()]); //edge with two vertices from cell
            bool new_edge = true;
            for (auto& edge_pair : edge_map) { //iterate through found edges
                if (edge == edge_pair.second) {
					cell_edge_indices.push_back(edge_pair.first); //tell cell that it has this edge
                    edge_pair.second.addCellJunction(); //tell edge that is is part of a new cell
                    vertex_map.at(cell_vertex_indices[i]).addEdgeContact(edge_pair.first); //tell vertices of this edge that they are part of this edge
                    vertex_map.at(cell_vertex_indices[(i+1)%cell_vertex_indices.size()]).addEdgeContact(edge_pair.first);
                    new_edge = false;
                    break;
                }
            }
            if (new_edge) {
                edge_map.emplace(edge_counter, edge);
                edge_map.at(edge_counter).addCellJunction();
                cell_edge_indices.push_back(edge_counter);
                vertex_map.at(cell_vertex_indices[i]).addEdgeContact(edge_counter);
                vertex_map.at(cell_vertex_indices[(i+1)%cell_vertex_indices.size()]).addEdgeContact(edge_counter);
                edge_counter++;
            }
        }

        cell_map.emplace(cell_counter, Cell(&vertex_counter, &vertex_map, cell_vertex_indices, &edge_counter, &edge_map, cell_edge_indices, &cell_counter, &cell_map));
        cell_counter++;
    }

	std::set<int> cells_to_remove;
	for (const auto& cell : cell_map) {
        if (cell.second.getVertices().size() <= 2) { cells_to_remove.insert(cell.first); }
    }
	for (const auto& vertex : vertex_map) {
		if ((vertex.second.getR().x()*vertex.second.getR().x() + vertex.second.getR().y()*vertex.second.getR().y() < 0.075) || (vertex.second.getR().x()*vertex.second.getR().x() + vertex.second.getR().y()*vertex.second.getR().y() > 0.25) ) {
			cells_to_remove.insert(vertex.second.getCellContacts().begin(), vertex.second.getCellContacts().end());
		}
	}
    for (int c : cells_to_remove) { cell_map.at(c).removeEdges(); }
	for (int c : cells_to_remove) { cell_map.at(c).removeVertices(); }
	for (int c : cells_to_remove) { cell_map.erase(c); }

}


void runSimulation(std::unordered_map<int, Vertex>& vertex_map, std::unordered_map<int, Edge>& edge_map, std::unordered_map<int, Cell>& cell_map, int time_steps)
{
    for (int step = 0; step < time_steps; step++)
    { 
		WriteVTKFile(vertex_map, edge_map, cell_map, "graph" + std::to_string(step) + ".vtk");
        for (auto& cell : cell_map) { //calculate forces
            cell.second.calcCentroid();
            for (int vertex_id : cell.second.getVertices()) {
                vertex_map.at(vertex_id).calcForce(cell.second.getCentroid());
            }
        }
        for (auto& vertex : vertex_map) {
            vertex.second.applyForce();
        }
        /*for (auto& cell : cell_map) {
			cell.second.calcArea();
			if (cell.first == 50) {
				cell.second.extrude();
			}
		}*/
		for (auto& cell : cell_map) { cell.second.calcG(); }
    }
    WriteVTKFile(vertex_map, edge_map, cell_map, "graph" + std::to_string(timesteps) + ".vtk");
}




void outputData(const std::unordered_map<int, Vertex>& vertex_map, const std::unordered_map<int, Edge>& edge_map, const std::unordered_map<int, Cell>& cell_map)
{
    for (const auto& cell : cell_map) 
    {
        std::cout << "Cell (" << cell.first << ") : ";
        for (int vertex_id : cell.second.getVertices()) { std::cout << vertex_map.at(vertex_id).getR() << ", "; }
        std::cout << '\n';
        for (int edge_id : cell.second.getEdges()) { std::cout << edge_map.at(edge_id).getE().first << ' ' << edge_map.at(edge_id).getE().second << '\n'; }
        std::cout << '\n';
    }
}


void WriteVTKFile(const std::unordered_map<int, Vertex>& vertex_map, const std::unordered_map<int, Edge>& edge_map, const std::unordered_map<int, Cell>& cell_map, const std::string& filename) {

    // Create a mapping from the original vertex IDs to zero-based indices
    std::unordered_map<int, int> index_map;
    int index = 0;
    for (const auto& vertex : vertex_map) { index_map[vertex.first] = index++; }

    std::ofstream vtkFile(filename);
    vtkFile << "# vtk DataFile Version 3.0\nGraph\nASCII\nDATASET POLYDATA\nPOINTS " << vertex_map.size() << " float\n";
    for (const auto& vertex : vertex_map) { vtkFile << vertex.second.getR().x() << " " << vertex.second.getR().y() << " " << "0\n"; }

    vtkFile << "LINES " << edge_map.size() << " " << edge_map.size() * 3 << '\n';
    for (const auto& edge : edge_map) { vtkFile << "2 " << index_map[edge.second.getE().first] << " " << index_map[edge.second.getE().second] << '\n'; }   
	
    vtkFile.close();
}
