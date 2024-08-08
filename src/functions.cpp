#include "functions.h"


void getInitialData(VD& vd, Global& global, bool (*in)(const Point&))
{
    std::cout << "COLLECTING INITIAL DATA\n";
	for (VD::Vertex_iterator vit = vd.vertices_begin(); vit != vd.vertices_end(); vit++)
	{       
        global.createVertex(vit->point());
    }
    
    for (VD::Face_iterator fi = vd.faces_begin(); fi != vd.faces_end(); fi++) 
    {
        std::vector<int> cell_vertex_indices;
        VD::Ccb_halfedge_circulator ec_start = fi->ccb();
        VD::Ccb_halfedge_circulator ec = ec_start;

        do { 
			if (!ec->is_unbounded()) 
			{
				for (const auto& vertex : global.vertexMap())
				{
					if (ec->source()->point() == vertex.second.getR())
					{
						cell_vertex_indices.push_back(vertex.first);
						break;
					}
				}
			}
			++ec;
        } while (ec != ec_start);
		
        std::vector<int> cell_edge_indices; 
        if (cell_vertex_indices.size() >= 3)
        {
			for (int i = 0; i < cell_vertex_indices.size(); i++)
			{
				int v1 = cell_vertex_indices[i]; int v2 = cell_vertex_indices[(i+1)%cell_vertex_indices.size()];
				bool found = false;
				for (const auto& edge : global.edgeMap())
				{
					if ( (edge.second.getE().first == v1 && edge.second.getE().second == v2) || (edge.second.getE().first == v2 && edge.second.getE().second == v1) )
					{
						cell_edge_indices.push_back(edge.first);
						found = true;
						break;
					}
				}
				if (!found)
				{
					cell_edge_indices.push_back(global.createEdge(v1, v2));
				}		
			}
		}
        global.createCell(cell_vertex_indices, cell_edge_indices);
    } 
    
	std::unordered_set<int> cells_to_remove;
	for (const auto& vertex : global.vertexMap())  
	{ 
		if ( !in(vertex.second.getR()) ) 
		{ 
			cells_to_remove.insert(vertex.second.getCellContacts().begin(), vertex.second.getCellContacts().end()); 
		}
	}
	for (int c : cells_to_remove) { global.destroyCell(c); }
	int V = global.vertexMap().size(); int E = global.edgeMap().size(); int C = global.cellMap().size();
	int Euler = V-E+C;
    std::cout << "V=" << V << "\nE=" << E << "\nC=" << C << "\nV-E+C=" << Euler << '\n';	
}



void runSimulation(Global& global, int time_steps)
{
	for (auto& cell : global.cellMap()) { cell.second.calcm(); }
    for (int step = 0; step < time_steps; step++)
    { 
		for (auto& edge : global.edgeMap()) { edge.second.calcLength(); }
		for (auto& cell : global.cellMap()) 
		{ 
			cell.second.calcL();
			cell.second.calcA(); 
			cell.second.calcT_A();
			cell.second.calcG();
		} 
		WriteVTKFile(global, "graph" + std::to_string(step) + ".vtk", "directors" + std::to_string(step) + ".vtk");
		
		for (auto& edge : global.edgeMap()) { edge.second.calcT_l(); }
		for (auto& vertex : global.vertexMap()) { vertex.second.calcForce(); }		
        for (auto& vertex : global.vertexMap()) { vertex.second.applyForce(); }
        	
    }
    WriteVTKFile(global, "graph" + std::to_string(timesteps) + ".vtk", "directors" + std::to_string(timesteps) + ".vtk");
}



void outputData(Global& global)
{
    for (auto& cell : global.cellMap()) 
    {
        std::cout << "Cell (" << cell.first << ") : ";
        for (int vertex_id : cell.second.getVertices()) { std::cout << global.vertexMap().at(vertex_id).getR() << ", "; }
        std::cout << '\n';
        for (int edge_id : cell.second.getEdges()) { std::cout << global.edgeMap().at(edge_id).getE().first << ' ' << global.edgeMap().at(edge_id).getE().second << '\n'; }
        std::cout << '\n';
    }
}



void WriteVTKFile(Global& global, const std::string& filename_graph, const std::string& filename_director) {

    std::unordered_map<int, int> index_map;
    int index = 0;
    for (const auto& vertex : global.vertexMap()) { index_map[vertex.first] = index++; }

    std::ofstream graphFile(filename_graph);
    graphFile << "# vtk DataFile Version 2.0\nGraph\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << global.vertexMap().size() << " float\n";
    for (const auto& vertex : global.vertexMap()) { graphFile << vertex.second.getR().x() << " " << vertex.second.getR().y() << " 0\n"; }
    
    int n = 0; for (const auto& cell : global.cellMap()) { n += cell.second.getVertices().size(); }
    n += global.cellMap().size();
    graphFile << "CELLS " << global.cellMap().size() << " " << n << '\n';
    for (const auto& cell : global.cellMap()) 
    {
		graphFile << cell.second.getVertices().size() << " ";
		for (int v : cell.second.getVertices()) { graphFile << index_map[v] << " "; }
		graphFile << '\n';
	}
	
	graphFile << "CELL_TYPES " << global.cellMap().size() << '\n';
	for (int c = 0; c < global.cellMap().size(); c++) { graphFile << "7\n"; }
	   
    graphFile.close();
    
    std::ofstream directorFile(filename_director);
    directorFile << "# vtk DataFile Version 2.0\nDirectors\nASCII\nDATASET POLYDATA\nPOINTS " << 2*global.cellMap().size() << " float\n";
    for (const auto& cell : global.cellMap()) 
    { 
		directorFile << (cell.second.getCentroid()-cell.second.getDirector()).x() << " " << (cell.second.getCentroid()-cell.second.getDirector()).y() << " 0\n";
		directorFile << (cell.second.getCentroid()+cell.second.getDirector()).x() << " " << (cell.second.getCentroid()+cell.second.getDirector()).y() << " 0\n";
	}
	
	directorFile << "LINES " << global.cellMap().size() << " " << 3*global.cellMap().size() << "\n";
	for (int i = 0; i < global.cellMap().size(); i++) { directorFile << "2 " << 2*i << " " << 2*i+1 << "\n"; }
	
	directorFile.close();
}

