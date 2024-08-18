#include "functions.h"

void removeDuplicates(std::vector<std::pair<int, int>>& vec) {
    std::unordered_set<int> seen;   // To track seen elements
    auto it = vec.begin();

    while (it != vec.end()) {
        // If the element is seen for the first time, keep it
        if (seen.insert(it->first).second) {
            ++it;
        } else {
            // Otherwise, erase the duplicate
            it = vec.erase(it);
        }
    }
}


void getInitialData(VD& vd, Global& global, bool (*in)(const Point&))
{
    std::cout << "COLLECTING INITIAL DATA\n";
	for (VD::Vertex_iterator vit = vd.vertices_begin(); vit != vd.vertices_end(); vit++)
	{       
        global.createVertex(vit->point());
    }
    
    for (VD::Face_iterator fi = vd.faces_begin(); fi != vd.faces_end(); fi++) 
    {
        std::vector<int> cell_vertices;
        VD::Ccb_halfedge_circulator ec_start = fi->ccb();
        VD::Ccb_halfedge_circulator ec = ec_start;

        do { 
			if (!ec->is_unbounded()) 
			{
				for (const auto& vertex : global.vertexMap())
				{
					if (ec->source()->point() == vertex.second.R())
					{
						cell_vertices.push_back(vertex.first);
						break;
					}
				}
			}
			++ec;
        } while (ec != ec_start);
		
        std::vector<std::pair<int, int>> cell_edges; 
        if (true)
        {
			for (int i = 0; i < cell_vertices.size(); i++)
			{
				int v1 = cell_vertices[i]; int v2 = cell_vertices[(i+1)%cell_vertices.size()];
				bool found = false;
				for (const auto& edge : global.edgeMap())
				{
					if (edge.second.v1() == v1 && edge.second.v2() == v2)
					{
						cell_edges.push_back(std::make_pair(edge.first, 0));
						found = true;
						break;
					}
					else if (edge.second.v1() == v2 && edge.second.v2() == v1)
					{
						cell_edges.push_back(std::make_pair(edge.first, 1));
						found = true;
						break;
					}
				}
				if (!found)
				{
					cell_edges.push_back(std::make_pair(global.createEdge(v1, v2), 0));
				}		
			}
		}
		//for (int v : cell_edges) { std::cout << v << '\n'; } std::cout << '\n';
		removeDuplicates(cell_edges); //temporary fix
        global.createCell(cell_vertices, cell_edges);
        //for (int v : cell_vertices) { std::cout << v << ' '; } std::cout << '\n';
		//for (auto e : cell_edges) { std::cout << global.edge(e.first).v1() << ' ' << global.edge(e.first).v2() << '\t'; } std::cout << "\n\n";
    } 
    
    
    
	std::unordered_set<int> cells_to_remove;
	for (const auto& vertex : global.vertexMap())  
	{ 
		if ( !in(vertex.second.R()) ) 
		{ 
			cells_to_remove.insert(vertex.second.cellContacts().begin(), vertex.second.cellContacts().end()); 
		}
	}
	for (const auto& cell : global.cellMap())
	{
		if (cell.second.Edges().size() < 3) { cells_to_remove.insert(cell.first); }
	}
	for (int c : cells_to_remove) { global.destroyCell(c); }
	int V = global.vertexMap().size(); int E = global.edgeMap().size(); int C = global.cellMap().size();
	int Euler = V-E+C;
    std::cout << "V=" << V << "\nE=" << E << "\nC=" << C << "\nV-E+C=" << Euler << '\n';	
}



void runSimulation(Global& global, int time_steps)
{
    for (int step = 0; step < time_steps; step++)
    { 
		//std::cout << "step " << step << ":\n";			
		//std::cout << "cell map size: " << global.cellMap().size() << "\n";
		//std::cout << "edge map size: " << global.edgeMap().size() << "\n";
		//std::cout << "vertex map size: " << global.vertexMap().size() << "\n\n";
		global.transitions();
		
		for (auto& edge : global.edgeMap()) { edge.second.calcLength(); }
		for (auto& cell : global.cellMap()) 
		{ 
			cell.second.calcL();
			cell.second.calcA(); 
			cell.second.calcT_A();
			cell.second.calcG();
		}
		for (auto& cell : global.cellMap()) { cell.second.calcm(); }
		global.addDefects();
		
		for (auto& edge : global.edgeMap()) { edge.second.calcT_l(); }		
		for (auto& vertex : global.vertexMap()) { vertex.second.calcForce(); }		
        for (auto& vertex : global.vertexMap()) { vertex.second.applyForce(); }
        WriteVTKFile(global, "graph" + std::to_string(step) + ".vtk", "defects" + std::to_string(step) + ".vtk", "directors" + std::to_string(step) + ".vtk");
		global.nextStep();
        	
    }
}



void outputData(Global& global)
{
    for (auto& cell : global.cellMap()) 
    {
        std::cout << "Cell (" << cell.first << ") : ";
        for (int vertex_id : cell.second.Vertices()) { std::cout << global.vertexMap().at(vertex_id).R() << ", "; }
        std::cout << '\n';
        for (int edge_id : cell.second.Edges()) { std::cout << global.edgeMap().at(edge_id).v1() << ' ' << global.edgeMap().at(edge_id).v2() << '\n'; }
        std::cout << '\n';
    }
}



void WriteVTKFile(Global& global, const std::string& filename_graph, const std::string& filename_defect, const std::string& filename_director)
{
	//graph
    std::unordered_map<int, int> index_map;
    int index = 0;
    for (const auto& vertex : global.vertexMap()) { index_map[vertex.first] = index++; }

    std::ofstream graphFile(filename_graph);
    graphFile << "# vtk DataFile Version 2.0\nGraph\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << global.vertexMap().size() << " float\n";
    for (const auto& vertex : global.vertexMap()) { graphFile << vertex.second.R().x() << " " << vertex.second.R().y() << " 0\n"; }
    
    int n = 0; for (const auto& cell : global.cellMap()) { n += cell.second.Vertices().size(); }
    n += global.cellMap().size();
    graphFile << "CELLS " << global.cellMap().size() << " " << n << '\n';
    for (const auto& cell : global.cellMap()) 
    {
		graphFile << cell.second.Vertices().size() << " ";
		for (int v : cell.second.Vertices()) { graphFile << index_map[v] << " "; }
		graphFile << '\n';
	}
	
	graphFile << "CELL_TYPES " << global.cellMap().size() << '\n';
	for (int c = 0; c < global.cellMap().size(); c++) { graphFile << "7\n"; }
	   
    graphFile.close();
 

    //defects
	std::unordered_set<int> defect_vertices;
	for (int c : global.stepDefects(global.Step())) { defect_vertices.insert(global.cell(c).Vertices().begin(), global.cell(c).Vertices().end()); }
			
	std::unordered_map<int, int> index_map2;
	int index2 = 0;
	for (int v : defect_vertices) { index_map2[v] = index2++; }
		
	std::ofstream defectFile(filename_defect);
	defectFile << "# vtk DataFile Version 2.0\nDefect\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << defect_vertices.size() << " float\n";
	for (int v : defect_vertices) { defectFile << global.vert(v).R().x() << " " << global.vert(v).R().y() << " 0\n"; }

	int n2 = 0; 
	for (int c : global.stepDefects(global.Step())) { n2 += global.cell(c).Vertices().size(); }
	n2 += global.stepDefects(global.Step()).size();

	defectFile << "CELLS " << global.stepDefects(global.Step()).size() << " " << n2 << '\n';
	for (int c : global.stepDefects(global.Step())) 
	{
		defectFile << global.cell(c).Vertices().size() << " ";
		for (int v : global.cell(c).Vertices()) { defectFile << index_map2[v] << " "; }
		defectFile << '\n'; 
	}
		
	defectFile << "CELL_TYPES " << global.stepDefects(global.Step()).size() << '\n';
	for (int c = 0; c < global.stepDefects(global.Step()).size(); c++) { defectFile << "7\n"; }
		   
	defectFile.close();

    
    
    //directors
    std::ofstream directorFile(filename_director);
    directorFile << "# vtk DataFile Version 2.0\nns\nASCII\nDATASET POLYDATA\nPOINTS " << 2*global.cellMap().size() << " float\n";
    for (const auto& cell : global.cellMap()) 
    { 
		directorFile << (cell.second.R_0()-0.003*cell.second.N()).x() << " " << (cell.second.R_0()-0.003*cell.second.N()).y() << " 0\n";
		directorFile << (cell.second.R_0()+0.003*cell.second.N()).x() << " " << (cell.second.R_0()+0.003*cell.second.N()).y() << " 0\n";
	}
	
	directorFile << "LINES " << global.cellMap().size() << " " << 3*global.cellMap().size() << "\n";
	for (int i = 0; i < global.cellMap().size(); i++) { directorFile << "2 " << 2*i << " " << 2*i+1 << "\n"; }
	
	directorFile.close();
}

