#include "functions.h"


void removeDuplicates(std::vector<int>& vec) {
    std::unordered_set<int> seen;   // To track seen elements
    auto it = vec.begin();

    while (it != vec.end()) {
        // If the element is seen for the first time, keep it
        if (seen.insert(*it).second) {
            ++it;
        } else {
            // Otherwise, erase the duplicate
            it = vec.erase(it);
        }
    }
}


void getInitialData(VD& vd, Tissue& Tissue, bool (*in)(const Point&))
{
    std::cout << "COLLECTING INITIAL DATA\n";
	for (VD::Vertex_iterator vit = vd.vertices_begin(); vit != vd.vertices_end(); vit++)
	{       
        Tissue.createVertex(vit->point());
    }
    
    for (VD::Face_iterator fi = vd.faces_begin(); fi != vd.faces_end(); fi++) 
    {
        std::vector<int> cell_vertices;
        VD::Ccb_halfedge_circulator ec_start = fi->ccb();
        VD::Ccb_halfedge_circulator ec = ec_start;

        do { 
			if (!ec->is_unbounded()) 
			{
				for (const auto& vertex : Tissue.vertexMap())
				{
					if (ec->source()->point() == vertex.second.r())
					{
						cell_vertices.push_back(vertex.first);
						break;
					}
				}
			}
			++ec;
        } while (ec != ec_start);
		
        std::vector<int> cell_edges; 
        if (true)
        {
			for (int i = 0; i < cell_vertices.size(); i++)
			{
				int v1 = cell_vertices[i]; int v2 = cell_vertices[(i+1)%cell_vertices.size()];
				bool found = false;
				for (const auto& edge : Tissue.edgeMap())
				{
					if (edge.second.v1() == v1 && edge.second.v2() == v2)
					{
						cell_edges.push_back(edge.first);
						found = true;
						break;
					}
					else if (edge.second.v1() == v2 && edge.second.v2() == v1)
					{
						cell_edges.push_back(edge.first);
						found = true;
						break;
					}
				}
				if (!found)
				{
					cell_edges.push_back(Tissue.createEdge(v1, v2));
				}		
			}
		}
		//for (int v : cell_edges) { std::cout << v << '\n'; } std::cout << '\n';
		removeDuplicates(cell_edges); //temporary fix
        Tissue.createCell(cell_vertices, cell_edges);
    } 
    
    
    
	std::unordered_set<int> cells_to_remove;
	for (const auto& vertex : Tissue.vertexMap())  
	{ 
		if ( !in(vertex.second.r()) ) 
		{ 
			cells_to_remove.insert(vertex.second.cellContacts().begin(), vertex.second.cellContacts().end()); 
		}
	}
	for (const auto& cell : Tissue.cellMap())
	{
		//if (cell.second.Edges().size() < 3) { cells_to_remove.insert(cell.first); }
	}
	for (int c : cells_to_remove) { Tissue.destroyCell(c); }
	int V = Tissue.vertexMap().size(); int E = Tissue.edgeMap().size(); int C = Tissue.cellMap().size();
	int Euler = V-E+C;
    std::cout << "V=" << V << "\nE=" << E << "\nC=" << C << "\nV-E+C=" << Euler << '\n';	
	Tissue.cellFindNeighbours();
}


void outputData(const Tissue& Tissue)
{
    for (auto& cell : Tissue.cellMap()) 
    {
        std::cout << "Cell (" << cell.first << ") : ";
        for (int vertex_id : cell.second.Vertices()) { std::cout << Tissue.vertexMap().at(vertex_id).r() << ", "; }
        std::cout << '\n';
        for (int edge_id : cell.second.Edges()) { std::cout << Tissue.edgeMap().at(edge_id).v1() << ' ' << Tissue.edgeMap().at(edge_id).v2() << '\n'; }
        std::cout << '\n';
    }
}

void writeCellsFile(Tissue* T, const std::string& filename_cells)
{
	std::unordered_map<int, int> index_map;
    int index = 0;
    for (const auto& vertex : T->vertexMap()) { index_map[vertex.first] = index++; }

    std::ofstream graphFile(filename_cells);
    graphFile << "# vtk DataFile Version 2.0\nGraph\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << T->vertexMap().size() << " float\n";
    for (const auto& vertex : T->vertexMap()) { graphFile << vertex.second.r().x() << " " << vertex.second.r().y() << " 0\n"; }
    
    int n = 0; for (const auto& cell : T->cellMap()) { n += cell.second.Vertices().size(); }
    n += T->cellMap().size();
    graphFile << "CELLS " << T->cellMap().size() << " " << n << '\n';
    for (const auto& cell : T->cellMap()) 
    {
		graphFile << cell.second.Vertices().size() << " ";
		for (int v : cell.second.Vertices()) { graphFile << index_map[v] << " "; }
		graphFile << '\n';
	}
	
	graphFile << "CELL_TYPES " << T->cellMap().size() << '\n';
	for (int c = 0; c < T->cellMap().size(); c++) { graphFile << "7\n"; }
	   
    graphFile.close();
}

void writeDirectorsFile(Tissue* T, const std::string& filename_directors)
{
	std::ofstream directorFile(filename_directors);
    directorFile << "# vtk DataFile Version 2.0\nns\nASCII\nDATASET POLYDATA\nPOINTS " << 2*T->cellMap().size() << " float\n";
    for (const auto& cell : T->cellMap()) 
    { 
		directorFile << (cell.second.r_0()-0.5*cell.second.n()).x() << " " << (cell.second.r_0()-0.5*cell.second.n()).y() << " 0\n";
		directorFile << (cell.second.r_0()+0.5*cell.second.n()).x() << " " << (cell.second.r_0()+0.5*cell.second.n()).y() << " 0\n";
	}
	
	directorFile << "LINES " << T->cellMap().size() << " " << 3*T->cellMap().size() << "\n";
	for (int i = 0; i < T->cellMap().size(); i++) { directorFile << "2 " << 2*i << " " << 2*i+1 << "\n"; }
	
	directorFile.close();
}

void writeCellDefectsFile(Tissue* T, const std::string& filename_cell_defects)
{
	std::unordered_set<int> defect_vertices;
	const std::vector<int>& cell_defects = T->cellStepDefects();
	for (int c : cell_defects) { defect_vertices.insert(T->cell(c).Vertices().begin(), T->cell(c).Vertices().end()); }
			
	std::unordered_map<int, int> index_map;
	int index = 0;
	for (int v : defect_vertices) { index_map[v] = index++; }
		
	std::ofstream defectFile(filename_cell_defects);
	defectFile << "# vtk DataFile Version 2.0\nDefect\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << defect_vertices.size() << " float\n";
	for (int v : defect_vertices) { defectFile << T->vert(v).r().x() << " " << T->vert(v).r().y() << " 0\n"; }

	int n = 0; 
	for (int c : cell_defects) { n += T->cell(c).Vertices().size(); }
	n += cell_defects.size();

	defectFile << "CELLS " << cell_defects.size() << " " << n << '\n';
	for (int c : cell_defects) 
	{
		defectFile << T->cell(c).Vertices().size() << " ";
		for (int v : T->cell(c).Vertices()) { defectFile << index_map[v] << " "; }
		defectFile << '\n'; 
	}
		
	defectFile << "CELL_TYPES " << cell_defects.size() << '\n';
	for (int c = 0; c < cell_defects.size(); c++) { defectFile << "7\n"; }
		   
	defectFile.close();
}

void writeVertexDefectsFile(Tissue* T, const std::string& filename_vertex_defects)
{
	std::ofstream file(filename_vertex_defects);
    file << "# vtk DataFile Version 3.0\nPoint data\nASCII\nDATASET POLYDATA\n";
    file << "POINTS " << T->vertexStepDefects().size() << " float\n";
    for (int v : T->vertexStepDefects()) file << T->vert(v).r().x() << " " << T->vert(v).r().y() << " 0\n";
    file.close();
}



