#include "functions.h"


double random(double min, double max, unsigned int seed)
{
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> dist(min, max);
    return dist(generator);
}


/*void outputData(const Tissue& T)
{
    for (const auto& cell : T.cellMap()) 
    {
        std::cout << "Cell (" << cell.first << ") : ";
        for (int vertex_id : cell.second.Vertices()) { std::cout << T.vertexMap().at(vertex_id).r() << ", "; }
        std::cout << '\n';
        for (int edge_id : cell.second.Edges()) { std::cout << T.edgeMap().at(edge_id).v1() << ' ' << T.edgeMap().at(edge_id).v2() << '\n'; }
        std::cout << '\n';
    }
}*/

void writeCellsFile(Tissue* T, const std::string& filename_cells)
{
	std::vector<Vertex*> vertices = T->vertices();
	Vertex* v_0 = T->v_0();
	std::unordered_map<int, int> index_map;
    int index = 0;
    for (Vertex* v : vertices) { index_map[v-v_0] = index++; }

    std::ofstream graphFile(filename_cells);
    graphFile << "# vtk DataFile Version 2.0\nGraph\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS " << vertices.size() << " float\n";
    for (Vertex* v : vertices) { graphFile << v->r().x() << " " << v->r().y() << " 0\n"; }
    
    std::vector<Cell*> cells = T->cells();
    Cell* c_0 = T->c_0();
    int n = 0; for (Cell* c : cells) { n += c->vertices().size(); }
    n += cells.size();
    graphFile << "CELLS " << cells.size() << " " << n << '\n';
    for (Cell* c : cells) 
    {
		graphFile << c->vertices().size() << " ";
		for (Vertex* v : c->vertices()) { graphFile << index_map[v-v_0] << " "; }
		graphFile << '\n';
	}
	
	graphFile << "CELL_TYPES " << cells.size() << '\n';
	for (int c = 0; c < cells.size(); c++) { graphFile << "7\n"; }
	
    graphFile.close();
}

void writeDirectorsFile(Tissue* T, const std::string& filename_directors)
{
	std::vector<Cell*> cells = T->cells();
	size_t n = cells.size();
	
	std::ofstream directorFile(filename_directors);
    directorFile << "# vtk DataFile Version 2.0\nns\nASCII\nDATASET POLYDATA\nPOINTS " << 2*n << " float\n";
    for (Cell* c : cells) 
    { 
		directorFile << (c->r_0()-0.5*c->n()).x() << " " << (c->r_0()-0.5*c->n()).y() << " 0\n";
		directorFile << (c->r_0()+0.5*c->n()).x() << " " << (c->r_0()+0.5*c->n()).y() << " 0\n";
	}
	
	directorFile << "LINES " << n << " " << 3*n << "\n";
	for (int i = 0; i < n; i++) { directorFile << "2 " << 2*i << " " << 2*i+1 << "\n"; }
	
	directorFile.close();
}

/*void writeCellDefectsFile(Tissue* T, const std::string& filename_cell_defects)
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
    file << "# vtk DataFile Version 2.0\nPoint data\nASCII\nDATASET POLYDATA\n";
    file << "POINTS " << T->vertexStepDefects().size() << " float\n";
    for (int v : T->vertexStepDefects()) file << T->vert(v).r().x() << " " << T->vert(v).r().y() << " 0\n";
    
    file << "VERTICES " << T->vertexStepDefects().size() << " " << 2 * T->vertexStepDefects().size() << "\n";
    for (int i = 0; i < T->vertexStepDefects().size(); i++) file << "1 " << i << "\n";
        
    file.close();
}*/
