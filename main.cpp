#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>

#include "functions.h"
#include "parameters.h"

#include "vertex.h"
#include "edge.h"
#include "cell.h"

void getInitialData(VD& vd, std::vector<Cell>& cells, std::vector<Edge>& edges, std::vector<Vertex>& all_vertices);
void runSimulation(std::vector<Cell>& cells, std::vector<Edge>& edges, std::vector<Vertex>& vertices, int time_steps);
void outputData(const std::vector<Vertex>& vertices, const std::vector<Cell>& cells);
void WriteVTKFile(const std::vector<Vertex>& vertices, const std::vector<Edge>& edges, const std::string& filename);


int main() 
{

    std::vector<Point> points;
    for (int i = 0; i < std::sqrt(cell_count); i++) {
        for (int j = 0; j < std::sqrt(cell_count); j++) {
            points.push_back(Point(i,j+(0.5*(i%2))));
        }
    }

    DT d_tri; 
    d_tri.insert(points.begin(), points.end()); //Delauney triangulation from points
    VD vd(d_tri); //Voronoi diagram dual to Delauney triangulation

    
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Cell> cells;
    getInitialData(vd, cells, edges, vertices);
    std::cout << "V=" << vertices.size() << "\nE=" << edges.size() << "\nF=" << cells.size() << '\n';
    WriteVTKFile(vertices, edges, "Initial.vtk");
    std::cout << "\nINITIAL DATA COLLECTED\nPRESS ENTER TO RUN SIMULATION"; std::cin.get();

    auto t_start = std::chrono::high_resolution_clock::now();

    runSimulation(cells, edges, vertices, 2000);

    auto t_end = std::chrono::high_resolution_clock::now();

    WriteVTKFile(vertices, edges, "Final.vtk");
    //outputData(vertices, cells);

    std::cout << "SIMULATION RAN IN " << std::chrono::duration<double, std::milli>(t_end - t_start).count()/1000 << "s\n";

    return 0;
}


