#include <iostream>

#include <vector>
#include <unordered_map>
#include <cmath>

#include <ctime>
#include <chrono>

#include "functions.h"
#include "parameters.h"

#include "vertex.h"
#include "edge.h"
#include "cell.h"


int main() 
{

    std::vector<Point> points;
    for (int i = 0; i < cell_count; i++) {
        points.push_back(Point((static_cast<double>(std::rand())/RAND_MAX), (static_cast<double>(std::rand())/RAND_MAX)));
    }

    DT d_tri; 
    d_tri.insert(points.begin(), points.end()); //Delauney triangulation from points
    VD vd(d_tri); //Voronoi diagram dual to Delauney triangulation

    std::unordered_map<int, Vertex> vertex_map; int vertex_counter = 0;
    std::unordered_map<int, Edge> edge_map; int edge_counter = 0;
    std::unordered_map<int, Cell> cell_map; int cell_counter = 0;


    getInitialData(vd, cell_map, edge_map, vertex_map, vertex_counter, edge_counter, cell_counter);
    for (auto& cell : cell_map) {
        for (int v : cell.second.getVertices()) {
            std::cout << v << ' ';
        } std::cout << '\n';
    }
    outputData(vertex_map, cell_map);
    std::cout << "\nV=" << vertex_map.size() << "\nE=" << edge_map.size() << "\nF=" << cell_map.size() << "\nV-E+F=" << vertex_map.size() - edge_map.size() + cell_map.size() << '\n';
    WriteVTKFile(vertex_map, edge_map, "Initial.vtk");
    std::cout << "\nINITIAL DATA COLLECTED\nPRESS ENTER TO RUN SIMULATION"; std::cin.get();

    auto t_start = std::chrono::high_resolution_clock::now();
    runSimulation(cell_map, edge_map, vertex_map, timesteps);
    auto t_end = std::chrono::high_resolution_clock::now();

    WriteVTKFile(vertex_map, edge_map, "Final.vtk");


    std::cout << "SIMULATION RAN IN " << std::chrono::duration<double, std::milli>(t_end - t_start).count()/1000 << "s\n";
    return 0;
}


