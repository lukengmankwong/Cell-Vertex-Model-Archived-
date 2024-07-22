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

#include <Eigen/Dense>

int main() 
{

    std::vector<Point> points;
    for (int i = 0; i < cell_count; i++) {
        points.push_back(Point((static_cast<double>(std::rand())/RAND_MAX), (static_cast<double>(std::rand())/RAND_MAX)));
    }
    /*for (int i = 0; i < cell_count; i++) {
		for (int j = 0; j < cell_count; j++) {
			points.push_back(Point(i,j));
		}
	}*/
	
    DT d_tri; 
    d_tri.insert(points.begin(), points.end()); //Delauney triangulation from points
    VD vd(d_tri); //Voronoi diagram dual to Delauney triangulation

    std::unordered_map<int, Vertex> vertex_map; int vertex_counter = 0;
    std::unordered_map<int, Edge> edge_map; int edge_counter = 0;
    std::unordered_map<int, Cell> cell_map; int cell_counter = 0;

	auto t_start1 = std::chrono::high_resolution_clock::now();
    getInitialData(vd, vertex_map, edge_map, cell_map, vertex_counter, edge_counter, cell_counter);
    auto t_end1 = std::chrono::high_resolution_clock::now();
    std::cout << "DATA COLLECTED IN " << std::chrono::duration<double, std::milli>(t_end1 - t_start1).count()/1000 << "s\n";
    
    //outputData(vertex_map, edge_map, cell_map);
	//for (auto edge : edge_map) { std::cout << edge.second.getID() << ' ' << edge.first << '\n';}
	//for (auto vertex : vertex_map) {std::cout << vertex.second.getID() << ' ' << vertex.first << '\n'; }
    int Euler = vertex_map.size() - edge_map.size() + cell_map.size();
    std::cout << "\nV=" << vertex_map.size() << "\nE=" << edge_map.size() << "\nF=" << cell_map.size() << "\nV-E+F=" <<  Euler << '\n';
    WriteVTKFile(vertex_map, edge_map, "Initial.vtk");
    std::cout << "\nPRESS ENTER TO RUN SIMULATION"; std::cin.get();

    auto t_start2 = std::chrono::high_resolution_clock::now();
    runSimulation(vertex_map, edge_map, cell_map, timesteps);
    auto t_end2 = std::chrono::high_resolution_clock::now();

    WriteVTKFile(vertex_map, edge_map, "Final.vtk");


    std::cout << "SIMULATION RAN IN " << std::chrono::duration<double, std::milli>(t_end2 - t_start2).count()/1000 << "s\n";
    return 0;
}


