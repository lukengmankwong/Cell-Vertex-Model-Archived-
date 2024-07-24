#include <iostream>

#include <cmath>
#include <vector>
#include <unordered_map>

#include <ctime>
#include <chrono>

#include "functions.h"
#include "globals.h"

#include "vertex.h"
#include "edge.h"
#include "cell.h"

bool in(const Point& p)
{
	return p.x()*p.x() + p.y()*p.y() > 0.075 && p.x()*p.x() + p.y()*p.y() < 0.25;
}

int main() 
{

    std::vector<Point> points;
    for (int i = 0; i < cell_count; i++) {
        points.push_back( Point( ((static_cast<double>(std::rand())/RAND_MAX)-0.5), ((static_cast<double>(std::rand())/RAND_MAX)-0.5) ) );
    }
	
    DT d_tri; 
    d_tri.insert(points.begin(), points.end()); //Delauney triangulation from points
    VD vd(d_tri); //Voronoi diagram dual to Delauney triangulation


	auto t_start1 = std::chrono::high_resolution_clock::now();
    getInitialData(vd, in);
    auto t_end1 = std::chrono::high_resolution_clock::now();
    std::cout << "DATA COLLECTED IN " << std::chrono::duration<double, std::milli>(t_end1 - t_start1).count()/1000 << "s\n";
    

    int Euler = vertex_map.size() - edge_map.size() + cell_map.size();
    std::cout << "\nV=" << vertex_map.size() << "\nE=" << edge_map.size() << "\nF=" << cell_map.size() << "\nV-E+F=" <<  Euler << '\n';
    WriteVTKFile(vertex_map, edge_map, cell_map, "Initial.vtk");
    std::cout << "\nPRESS ENTER TO RUN SIMULATION"; std::cin.get();

    auto t_start2 = std::chrono::high_resolution_clock::now();
    runSimulation(timesteps);
    auto t_end2 = std::chrono::high_resolution_clock::now();

    WriteVTKFile(vertex_map, edge_map, cell_map, "Final.vtk");


    std::cout << "SIMULATION RAN IN " << std::chrono::duration<double, std::milli>(t_end2 - t_start2).count()/1000 << "s\n";
    return 0;
}


