#include <iostream>

#include <cmath>
#include <vector>
#include <unordered_map>

#include <ctime>
#include <chrono>

#include "functions.h"
#include "global.h"

#include "vertex.h"
#include "edge.h"
#include "cell.h"

bool doughnut(const Point& p) { return p.x()*p.x() + p.y()*p.y() > 0.03 && p.x()*p.x() + p.y()*p.y() < 0.25; }
bool square(const Point& p) { return std::fabs(p.x()) < 0.5 && std::fabs(p.y()) < 0.5; }
bool unbound(const Point& p) { return true; }

int main() 
{
	Global& global = Global::get();
	
    std::vector<Point> points; points.reserve(cell_count);
    //for (int i = 0; i < cell_count; i++) { points.push_back( Point( ((static_cast<double>(std::rand())/RAND_MAX)-0.5), ((static_cast<double>(std::rand())/RAND_MAX)-0.5) ) ); }	
	int k = std::sqrt(cell_count);
	for (int i = -k/2; i < k/2; i++)
	{
		for (int j = -k/2; j < k/2; j++)
		{
			Point p((i+((static_cast<double>(std::rand())/RAND_MAX)-0.5))/k, (j+0.5*(i%2)+((static_cast<double>(std::rand())/RAND_MAX)-0.5))/k);
			points.push_back(p);
		}
	}
	
    DT d_tri; 
    d_tri.insert(points.begin(), points.end()); //Delauney triangulation from points
    VD vd(d_tri); //Voronoi diagram dual to Delauney triangulation


	auto t_start1 = std::chrono::high_resolution_clock::now();
    getInitialData(vd, global, square);
    auto t_end1 = std::chrono::high_resolution_clock::now();
    //outputData(global);
    std::cout << "DATA COLLECTED IN " << std::chrono::duration<double, std::milli>(t_end1 - t_start1).count()/1000 << "s\n";

    std::cout << "\nPRESS ENTER TO RUN SIMULATION"; std::cin.get();
    auto t_start2 = std::chrono::high_resolution_clock::now();
    runSimulation(global, timesteps);
    auto t_end2 = std::chrono::high_resolution_clock::now();

    std::cout << "SIMULATION RAN IN " << std::chrono::duration<double, std::milli>(t_end2 - t_start2).count()/1000 << "s\n";
    return 0;
}

