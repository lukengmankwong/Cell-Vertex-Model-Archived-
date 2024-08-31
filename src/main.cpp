#include <iostream>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <ctime>
#include <chrono>


#include "tissue.h"
#include "vertex.h"
#include "edge.h"
#include "cell.h"

#include "functions.h"
#include "parameters.h"

bool doughnut(const Point& p) { return p.x()*p.x() + p.y()*p.y() > 0.03 && p.x()*p.x() + p.y()*p.y() < 0.25; }
bool square(const Point& p) { return std::fabs(p.x()) < 15 && std::fabs(p.y()) < 15; }
bool unbound(const Point& p) { return true; }
bool circle(const Point& p) { return (p.x()-0)*(p.x()-0) + (p.y()-0)*(p.y()-0) < 700; }

int main() 
{
	
	Tissue tissue = Tissue();
	int cell_count = 2800;
	
    std::vector<Point> points; points.reserve(cell_count);
    //for (int i = 0; i < 1000; i++) { points.push_back( Point( 30*((static_cast<double>(std::rand())/RAND_MAX)-0.5), 30*((static_cast<double>(std::rand())/RAND_MAX)-0.5) ) ); }	
	int k = std::sqrt(cell_count);
	for (int i = -k/2; i < k/2; i++)
	{
		for (int j = -k/2; j < k/2; j++)
		{
			Point p(i+0.15*((static_cast<double>(std::rand())/RAND_MAX)-0.5), j+0.5*(i%2)+0.15*((static_cast<double>(std::rand())/RAND_MAX)-0.5));
			//Point p(i, j+0.5*(i%2));
			points.push_back(p);
		}
	}
	
	/*for (int i = 0; i*i < cell_count; i++)
	{
		for (int j = 0; j*j < cell_count; j++)
		{
			Point p(i, j);
			points.push_back(p);
		}
	}*/
	
    DT d_tri; 
    d_tri.insert(points.begin(), points.end()); //Delauney triangulation from points
    VD vd(d_tri); //Voronoi diagram dual to Delauney triangulation


	auto t_start1 = std::chrono::high_resolution_clock::now();
    getInitialData(vd, tissue, circle);
    auto t_end1 = std::chrono::high_resolution_clock::now();
    //outputData(Tissue);
    std::cout << "DATA COLLECTED IN " << std::chrono::duration<double, std::milli>(t_end1 - t_start1).count()/1000 << "s\n";

    std::cout << "\nPRESS ENTER TO RUN SIMULATION"; std::cin.get();
    auto t_start2 = std::chrono::high_resolution_clock::now();
    tissue.run(30000);
    auto t_end2 = std::chrono::high_resolution_clock::now();

    std::cout << "SIMULATION RAN IN " << std::chrono::duration<double, std::milli>(t_end2 - t_start2).count()/1000 << "s\n";
    return 0;
}

