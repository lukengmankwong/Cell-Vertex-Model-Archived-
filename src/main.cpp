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


//some example indicator functions for the initial shape
bool doughnut(const Point& p) 	{ return p.x()*p.x() + p.y()*p.y() > 100 && p.x()*p.x() + p.y()*p.y() < 700; }
bool square(const Point& p)		{ return std::fabs(p.x()) < std::sqrt(1000)/2 && std::fabs(p.y()) < std::sqrt(1000)/2; }
bool unbound(const Point& p) 	{ return true; }
bool circle(const Point& p) 	{ return (p.x()-0)*(p.x()-0) + (p.y()-0)*(p.y()-0) < 700; }

//some example voronoi seeds
std::vector<Point> randomPoints(unsigned int n) 
{ 
	std::vector<Point> points;
	double w = std::sqrt(n);
	for (int i = 0; i < n; i++) points.push_back(Point(w*((static_cast<double>(std::rand())/RAND_MAX)-0.5), w*((static_cast<double>(std::rand())/RAND_MAX)-0.5)));
	return points;
}
std::vector<Point> hexagonalWithNoise(unsigned int n, double g)
{
	std::vector<Point> points;
	double k = std::sqrt(n);
	for (int i = -k/2; i < k/2; i++)
	{
		for (int j = -k/2; j < k/2; j++)
		{
			Point p(i+g*((static_cast<double>(std::rand())/RAND_MAX)-0.5), j+0.5*(i%2)+g*((static_cast<double>(std::rand())/RAND_MAX)-0.5));
			points.push_back(p);
		}
	}
	return points;
}


int main() 
{
	unsigned int cell_count = 2800;
    std::vector<Point> points = hexagonalWithNoise(cell_count, 0.1);
    //std::vector<Point> points = randomPoints(cell_count);
    DT delauney_tri; 
    delauney_tri.insert(points.begin(), points.end()); 		//Delauney triangulation from points
    VD voronoi_diagram(delauney_tri); 						//Voronoi diagram dual to Delauney triangulation

	/*auto t_start1 = std::chrono::high_resolution_clock::now();
	Tissue T = Tissue(voronoi_diagram, circle);
    auto t_end1 = std::chrono::high_resolution_clock::now();
    std::cout << "DATA COLLECTED IN " << std::chrono::duration<double, std::milli>(t_end1 - t_start1).count()/1000 << "s\n";*/
    
    unsigned int timesteps = 100000;
    param::set_GAMMA(0.2);
    param::set_LAMBDA(-0.2);
    int i = 0;
    for (double LAMBDA = -0.5; LAMBDA < 0.21; LAMBDA += 0.1)
    {
		Tissue T = Tissue(voronoi_diagram, circle);
		param::set_LAMBDA(LAMBDA);
		T.run(timesteps, std::to_string(i));
		i++;
		std::cout << LAMBDA << " run success\n";
	}
    /*std::cout << "\nPRESS ENTER TO RUN SIMULATION"; std::cin.get();
    auto t_start2 = std::chrono::high_resolution_clock::now();
    T.run(timesteps, "");
    auto t_end2 = std::chrono::high_resolution_clock::now();
    std::cout << "SIMULATION RAN IN " << std::chrono::duration<double, std::milli>(t_end2 - t_start2).count()/1000 << "s\n";*/
    return 0;
}
