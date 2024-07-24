#ifndef GLOBALS_H
#define GLOBALS_H

#include <utility>
#include <unordered_map>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_2.h>
#include <CGAL/Vector_2.h>
typedef CGAL::Exact_predicates_inexact_constructions_kernel                   K;
typedef CGAL::Vector_2<K>                                                   Vec;
typedef CGAL::Point_2<K>												  Point;

#include "vertex.h"
#include "edge.h"
#include "cell.h"

class Vertex;
class Edge;
class Cell;

extern std::unordered_map<int, Vertex> vertex_map; extern int vertex_counter;
extern std::unordered_map<int, Edge> edge_map; extern int edge_counter;
extern std::unordered_map<int, Cell> cell_map; extern int cell_counter;


//parameters
extern const int cell_count;

extern const double dt;
extern const int timesteps;

extern const double A_0;

extern const double a;
extern const double k;

#endif // GLOBALS_H
