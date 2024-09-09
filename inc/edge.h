#ifndef EDGE_H
#define EDGE_H

#include <cmath>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <array>
#include <iostream>

#include "parameters.h"
#include "libraries.h"
#include "vertex.h"

class Tissue;

class Cell;


class Edge
{
private:
	
	Tissue* T;
    int id;
    int v_1, v_2;
    Vertex* v_i, v_j;
    double l_; //length
    double T_l_; //line tension
    
    std::unordered_set<int> cell_junctions;
    std::unordered_set<Cell*> cell_junctions_;
  
public:

    Edge(Tissue* T, int v1, int v2);
    Edge();
    bool operator==(const Edge& other) const;

    const int v1() const; const int v2() const;
    const double l() const;
    const double T_l() const;
    const std::unordered_set<int>& cellJunctions() const;

    void addCellJunction(int cell_id);
    void removeCellJunction(int cell_id);
    
    const bool hasVertex(int v) const;
    bool swapVertex(int v_old, int v_new); //changes cell vertices
    bool swapVertex_noedit(int v_old, int v_new);
    
    void calcLength();
    void calcT_l();
    
    void T1();

};

#endif // EDGE_H
