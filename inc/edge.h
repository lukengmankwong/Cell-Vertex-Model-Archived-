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
    Vertex* v_1; Vertex* v_2;
    double l_; //length
    double T_l_; //line tension
    
    std::unordered_set<Cell*> cell_junctions_;
  
public:

    Edge(Tissue* T, Vertex* v_1, Vertex* v_2);
    Edge();
    bool operator==(const Edge& other) const;

    Vertex* const v1() const; Vertex* const v2() const;
    const double l() const;
    const double T_l() const;
    const std::unordered_set<Cell*>& cellJunctions() const;

    void addCellJunction(Cell* c);
    void removeCellJunction(Cell* c);
    
    const bool hasVertex(Vertex* v) const;
    bool swapVertex(Vertex* v_old, Vertex* v_new); //changes cell vertices
    bool swapVertex_noedit(Vertex* v_old, Vertex* v_new);
    
    void calcLength();
    void calcT_l();
    
    void T1();

};

#endif // EDGE_H
