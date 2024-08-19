#ifndef EDGE_H
#define EDGE_H

#include <cmath>
#include <utility>

#include "global.h"
class Global;

class Edge
{
private:
	
	Global* g;
    const int id;
    int v_1, v_2;
    double l_; //length
    double T_l_; //line tension
    
    std::unordered_set<int> cell_junctions;
  
public:
    Edge(Global* g, int v1, int v2);
    bool operator==(const Edge& other) const;

    const int v1() const; const int v2() const;
    const std::unordered_set<int>& cellJunctions() const;
    const double l() const;
    const double T_l() const;

    void addCellJunction(int cell_id);
    void removeCellJunction(int cell_id);
    
    const bool hasVertex(int v) const;
    bool swapVertex_rep(int v_old, int v_new); //changes cell vertices
    
    void calcLength();
    void calcT_l();
    
    void T1merge();

};

#endif // EDGE_H
