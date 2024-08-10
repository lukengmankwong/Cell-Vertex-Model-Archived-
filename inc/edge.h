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
    std::pair<int, int> e;
    double l; //length
    double T_l; //line tension
    
    std::unordered_set<int> cell_junctions;
  
public:
    Edge(Global* g, int v1, int v2);
    bool operator==(const Edge& other) const;
    
    const int ID() const;
    const std::pair<int, int>& getE() const;
    const std::unordered_set<int>& getCellJunctions() const;
    const double getl() const;
    const double getT_l() const;

    void addCellJunction(int cell_id);
    void removeCellJunction(int cell_id);
    
    bool swapVertex(int v_old, int v_new);
    
    void calcLength();
    void calcT_l();

};

#endif // EDGE_H
