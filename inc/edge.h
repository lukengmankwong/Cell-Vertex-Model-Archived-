#ifndef EDGE_H
#define EDGE_H

#include <cmath>
#include <utility>

#include "globals.h"


class Edge
{
private:
    const int id;
    std::pair<int, int> e;
    double l; double dl; //length
    double T_l; //line tension
    
    std::unordered_set<int> cell_junctions;
  
public:
    Edge(int id, int v1, int v2);
    bool operator==(const Edge& other) const;
    
    const int getID() const;
    const std::pair<int, int>& getE() const;
    const std::unordered_set<int>& getCellJunctions() const;
    const double getl() const;
    const double getdl() const;
    const double getT_l() const;

    void addCellJunction(int cell_id);
    void removeCellJunction(int cell_id);
    
    bool swapVertex(int v_old, int v_new);
    
    void calcLength();
    void calcT_l();

};

#endif // EDGE_H
