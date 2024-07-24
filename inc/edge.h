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
    
    int cell_junction_count;
    double length;
    
public:
    Edge(int id, int v1, int v2);
    bool operator==(const Edge& other) const;
    
    const int getID() const;
    const std::pair<int, int>& getE() const;
    int getCellJunctions() const;

    void addCellJunction();
    void removeCellJunction();
    
    bool swapVertex(int v_old, int v_new);
    
    void calcLength();

};

#endif // EDGE_H
