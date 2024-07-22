#ifndef EDGE_H
#define EDGE_H

#include <cmath>
#include <utility>
#include <unordered_map>

class Edge
{
private:
    const int id;
    std::pair<int, int> e;
    
    int cell_junction_count;
    double length;
    
    std::unordered_map<int, Edge>* edge_map;
    
public:
    Edge(std::unordered_map<int, Edge>* edge_map, int id, int v1, int v2);
    bool operator==(const Edge& other) const;

    void addCellJunction();
    void removeCellJunction();
    int getCellJunctions() const;
    
    const int getID() const;

    const std::pair<int, int>& getE() const;
    
    bool swapVertex(int v_old, int v_new);

};

#endif // EDGE_H
