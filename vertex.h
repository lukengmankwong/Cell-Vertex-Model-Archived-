#ifndef VERTEX_H
#define VERTEX_H

#include "libraries.h"
#include "parameters.h"

class Vertex
{
private:
    Point r;
    Vec force;
    
public:
    Vertex(Point r);
    bool operator==(const Vertex& other) const;

    Point getR() const;

    void calcForce(Point centroid);
    void applyForce();
    
};

#endif // VERTEX_H