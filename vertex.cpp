#include "vertex.h"


Vertex::Vertex(Point r) : 
    r(r)
{ force = Vec(0,0); }

bool Vertex::operator==(const Vertex& other) const { return r == other.r; }


Point Vertex::getR() const { return r; }


void Vertex::calcForce(Point centroid)
{
    force += (static_cast<double>(std::rand())/RAND_MAX)*k*(centroid-r);
}

void Vertex::applyForce() 
{ 
    r += a*dt*force;
    force = Vec(0,0);
}


