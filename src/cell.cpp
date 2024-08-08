#include "cell.h"


Cell::Cell(Global* g, int id, std::vector<int>& vertex_keys, std::vector<int>& edge_keys) : 
	g(g), id(id), vertex_keys(vertex_keys), edge_keys(edge_keys)
{
	A = 0;
	for (int i = 0; i < vertex_keys.size(); i++) {
		A += g->vertexMap().at(vertex_keys[i]).getR().x()*g->vertexMap().at(vertex_keys[(i+1)%vertex_keys.size()]).getR().y()
			- g->vertexMap().at(vertex_keys[(i+1)%vertex_keys.size()]).getR().x()*g->vertexMap().at(vertex_keys[i]).getR().y();
	} A *= 0.5;
	S = A/std::fabs(A);
}

const int Cell::getID() const { return id; }
const Point& Cell::getCentroid() const { return centroid; }
const Vec& Cell::getDirector() const { return director; }

const double Cell::getA() const { return S*A; }
const double Cell::getS() const { return S; }
const double Cell::getL() const { return L; }
const double Cell::getT_A() const { return T_A; }

const std::vector<int>& Cell::getVertices() const { return vertex_keys; }
const std::vector<int>& Cell::getEdges() const { return edge_keys; }


void Cell::addVertex(int v, int i) { vertex_keys.insert(vertex_keys.begin()+i, v); }
void Cell::addEdge(int e, int i) { edge_keys.insert(vertex_keys.begin()+i, e); }


void Cell::removeVertices()
{
    for (int v : vertex_keys) { g->vertexMap().at(v).removeCellContact(id); }
}
void Cell::removeEdges()
{
    for (int e : edge_keys) { g->edgeMap().at(e).removeCellJunction(id); }
}

const int Cell::longestEdge() const //for cell division along longest edges
{
	double longest_l = 0;
	int e_l;
	for (int e : edge_keys) {
		if (g->edgeMap().at(e).getl() > longest_l) {
			longest_l = g->edgeMap().at(e).getl();
			e_l = e;
		}
	} return e_l;
}


void Cell::divide()
{
	
	
}

void Cell::extrude()
{

}


void Cell::calcCentroid()
{
    K::FT x_sum = 0;
    K::FT y_sum = 0;

    for (int v : vertex_keys) {
        x_sum += g->vertexMap().at(v).getR().x();
        y_sum += g->vertexMap().at(v).getR().y();
    }

    centroid = Point(x_sum/vertex_keys.size(), y_sum/vertex_keys.size());
}

void Cell::calcG()
{
	G[0]=0;	G[1]=0;
			G[2]=0;
	
	this->calcCentroid();
	double x_0 = centroid.x(); double y_0 = centroid.y();
	for (int v : vertex_keys)
	{
		double x_v = g->vertexMap().at(v).getR().x();
		double y_v = g->vertexMap().at(v).getR().y();
		G[0]+=(x_v-x_0)*(x_v-x_0);	G[1]+=(x_v-x_0)*(y_v-y_0);
									G[2]+=(y_v-y_0)*(y_v-y_0);
	}
	
	double f = 1.0/vertex_keys.size();
	G[0]*=f; 	G[1]*=f;
				G[2]*=f;
	
	lambda = 0.5*( G[0]+G[2] + std::sqrt( (G[0]+G[2])*(G[0]+G[2]) - 4*(G[0]*G[2]-G[1]*G[1]) ) );
	director = 1e-2*Vec(1, (lambda-G[0])/G[1]) / std::sqrt( 1 + ((lambda-G[0])/G[1])*((lambda-G[0])/G[1]) );
	
	//calculation of traceless of gyration tensor
	double tr_less = 0.5*(G[0]+G[2]);
	TL[0]=G[0]-tr_less;
	TL[1]=G[1];
	TL[2]=G[2]-tr_less;
}

void Cell::calcA()
{
	A = 0;
	for (int i = 0; i < vertex_keys.size(); i++) {
		A += g->vertexMap().at(vertex_keys[i]).getR().x()*g->vertexMap().at(vertex_keys[(i+1)%vertex_keys.size()]).getR().y()
			- g->vertexMap().at(vertex_keys[(i+1)%vertex_keys.size()]).getR().x()*g->vertexMap().at(vertex_keys[i]).getR().y();
	} A *= 0.5;
}

void Cell::calcL()
{
	L = 0;
	//edge lengths must already be calculated
	for (int e : edge_keys) { L += g->edgeMap().at(e).getl(); }
}

void Cell::calcT_A() { T_A = k_A*(A-A_0); }

void Cell::calcm()
{
	double w = 0;
	/*std::cout << '(' << id << ") ";
	for (int v : vertex_keys)
	{
		for (int c : g->vertexMap().at(v).getCellContacts()) { std::cout << c << ' '; }
		std::cout << '\t';
	}
	std::cout << '\n';*/
}





