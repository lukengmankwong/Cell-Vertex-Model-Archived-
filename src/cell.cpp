#include "cell.h"


Cell::Cell(Global* g, std::vector<int>& vertex_keys, std::vector<int>& edge_keys) : 
	g(g), id(g->cellCounter()), vertex_keys(vertex_keys), edge_keys(edge_keys)
{
	A = 0;
	for (int i = 0; i < vertex_keys.size(); i++) {
		A += g->vertexMap().at(vertex_keys[i]).R().x()*g->vertexMap().at(vertex_keys[(i+1)%vertex_keys.size()]).R().y()
			- g->vertexMap().at(vertex_keys[(i+1)%vertex_keys.size()]).R().x()*g->vertexMap().at(vertex_keys[i]).R().y();
	} A *= 0.5;
	S = A/std::fabs(A);
}

const int Cell::ID() const { return id; }
const std::vector<int>& Cell::getVertices() const { return vertex_keys; }
const std::vector<int>& Cell::getEdges() const { return edge_keys; }

const double Cell::getA() const { return S*A; }
const double Cell::getS() const { return S; }
const double Cell::getL() const { return L; }
const double Cell::getT_A() const { return T_A; }
const Point& Cell::getCentroid() const { return centroid; }
const Vec& Cell::getDirector() const { return director; }


void Cell::addVertex(int v, int i) { vertex_keys.insert(vertex_keys.begin()+i, v); }
void Cell::removeVertex(int v) { vertex_keys.erase(std::find(vertex_keys.begin(), vertex_keys.end(), v)); }
void Cell::addEdge(int e, int i) { edge_keys.insert(vertex_keys.begin()+i, e); }
void Cell::removeEdge(int e) 
{ 
	edge_keys.erase(std::find(edge_keys.begin(), edge_keys.end(), e)); 
	g->edgeMap().at(e).removeCellJunction(id);
}


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
	//simply destroy cell if it is on a boundary
	for (int e : edge_keys) 
	{
		if (g->edgeMap().at(e).getCellJunctions().size() < 2)
		{
			g->destroyCell(id);
			return;
		}
	}
	
	for (int v: vertex_keys) { 	if (g->vertexMap().at(v).edgeContacts().size() > 3) { return; } }
	calcCentroid();
	//detatch cell vertices and edges from cell
	for (int v : vertex_keys) { g->vertexMap().at(v).removeCellContact(id); }
	for (int e : edge_keys) { g->edgeMap().at(e).removeCellJunction(id); }
	
	int v_new = g->createVertex(centroid);
	//std::cout << "deleting cell: " << id << '\n';
	//std::cout << "v_new=" << v_new << '\n';
	
	std::vector<int> neighbour_cells;
	std::vector<int> v_new_index;
	
	
	for (int e : edge_keys)
	{
		//cell outside edge
		int c = *(g->edgeMap().at(e).getCellJunctions().begin());
		neighbour_cells.push_back(c);
		
		//tell vertices that it is no longer in contact with the removed edges
		for (int v : vertex_keys) 
		{ 
			g->vertexMap().at(v).removeEdgeContact(e); 
		}
		
		//disconnect edge from outside cell, edge will delete itself because it has no cell junctions
		g->cellMap().at(c).removeEdge(e);
		
		for (int v : vertex_keys)
		{
			auto it = std::find(g->cellMap().at(c).getVertices().begin(), g->cellMap().at(c).getVertices().end(),v);
			if (it != g->cellMap().at(c).getVertices().end())
			{
				v_new_index.push_back(std::distance(g->cellMap().at(c).getVertices().begin(), it));
				break;
			}

		}
	}
	for (int i = 0; i < neighbour_cells.size(); i++)
	{
		g->cellNewVertex(neighbour_cells[i], v_new, v_new_index[i]%(g->cellMap().at(neighbour_cells[i]).getVertices().size()));
	}
	
	for (int v : vertex_keys)
	{
		int e = *(g->vertexMap().at(v).edgeContacts().begin());		//only element left in the vertex edge_contacts (incident edge)

		g->edgeMap().at(e).swapVertex(v, v_new);	//reconnect incident edges so that they meet at centroid, this deletes the old vertex and tells cells it no longer has this vertex
		g->vertexMap().at(v_new).addEdgeContact(e);		//tell new vertex it has e as an edge contact
	}
	
	//std::cout << neighbour_cells.size() << ' ' << v_new_index.size() << '\n';
	/*for (int c : neighbour_cells)
	{
		std::cout << "c=" << c << "\tvertices: ";
		for (int v : g->cellMap().at(c).getVertices()) { std::cout << v << ' '; }
		std::cout << '\n';  
	}
	std::cout << '\n';*/
	g->cellMap().erase(id);
}


void Cell::calcCentroid()
{
    K::FT x_sum = 0;
    K::FT y_sum = 0;

    for (int v : vertex_keys) {
        x_sum += g->vertexMap().at(v).R().x();
        y_sum += g->vertexMap().at(v).R().y();
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
		double x_v = g->vertexMap().at(v).R().x();
		double y_v = g->vertexMap().at(v).R().y();
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
		A += g->vertexMap().at(vertex_keys[i]).R().x()*g->vertexMap().at(vertex_keys[(i+1)%vertex_keys.size()]).R().y()
			- g->vertexMap().at(vertex_keys[(i+1)%vertex_keys.size()]).R().x()*g->vertexMap().at(vertex_keys[i]).R().y();
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

	/*std::cout << '(' << id << ") ";
	for (int v : vertex_keys)
	{
		for (int c : g->vertexMap().at(v).cellContacts()) { std::cout << c << ' '; }
		std::cout << '\t';
	}
	std::cout << '\n';*/
	
	/*std::cout << '(' << id << ") ";
	for (int e : edge_keys)
	{
		for (int c : g->edgeMap().at(e).getCellJunctions()) { std::cout << c << ' '; }
		std::cout << '\t';
	}
	std::cout << '\n';*/
	
	//inefficient but ok for now
	//needs more testing
	std::vector<int> nearest_neighbour_cells;
	for (int e : edge_keys)
	{
		int c_ = -1;
		for (int c : g->edgeMap().at(e).getCellJunctions()) { if (c != id) { c_ = c; break; } }
		if (c_ != -1) 
		{
			nearest_neighbour_cells.push_back(c_);
		}
		else { m = 0; return; }
		//std::cout << c_ << ' ';
	} //std::cout << '\n';
	
	double w = 0;
	Vec D_old;
	Vec D = g->cellMap().at(nearest_neighbour_cells[0]).getDirector();
	
	for (int i = 1; i < nearest_neighbour_cells.size()+1; i++)
	{
		D_old = D;
		D = g->cellMap().at(nearest_neighbour_cells[i%nearest_neighbour_cells.size()]).getDirector();
		w += CGAL::angle(D_old, D);
	}
	m = 0.5*w/pi;
	std::cout << m << '\n';
}





