#include "cell.h"
#include <iomanip>

Cell::Cell(Global* g, std::vector<int>& vertices, std::vector<std::pair<int,int>>& edges) : 
	g(g), id(g->cellCounter()), vertices(vertices)
{
	for (std::pair<int, int> edge : edges)
	{
		this->edges.push_back(edge.first);
		this->edge_directions.push_back(edge.second);
	}
	
	A = 0;
	for (int i = 0; i < vertices.size(); i++) {
		A += g->vert(vertices[i]).R().x()*g->vert(vertices[(i+1)%vertices.size()]).R().y()
			- g->vert(vertices[(i+1)%vertices.size()]).R().x()*g->vert(vertices[i]).R().y();
	} A *= 0.5;
	S = A/std::fabs(A);
}

const std::vector<int>& Cell::Vertices() const { return vertices; }
const std::vector<int>& Cell::Edges() const { return edges; }

const double Cell::getA() const { return S*A; }
const double Cell::getS() const { return S; }
const double Cell::getL() const { return L; }
const double Cell::getT_A() const { return T_A; }
const Point& Cell::getCentroid() const { return centroid; }
const Vec& Cell::getDirector() const { return director; }
const double Cell::getm() const { return m; }

const bool Cell::hasEdge(int e) const
{
	auto it = std::find(edges.begin(), edges.end(), e);
	return it != edges.end();
}


void Cell::addVertex(int v, int i) { vertices.insert(vertices.begin()+i, v); }
void Cell::removeVertex(int v) { vertices.erase(std::find(vertices.begin(), vertices.end(), v)); }
void Cell::addEdge(int e, int i) { edges.insert(edges.begin()+i, e); }
void Cell::removeEdge(int e) 
{ 
	edges.erase(std::find(edges.begin(), edges.end(), e)); 
	g->edge(e).removeCellJunction(id);
}

void Cell::removeVertices() { for (int v : vertices) { g->vert(v).removeCellContact(id); } }
void Cell::removeEdges() { for (int e : edges) { g->edge(e).removeCellJunction(id); } }

const int Cell::longestEdge_i() const //for cell division along longest edges
{
	double longest_l = 0;
	int i_l;
	for (int i= 0; i < edges.size(); i++) 
	{
		if (g->edge(edges[i]).getl() > longest_l) 
		{
			longest_l = g->edge(edges[i]).getl();
			i_l = i;
		}
	} return i_l;
}


void Cell::divide()
{
	if (vertices.size() == 3) { return; }
	int i_l = longestEdge_i();
	int e_a = edges[i_l];
	int n = edges.size();
	int e_b = edges[(i_l+n/2)%edges.size()]; //edge opposite longest edge
	
	//midpoints of above edges
	Point a = CGAL::midpoint(g->vert(g->edge(e_a).E().first).R(), g->vert(g->edge(e_a).E().second).R());
	Point b = CGAL::midpoint(g->vert(g->edge(e_b).E().first).R(), g->vert(g->edge(e_b).E().second).R());
	
	int v_a = g->createVertex(a); int v_b = g->createVertex(b); //create vertices at midpoints of edges
	
	int c_a = -1; int c_b = -1;
	int i_a, i_b;
	if (g->edge(e_a).cellJunctions().size()==2)
	{
		auto it = g->edge(e_a).cellJunctions().begin();
		c_a = (*it == id) ? *std::next(it) : *it;
		/*std::cout << "c_a vertices: ";
		for (int v : g->cell(c_a).Vertices()) { std::cout << v << ' '; } std::cout <<'\n';*/
		
		std::vector<int>::const_iterator it1 = std::find(g->cell(c_a).Vertices().begin(), g->cell(c_a).Vertices().end(), g->edge(e_a).E().first);
		std::vector<int>::const_iterator it2 = std::find(g->cell(c_a).Vertices().begin(), g->cell(c_a).Vertices().end(), g->edge(e_a).E().second);
		if ( (it1 == g->cell(c_a).Vertices().begin() && it2 != (it1+1)) || (it2 == g->cell(c_a).Vertices().begin() && it1 != (it2+1)) ) i_a = 0;
		else i_a = (it1 < it2) ? std::distance(g->cell(c_a).Vertices().begin(), it2) : std::distance(g->cell(c_a).Vertices().begin(), it1);
		//std::cout << std::distance(g->cell(c_a).Vertices().begin(), it1) << ' ' << std::distance(g->cell(c_a).Vertices().begin(), it2) << ' ' << i_a <<  '\n';
		g->cellNewVertex(c_a, v_a, i_a);	
		
		
		
	}
	if (g->edge(e_b).cellJunctions().size()==2)
	{
		auto it = g->edge(e_b).cellJunctions().begin();
		c_b = (*it == id) ? *std::next(it) : *it;
		
		std::vector<int>::const_iterator it1 = std::find(g->cell(c_b).Vertices().begin(), g->cell(c_b).Vertices().end(), g->edge(e_b).E().first);
		std::vector<int>::const_iterator it2 = std::find(g->cell(c_b).Vertices().begin(), g->cell(c_b).Vertices().end(), g->edge(e_b).E().second);
		if ( (it1 == g->cell(c_b).Vertices().begin() && it2 != (it1+1)) || (it2 == g->cell(c_b).Vertices().begin() && it1 != (it2+1)) ) i_b = 0;
		else i_b = (it1 < it2) ? std::distance(g->cell(c_b).Vertices().begin(), it2) : std::distance(g->cell(c_b).Vertices().begin(), it1);
		//std::cout << std::distance(g->cell(c_b).Vertices().begin(), it1) << ' ' << std::distance(g->cell(c_b).Vertices().begin(), it2) << ' ' << i_b <<  '\n';
		g->cellNewVertex(c_b, v_b, i_b);
	}

	//find vertices of the two new cells in rotaional order
	std::vector<int> cell_p_vertices, cell_q_vertices;
	bool f_q = false;
	for (int i = 0; i < vertices.size(); i++)
	{
		if (!f_q)
		{ 
			cell_p_vertices.push_back(vertices[i]);
			if ( (vertices[i]==g->edge(e_a).E().first && vertices[(i+1)%vertices.size()]==g->edge(e_a).E().second) ||
				(vertices[i]==g->edge(e_a).E().second && vertices[(i+1)%vertices.size()]==g->edge(e_a).E().first) )
			{
				cell_p_vertices.push_back(v_a);
				cell_p_vertices.push_back(v_b);
				f_q = true;
			}
			else if ( (vertices[i]==g->edge(e_b).E().first && vertices[(i+1)%vertices.size()]==g->edge(e_b).E().second) ||
				(vertices[i]==g->edge(e_b).E().second && vertices[(i+1)%vertices.size()]==g->edge(e_b).E().first) )
			{
				cell_p_vertices.push_back(v_b);
				cell_p_vertices.push_back(v_a);
				f_q = true;			
			}
		}	
		else
		{
			cell_q_vertices.push_back(vertices[i]);
			if ( (vertices[i]==g->edge(e_a).E().first && vertices[(i+1)%vertices.size()]==g->edge(e_a).E().second) ||
				(vertices[i]==g->edge(e_a).E().second && vertices[(i+1)%vertices.size()]==g->edge(e_a).E().first) )
			{
				cell_q_vertices.push_back(v_a);
				cell_q_vertices.push_back(v_b);
				f_q = false;
			}
			else if ( (vertices[i]==g->edge(e_b).E().first && vertices[(i+1)%vertices.size()]==g->edge(e_b).E().second) ||
				(vertices[i]==g->edge(e_b).E().second && vertices[(i+1)%vertices.size()]==g->edge(e_b).E().first) )
			{
				cell_q_vertices.push_back(v_b);
				cell_q_vertices.push_back(v_a);
				f_q = false;
			}
		}
	}
	/*std::cout << "cell vertices: "; for (int v : vertices) { std::cout << v << ' '; } std::cout << '\n';
	std::cout << "p vertices: "; for (int v : cell_p_vertices) { std::cout << v << ' '; } std::cout << '\n';
	std::cout << "q vertices: "; for (int v : cell_q_vertices) { std::cout << v << ' '; } std::cout << "\n";*/
	
	int e_new = g->createEdge(v_a, v_b); //edge dividng cell
	int e_a1 = g->createEdge(g->edge(e_a).E().first, v_a); int e_a2 = g->createEdge(g->edge(e_a).E().second, v_a);
	int e_b1 = g->createEdge(g->edge(e_b).E().first, v_b); int e_b2 = g->createEdge(g->edge(e_b).E().second, v_b);
	
	std::vector<std::pair<int, int>> cell_p_edges, cell_q_edges;
	f_q = false;
	bool p_for, q_for;
	for (int i = 0; i < edges.size(); i++)
	{
		if (!f_q)
		{
			if (edges[i] == e_a)
			{
				if (vertices[i] == g->edge(e_a).E().first) 
				{ 
					cell_p_edges.push_back(std::make_pair(e_a1, 0));
					cell_p_edges.push_back(std::make_pair(e_new, 0));
					cell_q_edges.push_back(std::make_pair(e_a2, 0));
				}
				else if (vertices[i] == g->edge(e_a).E().second)
				{ 
					cell_p_edges.push_back(std::make_pair(e_a2, 1));
					cell_p_edges.push_back(std::make_pair(e_new, 0));
					cell_q_edges.push_back(std::make_pair(e_a1, 1));
				}
				f_q = true;
			}
			else if (edges[i] == e_b)
			{
				if (vertices[i] == g->edge(e_b).E().first) 
				{ 
					cell_p_edges.push_back(std::make_pair(e_b1, 0));
					cell_p_edges.push_back(std::make_pair(e_new, 1));
					cell_q_edges.push_back(std::make_pair(e_b2, 0));
					
				}
				else if (vertices[i] == g->edge(e_b).E().second)
				{ 
					cell_p_edges.push_back(std::make_pair(e_b2, 1));
					cell_p_edges.push_back(std::make_pair(e_new, 1));
					cell_q_edges.push_back(std::make_pair(e_b1, 1));
				}
				f_q = true;
			}
			else { cell_p_edges.push_back(std::make_pair(edges[i], edge_directions[i])); }
		}
		else
		{
			if (edges[i] == e_a)
			{
				if (vertices[i] == g->edge(e_a).E().first) 
				{ 
					cell_q_edges.push_back(std::make_pair(e_a1, 0));
					cell_q_edges.push_back(std::make_pair(e_new, 1));
					cell_p_edges.push_back(std::make_pair(e_a2, 0));
				}
				else if (vertices[i] == g->edge(e_a).E().second)
				{ 
					cell_q_edges.push_back(std::make_pair(e_a2, 1));
					cell_q_edges.push_back(std::make_pair(e_new, 1));
					cell_p_edges.push_back(std::make_pair(e_a1, 1));
				}
				f_q = false;
			}
			else if (edges[i] == e_b)
			{
				if (vertices[i] == g->edge(e_b).E().first) 
				{ 
					cell_q_edges.push_back(std::make_pair(e_b1, 0));
					cell_q_edges.push_back(std::make_pair(e_new, 0));
					cell_p_edges.push_back(std::make_pair(e_b2, 0));
				}
				else if (vertices[i] == g->edge(e_b).E().second)
				{ 
					cell_q_edges.push_back(std::make_pair(e_b2, 1));
					cell_q_edges.push_back(std::make_pair(e_new, 0));
					cell_p_edges.push_back(std::make_pair(e_b1, 1));
				}
				f_q = false;
			}
			else { cell_q_edges.push_back(std::make_pair(edges[i], edge_directions[i])); }
		}
	}
	
	/*std::cout << "cell edges: "; for (int e : edges) { std::cout << e << ' '; } std::cout << '\n';
	std::cout << "p edges: "; for (auto e : cell_p_edges) { std::cout << e.first << ' '; } std::cout << '\n';
	std::cout << "q edges: "; for (auto e : cell_q_edges) { std::cout << e.first << ' '; } std::cout << "\n";*/

	
	if (c_a != -1)
	{

		if (g->edge(e_a).E().first == g->cell(c_a).Vertices()[(i_a+1)%g->cell(c_a).Vertices().size()])
		{
			//std::cout << "first: " <<  g->edge(e_a).E().first << "\n";
			g->cellNewEdge(c_a, e_a1, i_a);
			g->cellNewEdge(c_a, e_a2, i_a);
		}
		else if (g->edge(e_a).E().second == g->cell(c_a).Vertices()[(i_a+1)%g->cell(c_a).Vertices().size()])
		{
			//std::cout << "second: " << g->edge(e_a).E().second << "\n";
			g->cellNewEdge(c_a, e_a2, i_a);
			g->cellNewEdge(c_a, e_a1, i_a);
		}
		g->cellRemoveEdge(c_a, e_a);
	}
	if (c_b != -1)
	{

		if (g->edge(e_b).E().first == g->cell(c_b).Vertices()[(i_b+1)%g->cell(c_b).Vertices().size()])
		{
			//std::cout << "first: " <<  g->edge(e_b).E().first << "\n";
			g->cellNewEdge(c_b, e_b1, i_b);
			g->cellNewEdge(c_b, e_b2, i_b);
		}
		else if (g->edge(e_b).E().second == g->cell(c_b).Vertices()[(i_b+1)%g->cell(c_b).Vertices().size()])
		{
			//std::cout << "second: " << g->edge(e_b).E().second << "\n";
			g->cellNewEdge(c_b, e_b2, i_b);
			g->cellNewEdge(c_b, e_b1, i_b);
		}
		g->cellRemoveEdge(c_b, e_b);
	}
	
	
	/*if (c_a != -1)
	{
		std::cout << "c_a vertices: "; for (int v : g->cell(c_a).Vertices()) { std::cout << v << ' '; } std::cout << '\n';
		std::cout << "c_a edge-vertices: "; for (int e : g->cell(c_a).Edges()) { std::cout << g->edge(e).E().first << ' ' << g->edge(e).E().second << "   "; } std::cout << '\n';
	}
	std::cout << '\n';*/
	g->createCell(cell_p_vertices, cell_p_edges);
	g->createCell(cell_q_vertices, cell_q_edges);
	g->destroyCell(id);
	
	
}

void Cell::extrude()
{
	//simply destroy cell if it is on a boundary
	for (int e : edges) 
	{
		if (g->edge(e).cellJunctions().size() < 2)
		{
			g->destroyCell(id);
			return;
		}
	}
	
	for (int v: vertices) { if (g->vert(v).edgeContacts().size() > 3) { return; } }
	calcCentroid();
	//detatch cell vertices and edges from cell
	for (int v : vertices) { g->vert(v).removeCellContact(id); }
	for (int e : edges) { g->edge(e).removeCellJunction(id); }
	
	int v_new = g->createVertex(centroid);
	//std::cout << "deleting cell: " << id << '\n';
	//std::cout << "v_new=" << v_new << '\n';
	
	std::vector<int> neighbour_cells;
	std::vector<int> v_new_index;
	
	
	for (int e : edges)
	{
		//cell outside edge
		int c = *(g->edge(e).cellJunctions().begin());
		neighbour_cells.push_back(c);
		
		//tell vertices that it is no longer in contact with the removed edges
		g->vert(g->edge(e).E().first).removeEdgeContact(e);
		g->vert(g->edge(e).E().second).removeEdgeContact(e);

		//disconnect edge from outside cell, edge will delete itself because it has no cell junctions
		g->cell(c).removeEdge(e);
		
		for (int v : vertices)
		{
			auto it = std::find(g->cell(c).Vertices().begin(), g->cell(c).Vertices().end(),v);
			if (it != g->cell(c).Vertices().end())
			{
				v_new_index.push_back(std::distance(g->cell(c).Vertices().begin(), it));
				break;
			}
		}
	}
	
	for (int i = 0; i < neighbour_cells.size(); i++)
	{
		g->cellNewVertex(neighbour_cells[i], v_new, v_new_index[i]%(g->cell(neighbour_cells[i]).Vertices().size()));
	}
	
	for (int v : vertices)
	{
		int e = *(g->vert(v).edgeContacts().begin());		//only element left in the vertex edge_contacts (incident edge)

		g->edge(e).swapVertex(v, v_new);	//reconnect incident edges so that they meet at centroid, this deletes the old vertex and tells cells it no longer has this vertex
		g->vert(v_new).addEdgeContact(e);		//tell new vertex it has e as an edge contact
	}
	
	//std::cout << neighbour_cells.size() << ' ' << v_new_index.size() << '\n';
	/*for (int c : neighbour_cells)
	{
		std::cout << "c=" << c << "\tvertices: ";
		for (int v : g->cell(c).Vertices()) { std::cout << v << ' '; }
		std::cout << '\n';  
	}
	std::cout << '\n';*/
	g->cellMap().erase(id);
}


void Cell::calcCentroid()
{
    K::FT x_sum = 0;
    K::FT y_sum = 0;

    for (int v : vertices) {
        x_sum += g->vert(v).R().x();
        y_sum += g->vert(v).R().y();
    }

    centroid = Point(x_sum/vertices.size(), y_sum/vertices.size());
}

void Cell::calcG()
{
	G[0]=0;	G[1]=0;
			G[2]=0;
	
	this->calcCentroid();
	double x_0 = centroid.x(); double y_0 = centroid.y();
	for (int v : vertices)
	{
		double x_v = g->vert(v).R().x();
		double y_v = g->vert(v).R().y();
		G[0]+=(x_v-x_0)*(x_v-x_0);	G[1]+=(x_v-x_0)*(y_v-y_0);
									G[2]+=(y_v-y_0)*(y_v-y_0);
	}
	
	double f = 1.0/vertices.size();
	G[0]*=f; 	G[1]*=f;
				G[2]*=f;
	
	lambda = 0.5*( G[0]+G[2] + std::sqrt( (G[0]+G[2])*(G[0]+G[2]) - 4*(G[0]*G[2]-G[1]*G[1]) ) );
	director = Vec(1, (lambda-G[0])/G[1]) / std::sqrt( 1 + ((lambda-G[0])/G[1])*((lambda-G[0])/G[1]) );
	
	//calculation of traceless of gyration tensor
	double tr_less = 0.5*(G[0]+G[2]);
	TL[0]=G[0]-tr_less;
	TL[1]=G[1];
	TL[2]=G[2]-tr_less;
}

void Cell::calcA()
{
	A = 0;
	for (int i = 0; i < vertices.size(); i++) {
		A += g->vert(vertices[i]).R().x()*g->vert(vertices[(i+1)%vertices.size()]).R().y()
			- g->vert(vertices[(i+1)%vertices.size()]).R().x()*g->vert(vertices[i]).R().y();
	} A *= 0.5;
}

void Cell::calcL()
{
	L = 0;
	//edge lengths must already be calculated
	for (int e : edges) { L += g->edge(e).getl(); }
}

void Cell::calcT_A() { T_A = k_A*(A-A_0); }

std::vector<int> Cell::nearestNeighbours()
{
	std::vector<int> nearest_neighbours;
	for (int i = 0; i < vertices.size(); i++)
	{
		int v_prev = vertices[(i-1+vertices.size())%vertices.size()];
		int v = vertices[i];
		int v_next = vertices[(i+1)%vertices.size()];
		for (int c : g->vert(v).cellContacts())
		{
			if (c != id)
			{
				if (g->vert(v).cellContacts().size() <= 3)
				{
					//look for cell c in next vertex edge contacts
					auto it_next = std::find(g->vert(v_next).cellContacts().begin(), g->vert(v_next).cellContacts().end(), c);
					//add c to nearest neighbours if it is not in the next vertex edge contacts
					if (it_next == g->vert(v_next).cellContacts().end()) { nearest_neighbours.push_back(c); }
				}
				else 
				{
					auto it_prev = std::find(g->vert(v_prev).cellContacts().begin(), g->vert(v_prev).cellContacts().end(), c);
					auto it_next = std::find(g->vert(v_next).cellContacts().begin(), g->vert(v_next).cellContacts().end(), c);
					if (it_next == g->vert(v_next).cellContacts().end() && it_prev == g->vert(v_prev).cellContacts().end()) { nearest_neighbours.push_back(c); }
				}
			}		
		}
	}
	return nearest_neighbours;
}

void Cell::calcm()
{
	std::vector<int> nearest_neighbour_cells = nearestNeighbours();
	
	double w = 0;
}

