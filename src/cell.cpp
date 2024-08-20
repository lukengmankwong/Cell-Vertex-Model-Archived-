#include "cell.h"

Cell::Cell(Global* g, std::vector<int>& vertices, std::vector<int>& edges) : 
	g(g), id(g->c_c()), vertices(vertices), edges(edges)
{	
	A_ = 0;
	for (int i = 0; i < vertices.size(); i++) {
		A_ += g->vert(vertices[i]).r().x()*g->vert(vertices[(i+1)%vertices.size()]).r().y()
			- g->vert(vertices[(i+1)%vertices.size()]).r().x()*g->vert(vertices[i]).r().y();
	} A_ *= 0.5;
	S_ = A_/std::fabs(A_);
}


void Cell::outputVertices() const 
{ 
	std::cout << "cell (" << id << ") vertices: "; 
	for (int v : vertices) std::cout << v << ' '; 
	std::cout << '\n';
}
void Cell::outputEdges() const 
{ 
	std::cout << "cell (" << id << ") edges: ";
	for (int e : edges) std::cout << e << ' '; 
	std::cout << '\n';
}
void Cell::outputEdgeVertices() const 
{ 
	std::cout << "cell (" << id << ") vertices by edge: ";
	for (int e : edges) std::cout << g->edge(e).v1() << ' ' << g->edge(e).v2() << "   "; 
	std::cout << '\n';
}


const std::vector<int>& Cell::Vertices() const { return vertices; }
const std::vector<int>& Cell::Edges() 	 const { return edges; }

const Point& 	Cell::r_0() const { return r_0_; }
const double 	Cell::A() 	const { return S_*A_; }
const double 	Cell::S() 	const { return S_; }
const double 	Cell::L() 	const { return L_; }
const double 	Cell::T_A() const { return T_A_; }
const Vec& 	 	Cell::n() 	const { return n_; }
const double 	Cell::Z() 	const { return Z_; }
const double 	Cell::X() 	const { return X_; }
const double 	Cell::m() 	const { return m_; }

const bool Cell::hasEdge(int e) const
{
	std::vector<int>::const_iterator it = std::find(edges.begin(), edges.end(), e);
	return it != edges.end();
}

const bool Cell::onBoundary() const
{
	for (int e : edges) { if (g->edge(e).cellJunctions().size() < 2) return true; }
	return false;
}

void Cell::addVertex(int v, int i) { vertices.insert(vertices.begin()+i, v); }
void Cell::removeVertex(int v) 
{ 
	std::vector<int>::const_iterator it = std::find(vertices.begin(), vertices.end(), v);
	if (it != vertices.end())
	{
		vertices.erase(it); 
		g->vert(v).removeCellContact(id);
	}
	
}
void Cell::addEdge(int e, int i) { edges.insert(edges.begin()+i, e); }
int Cell::removeEdge(int e) 
{ 
	std::vector<int>::const_iterator it = std::find(edges.begin(), edges.end(),e);
	if (it != edges.end())
	{
		edges.erase(it); 
		g->edge(e).removeCellJunction(id);
	}
	return std::distance(edges.cbegin(), it);
}


void Cell::exchangeVertex(int v_old, int v_new)
{	
	std::vector<int>::iterator it_new = std::find(vertices.begin(), vertices.end(), v_new);
	if (it_new == vertices.end()) 
	{ 
		std::vector<int>::iterator it_old = std::find(vertices.begin(), vertices.end(), v_old);
		if (it_old != vertices.end())
		{
			*it_old = v_new;
			g->vert(v_new).addCellContact(id);
		}
	}
	removeVertex(v_old);
}

void Cell::rotateVertices() { std::rotate(vertices.rbegin(), vertices.rbegin() + 1, vertices.rend());}

const int Cell::longestEdge_i() const
{
	double longest_l = 0; int i_l;
	for (int i= 0; i < edges.size(); i++) 
	{
		if (g->edge(edges[i]).l() > longest_l) 
		{
			longest_l = g->edge(edges[i]).l();
			i_l = i;
		}
	} return i_l;
}


void Cell::divide()
{
	if (vertices.size() <= 3 || onBoundary()) { return; }
	
	int i_va = longestEdge_i();
	int e_a = edges[i_va];
	int i_vb = (i_va+edges.size()/2)%edges.size();
	int e_b = edges[i_vb]; //edge opposite longest edge
	
	//midpoints of above edges
	Point a = CGAL::midpoint(g->vert(g->edge(e_a).v1()).r(), g->vert(g->edge(e_a).v2()).r());
	Point b = CGAL::midpoint(g->vert(g->edge(e_b).v1()).r(), g->vert(g->edge(e_b).v2()).r());
	//outputVertices();
	int v_a = g->createVertex(a); int v_b = g->createVertex(b); 
	//std::cout << "v_a: " << v_a << "\t\tv_b: " << v_b << '\n';
	
	//add newly created vertices to relevent cells at correct index
	int c_a = -1; int c_b = -1; 
	int i_ea; int i_eb;
	auto addNewVertices = [this, c_a, c_b](int& c_x, const int& e_x, const int& v_x, int& i_ex) //happy
	{
		if (g->edge(e_x).cellJunctions().size() == 2)
		{
			auto it = g->edge(e_x).cellJunctions().begin();
			c_x = (*it == id) ? *std::next(it) : *it;
			
			//place vertex between e_xv1 and e_xv2
			const std::vector<int>& vertices_c_x = g->cell(c_x).Vertices();			
			std::vector<int>::const_iterator it1 = std::find(vertices_c_x.begin(), vertices_c_x.end(), g->edge(e_x).v1());
			std::vector<int>::const_iterator it2 = std::find(vertices_c_x.begin(), vertices_c_x.end(), g->edge(e_x).v2());
			if (it1 == vertices_c_x.begin() && it2 != (it1+1)) i_ex = 0;
			else if (it2 == vertices_c_x.begin() && it1 != (it2+1)) i_ex = 0;
			else i_ex = (it1 < it2) ? std::distance(vertices_c_x.begin(), it2) : std::distance(vertices_c_x.begin(), it1);
			g->cellNewVertex(c_x, v_x, i_ex);
		}
	};
	addNewVertices(c_a, e_a, v_a, i_ea); addNewVertices(c_b, e_b, v_b, i_eb);
			
	int e_new = g->createEdge(v_a, v_b); //edge dividng cell
	int e_a1, e_a2, e_b1, e_b2;
	//create edges so that e_x1 is always before e_x2 in vertex/edge ordering	
	auto updateAB = [this](const int c_x, const int e_x, const int v_x, int& e_x1, int& e_x2) //happy
	{
		const std::vector<int>& vertices_c_x = g->cell(c_x).Vertices();
		std::vector<int>::const_iterator it_v1 = std::find(vertices_c_x.begin(), vertices_c_x.end(), g->edge(e_x).v1());
		std::vector<int>::const_iterator it_v2 = std::find(vertices_c_x.begin(), vertices_c_x.end(), g->edge(e_x).v2());
		
		int i_vx1 = std::distance(vertices_c_x.begin(), it_v1); int i_vx2 = std::distance(vertices_c_x.begin(), it_v2);
		if (((i_vx2+2) % vertices_c_x.size()) == i_vx1) { std::swap(it_v1, it_v2); std::swap(i_vx1, i_vx2); }
		//std::cout << std::distance(vertices_c_x.begin(), it_v1) << ' ' << std::distance(vertices_c_x.begin(), it_v2) << '\n';
		e_x1 = g->createEdge(*it_v1, v_x);
		e_x2 = g->createEdge(*it_v2, v_x);
		if (c_x != -1)
		{
			const std::vector<int>& edges_c_x = g->cell(c_x).Edges();
			int i = g->cellRemoveEdge(c_x, e_x);
			if (i_vx2 != 1)
			{
				g->cellNewEdge(c_x, e_x2, i);
				g->cellNewEdge(c_x, e_x1, i);
			}
			else
			{
				g->cellNewEdge(c_x, e_x2, 0);
				g->cellNewEdge(c_x, e_x1, i+1);
			}
		}
	};
	updateAB(c_a, e_a, v_a, e_a1, e_a2); updateAB(c_b, e_b, v_b, e_b1, e_b2);
	/*std::cout << "e_a1: " << g->edge(e_a1).v1() << ' ' << g->edge(e_a1).v2() << '\n';
	std::cout << "e_a2: " << g->edge(e_a2).v1() << ' ' << g->edge(e_a2).v2() << '\n';
	std::cout << "e_b1: " << g->edge(e_b1).v1() << ' ' << g->edge(e_b1).v2() << '\n';
	std::cout << "e_b2: " << g->edge(e_b2).v1() << ' ' << g->edge(e_b2).v2() << '\n';*/
	
	auto isEdge = [this](int e, int v_1, int v_2) { return ( (v_1 == g->edge(e).v1() && v_2 == g->edge(e).v2()) || (v_1 == g->edge(e).v2() && v_2 == g->edge(e).v1()) ); };	
	//find vertices of the two new cells in rotaional order
	std::vector<int> cell_p_vertices, cell_q_vertices; std::vector<int> cell_p_edges, cell_q_edges;
	bool f_q = !(i_va < i_vb);
	for (int i = 0; i < vertices.size(); i++) //neeed to fix
	{
		if (!f_q)
		{ 
			cell_p_vertices.push_back(vertices[i]);
			if (isEdge(e_a, vertices[i], vertices[(i+1)%vertices.size()]))
			{
				cell_p_vertices.push_back(v_a);
				cell_p_vertices.push_back(v_b);
				f_q = true;
			}
		}	
		else
		{
			cell_q_vertices.push_back(vertices[i]);
			if (isEdge(e_b, vertices[i], vertices[(i+1)%vertices.size()]))
			{
				cell_q_vertices.push_back(v_b);
				cell_q_vertices.push_back(v_a);		
				f_q = false;
			}
		}
	}
	
	//update edges of two new cells in rotational order
	for (int i = 0; i < cell_p_vertices.size(); i++)
	{
		if (isEdge(e_a2, cell_p_vertices[i], cell_p_vertices[(i+1)%cell_p_vertices.size()])) cell_p_edges.push_back(e_a2);
		else if (isEdge(e_b1, cell_p_vertices[i], cell_p_vertices[(i+1)%cell_p_vertices.size()])) cell_p_edges.push_back(e_b1);
		else if (isEdge(e_new, cell_p_vertices[i], cell_p_vertices[(i+1)%cell_p_vertices.size()])) cell_p_edges.push_back(e_new);
		else for (int e : edges) if (isEdge(e, cell_p_vertices[i], cell_p_vertices[(i+1)%cell_p_vertices.size()])) cell_p_edges.push_back(e);
	}
	for (int i = 0; i < cell_q_vertices.size(); i++)
	{
		if (isEdge(e_b2, cell_q_vertices[i], cell_q_vertices[(i+1)%cell_q_vertices.size()])) cell_q_edges.push_back(e_b2);
		else if (isEdge(e_a1, cell_q_vertices[i], cell_q_vertices[(i+1)%cell_q_vertices.size()])) cell_q_edges.push_back(e_a1);
		else if (isEdge(e_new, cell_q_vertices[i], cell_q_vertices[(i+1)%cell_q_vertices.size()])) cell_q_edges.push_back(e_new);
		else for (int e : edges) if (isEdge(e, cell_q_vertices[i], cell_q_vertices[(i+1)%cell_q_vertices.size()])) cell_q_edges.push_back(e);
	}
	
	int c_p = g->createCell(cell_p_vertices, cell_p_edges);
	int c_q = g->createCell(cell_q_vertices, cell_q_edges);
	
	/*std::cout << "AFTER:\n";
	g->cell(c_a).outputVertices(); g->cell(c_a).outputEdgeVertices();
	g->cell(c_b).outputVertices(); g->cell(c_b).outputEdgeVertices();
	g->cell(c_p).outputVertices(); g->cell(c_p).outputEdgeVertices();
	g->cell(c_q).outputVertices(); g->cell(c_q).outputEdgeVertices();
	std::cout << '\n';*/
	
	/*if (!(g->cell(c_a).valid())) { std::cout << "c_a invalid\n"; std::cin.get(); }
	if (!(g->cell(c_b).valid())) { std::cout << "c_b invalid\n"; std::cin.get(); }
	if (!(g->cell(c_p).valid())) { std::cout << "c_p invalid\n"; std::cin.get(); }
	if (!(g->cell(c_q).valid())) { std::cout << "c_q invalid\n"; std::cin.get(); }*/
	
	g->destroyCell(id);
}

void Cell::extrude()
{
	//simply destroy cell if it is on a boundary or is a triangle
	if (onBoundary()) { g->destroyCell(id); return; }
	for (int v: vertices) { if (g->vert(v).edgeContacts().size() > 3) return; }
	
	calcR_0(); 														//calculate centroid, detatch cell vertices and edges from cell, and create vertex at centroid
	for (int v : vertices) { g->vert(v).removeCellContact(id); }
	for (int e : edges) { g->edge(e).removeCellJunction(id); }
	int v_new = g->createVertex(r_0_);

	for (int e : edges)
	{	
		int c = *(g->edge(e).cellJunctions().begin()); 				//cell outside edge
		
		//tell vertices that it is no longer in contact with the removed edges
		g->vert(g->edge(e).v1()).removeEdgeContact(e);
		g->vert(g->edge(e).v2()).removeEdgeContact(e);
	
		g->cell(c).removeEdge(e);									//disconnect edge from outside cell, edge will delete itself because it has no cell junctions	
	}
	for (int v : vertices)
	{
		int e = *(g->vert(v).edgeContacts().begin());				//only element left in the vertex edge_contacts (incident edge)
		g->edge(e).swapVertex_rep(v, v_new);							//reconnect incident edges so that they meet at r_0, this deletes the old vertex and tells cells it no longer has this vertex
	}
	for (int c : g->vert(v_new).cellContacts()) { if (!(g->cell(c).valid())) g->cell(c).rotateVertices(); } //{ g->cell(c).outputVertices(); g->cell(c).outputEdgeVertices(); }}
	//for (int c : g->vert(v_new).cellContacts()) { if (!(g->cell(c).valid())) { g->cell(c).outputVertices(); g->cell(c).outputEdgeVertices(); }}
	
	g->cellMap().erase(id);
}


void Cell::calcR_0()
{
    double x_sum = 0;
    double y_sum = 0;
    for (int v : vertices) 
    {
        x_sum += g->vert(v).r().x();
        y_sum += g->vert(v).r().y();
    }
    r_0_ = Point(x_sum/vertices.size(), y_sum/vertices.size());
}

void Cell::calcG()
{
	G[0]=0;	G[1]=0; G[2]=0;
	calcR_0();
	double x_0 = r_0_.x(); double y_0 = r_0_.y();
	for (int v : vertices)
	{
		double x_v = g->vert(v).r().x();
		double y_v = g->vert(v).r().y();
		G[0]+=(x_v-x_0)*(x_v-x_0);	
		G[1]+=(x_v-x_0)*(y_v-y_0);
		G[2]+=(y_v-y_0)*(y_v-y_0);
	}
	
	double f = 1.0/vertices.size();
	G[0]*=f; G[1]*=f; G[2]*=f;
	
	lambda = 0.5*( G[0]+G[2] + std::sqrt( (G[0]+G[2])*(G[0]+G[2]) - 4*(G[0]*G[2]-G[1]*G[1]) ) );
	n_ = Vec(1, (lambda-G[0])/G[1]) / std::sqrt( 1 + ((lambda-G[0])/G[1])*((lambda-G[0])/G[1]) );
	
	Z_ = 0.5*(G[0]-G[2]);
	X_ = G[1];
}

void Cell::calcA()
{
	A_ = 0;
	for (int i = 0; i < vertices.size(); i++) {
		A_ += g->vert(vertices[i]).r().x()*g->vert(vertices[(i+1)%vertices.size()]).r().y()
			- g->vert(vertices[(i+1)%vertices.size()]).r().x()*g->vert(vertices[i]).r().y();
	} A_ *= 0.5;
}

void Cell::calcL()
{
	L_ = 0;
	//edge lengths must already be calculated
	for (int e : edges) L_ += g->edge(e).l();
}

void Cell::calcT_A() { T_A_ = param::K_a*(A_-param::A_0); }

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
	std::vector<int> nn_cells = nearestNeighbours(); int n = nn_cells.size();
	double w = 0;
	for (int i = 0; i < nn_cells.size(); i++) w += g->D_angle(nn_cells[i], nn_cells[(i+1)%n]);
	m_ = w*param::HALF_INV_PI;
}


bool Cell::valid()
{
	if (vertices.size() != edges.size()) return false;
	for (int i = 0; i < vertices.size(); i++) if ((g->edge(edges[i]).hasVertex(vertices[i]) && g->edge(edges[i]).hasVertex(vertices[(i+1)%vertices.size()])) == false) return false;
	return true;
}










