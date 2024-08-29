#include "cell.h"


Cell::Cell(Tissue* T, std::vector<int>& vertices, std::vector<int>& edges) : 
	T(T), id(T->c_c()), vertices(vertices), edges(edges)
{	
	A_ = 0;
	for (int i = 0; i < vertices.size(); i++) {
		A_ += T->vert(vertices[i]).r().x()*T->vert(vertices[(i+1)%vertices.size()]).r().y()
			- T->vert(vertices[(i+1)%vertices.size()]).r().x()*T->vert(vertices[i]).r().y();
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
	for (int e : edges) std::cout << T->edge(e).v1() << ' ' << T->edge(e).v2() << "   "; 
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
	for (int e : edges) { if (T->edge(e).cellJunctions().size() < 2) return true; }
	return false;
}

void Cell::addVertex(int v, int i) { vertices.insert(vertices.begin()+i, v); }
void Cell::removeVertex(int v) 
{ 
	std::vector<int>::const_iterator it = std::find(vertices.begin(), vertices.end(), v);
	if (it != vertices.end())
	{
		vertices.erase(it); 
		T->vert(v).removeCellContact(id);
	}
	
}
void Cell::addEdge(int e, int i) { edges.insert(edges.begin()+i, e); }
int Cell::removeEdge(int e) 
{ 
	std::vector<int>::const_iterator it = std::find(edges.begin(), edges.end(),e);
	if (it != edges.end())
	{
		edges.erase(it); 
		T->edge(e).removeCellJunction(id);
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
			T->vert(v_new).addCellContact(id);
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
		if (T->edge(edges[i]).l() > longest_l) 
		{
			longest_l = T->edge(edges[i]).l();
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
	Point a = CGAL::midpoint(T->vert(T->edge(e_a).v1()).r(), T->vert(T->edge(e_a).v2()).r());
	Point b = CGAL::midpoint(T->vert(T->edge(e_b).v1()).r(), T->vert(T->edge(e_b).v2()).r());
	int v_a = T->createVertex(a); int v_b = T->createVertex(b); 
	
	//add newly created vertices to relevent cells at correct index
	int c_a = -1; int c_b = -1; 
	int i_ea; int i_eb;
	auto addNewVertices = [this, c_a, c_b](int& c_x, const int& e_x, const int& v_x, int& i_ex) //happy
	{
		if (T->edge(e_x).cellJunctions().size() == 2)
		{
			auto it = T->edge(e_x).cellJunctions().begin();
			c_x = (*it == id) ? *std::next(it) : *it;
			
			//place vertex between e_xv1 and e_xv2
			const std::vector<int>& vertices_c_x = T->cell(c_x).Vertices();			
			std::vector<int>::const_iterator it1 = std::find(vertices_c_x.begin(), vertices_c_x.end(), T->edge(e_x).v1());
			std::vector<int>::const_iterator it2 = std::find(vertices_c_x.begin(), vertices_c_x.end(), T->edge(e_x).v2());
			if (it1 == vertices_c_x.begin() && it2 != (it1+1)) i_ex = 0;
			else if (it2 == vertices_c_x.begin() && it1 != (it2+1)) i_ex = 0;
			else i_ex = (it1 < it2) ? std::distance(vertices_c_x.begin(), it2) : std::distance(vertices_c_x.begin(), it1);
			T->cellNewVertex(c_x, v_x, i_ex);
		}
	};
	addNewVertices(c_a, e_a, v_a, i_ea); addNewVertices(c_b, e_b, v_b, i_eb);
			
			
	int e_new = T->createEdge(v_a, v_b); //edge dividng cell
	int e_a1, e_a2, e_b1, e_b2;
	//create edges so that e_x1 is always before e_x2 in vertex/edge ordering	
	auto updateAB = [this](const int c_x, const int e_x, const int v_x, int& e_x1, int& e_x2) //happy
	{
		const std::vector<int>& vertices_c_x = T->cell(c_x).Vertices();
		std::vector<int>::const_iterator it_v1 = std::find(vertices_c_x.begin(), vertices_c_x.end(), T->edge(e_x).v1());
		std::vector<int>::const_iterator it_v2 = std::find(vertices_c_x.begin(), vertices_c_x.end(), T->edge(e_x).v2());
		
		int i_vx1 = std::distance(vertices_c_x.begin(), it_v1); int i_vx2 = std::distance(vertices_c_x.begin(), it_v2);
		if (((i_vx2+2) % vertices_c_x.size()) == i_vx1) { std::swap(it_v1, it_v2); std::swap(i_vx1, i_vx2); }
		//std::cout << std::distance(vertices_c_x.begin(), it_v1) << ' ' << std::distance(vertices_c_x.begin(), it_v2) << '\n';
		e_x1 = T->createEdge(*it_v1, v_x);
		e_x2 = T->createEdge(*it_v2, v_x);
		if (c_x != -1)
		{
			const std::vector<int>& edges_c_x = T->cell(c_x).Edges();
			int i = T->cellRemoveEdge(c_x, e_x);
			if (i_vx2 != 1)
			{
				T->cellNewEdge(c_x, e_x2, i);
				T->cellNewEdge(c_x, e_x1, i);
			}
			else
			{
				T->cellNewEdge(c_x, e_x2, 0);
				T->cellNewEdge(c_x, e_x1, i+1);
			}
		}
	};
	updateAB(c_a, e_a, v_a, e_a1, e_a2); updateAB(c_b, e_b, v_b, e_b1, e_b2);
	
	
	auto isEdge = [this](int e, int v_1, int v_2) { return ( (v_1 == T->edge(e).v1() && v_2 == T->edge(e).v2()) || (v_1 == T->edge(e).v2() && v_2 == T->edge(e).v1()) ); };	
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
	
	int c_p = T->createCell(cell_p_vertices, cell_p_edges);
	int c_q = T->createCell(cell_q_vertices, cell_q_edges);
	
	//std::vector<int> neighbour_copy = neighbours;
	T->destroyCell(id);
	//for (int c : neighbour_copy) T->cell(c).findNeighbours();
	std::cout << "cell divided\n";
}

void Cell::extrude()
{
	
	//simply destroy cell if it is on a boundary
	if (onBoundary()) 
	{ 
		//std::vector<int> neighbours_copy = neighbours;
		T->destroyCell(id);
		//for (int c : neighbours_copy) T->cell(c).findNeighbours();
		return; 
	}
	for (int v: vertices) { if (T->vert(v).edgeContacts().size() > 3) return; }
	
	//std::vector<int> neighbours_copy = neighbours;
	calcR_0(); 														//calculate centroid, create vertex at centroid, and detatch cell vertices and edges from cell
	int v_new = T->createVertex(r_0_);
	for (int v : vertices) { T->vert(v).removeCellContact(id); }
	for (int e : edges) { T->edge(e).removeCellJunction(id); }

	for (int e : edges)
	{	
		int c = *(T->edge(e).cellJunctions().begin()); 				//cell outside edge
		
		//tell vertices that it is no longer in contact with the removed edges
		T->vert(T->edge(e).v1()).removeEdgeContact(e);
		T->vert(T->edge(e).v2()).removeEdgeContact(e);
	
		T->cell(c).removeEdge(e);									//disconnect edge from outside cell, edge will delete itself because it has no cell junctions	
	}
	for (int v : vertices)
	{
		int e = *(T->vert(v).edgeContacts().begin());				//only element left in the vertex edge_contacts (incident edge)
		T->edge(e).swapVertex(v, v_new);							//reconnect incident edges so that they meet at r_0, this deletes the old vertex and tells cells it no longer has this vertex
	}
	for (int c : T->vert(v_new).cellContacts()) { if (!(T->cell(c).valid())) T->cell(c).rotateVertices(); } //{ T->cell(c).outputVertices(); T->cell(c).outputEdgeVertices(); }}
	//for (int c : T->vert(v_new).cellContacts()) { if (!(T->cell(c).valid())) { T->cell(c).outputVertices(); T->cell(c).outputEdgeVertices(); }}
	
	vertices = {}; edges = {};
	T->destroyCell(id);
	//for (int c : neighbours_copy) T->cell(c).findNeighbours();
	std::cout << "cell extruded\n";
}


void Cell::calcR_0()
{
    double x_sum = 0;
    double y_sum = 0;
    for (int v : vertices) 
    {
        x_sum += T->vert(v).r().x();
        y_sum += T->vert(v).r().y();
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
		double x_v = T->vert(v).r().x();
		double y_v = T->vert(v).r().y();
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
		A_ += T->vert(vertices[i]).r().x()*T->vert(vertices[(i+1)%vertices.size()]).r().y()
			- T->vert(vertices[(i+1)%vertices.size()]).r().x()*T->vert(vertices[i]).r().y();
	} A_ *= 0.5;
}

void Cell::calcL()
{
	L_ = 0;
	//edge lengths must already be calculated
	for (int e : edges) L_ += T->edge(e).l();
}

void Cell::calcT_A() { T_A_ = param::K_a*(A_-param::A_0); }

void Cell::findNeighbours()
{
	neighbours = {};
	for (int i = 0; i < vertices.size(); i++)
	{
		int v_prev = vertices[(i-1+vertices.size())%vertices.size()];
		int v = vertices[i];
		int v_next = vertices[(i+1)%vertices.size()];
		
		const std::unordered_set<int>& cell_contacts_prev = T->vert(v_prev).cellContacts();
		const std::unordered_set<int>& cell_contacts = T->vert(v).cellContacts();
		const std::unordered_set<int>& cell_contacts_next = T->vert(v_next).cellContacts();
		size_t n = cell_contacts.size();
		
		for (int c : cell_contacts)
		{
			if (c != id)
			{
				if (n <= 3)
				{
					std::unordered_set<int>::const_iterator it_next = std::find(cell_contacts_next.begin(), cell_contacts_next.end(), c); //look for cell c in next vertex cell contacts
					if (it_next == cell_contacts_next.end()) neighbours.push_back(c); //add c to nearest neighbours if it is not in the next vertex edge contacts
				}
				else 
				{
					std::unordered_set<int>::const_iterator it_prev = std::find(cell_contacts_prev.begin(), cell_contacts_prev.end(), c);
					std::unordered_set<int>::const_iterator it_next = std::find(cell_contacts_next.begin(), cell_contacts_next.end(), c);
					if (it_next == cell_contacts_next.end() && it_prev == cell_contacts_prev.end()) neighbours.push_back(c);
				}
			}		
		}
	}
}

void Cell::calcm()
{
	findNeighbours(); 
	int n = neighbours.size();
	double w = 0;
	for (int i = 0; i < n; i++) w += T->D_angle(neighbours[i], neighbours[(i+1)%n]);
	m_ = w*boost::math::constants::one_div_two_pi <double>();
}


bool Cell::valid()
{
	if (vertices.size() != edges.size()) return false;
	for (int i = 0; i < vertices.size(); i++) if ((T->edge(edges[i]).hasVertex(vertices[i]) && T->edge(edges[i]).hasVertex(vertices[(i+1)%vertices.size()])) == false) return false;
	return true;
}










