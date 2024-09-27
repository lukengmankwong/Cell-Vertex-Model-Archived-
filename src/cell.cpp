#include "cell.h"
#include "tissue.h"


Cell::Cell(Tissue* T, std::vector<Vertex*>& vertices_, std::vector<Edge*>& edges_) : 
	T(T), vertices_(vertices_), edges_(edges_)
{	
	vertices_.reserve(8);
	edges_.reserve(8);
	neighbours_.reserve(8);
	
	A_ = 0;
	size_t n = vertices_.size();
	for (int i = 0; i < n; i++)
	{
		Point r_i = vertices_[i]->r();
		Point r_j = vertices_[(i+1)%n]->r();
		A_ += r_i.x()*r_j.y() - r_j.x()*r_i.y();
	} A_*= 0.5; S_ = A_/std::fabs(A_);
}
Cell::Cell() = default;


const Point& 	Cell::r_0() const { return r_0_; }
const double 	Cell::A() 	const { return S_*A_; }
const double 	Cell::S() 	const { return S_; }
const double 	Cell::L() 	const { return L_; }
const double 	Cell::T_A() const { return T_A_; }
const Vec& 	 	Cell::n() 	const { return n_; }
const double 	Cell::Z() 	const { return Z_; }
const double 	Cell::X() 	const { return X_; }
const double 	Cell::m() 	const { return m_; }

const std::vector<Vertex*>& Cell::vertices() 	const { return vertices_; }
const std::vector<Edge*>& Cell::edges()			const { return edges_; }
const std::vector<Cell*>& Cell::neighbours()	const { return neighbours_; }


void Cell::addVertex(Vertex* v, int i) { vertices_.insert(vertices_.begin()+i, v); }
void Cell::removeVertex(Vertex* v) 
{ 
	std::vector<Vertex*>::const_iterator it = std::find(vertices_.begin(), vertices_.end(), v);
	if (it != vertices_.end())
	{
		vertices_.erase(it); 
		v->removeCellContact(this);
	}
	
}
void Cell::exchangeVertex(Vertex* v_old, Vertex* v_new)
{	
	std::vector<Vertex*>::iterator it_new = std::find(vertices_.begin(), vertices_.end(), v_new);
	if (it_new == vertices_.end()) 
	{ 
		std::vector<Vertex*>::iterator it_old = std::find(vertices_.begin(), vertices_.end(), v_old);
		if (it_old != vertices_.end())
		{
			*it_old = v_new;
			v_new->addCellContact(this);
		}
	}
	removeVertex(v_old);
}
void Cell::rotateVertices() { std::rotate(vertices_.rbegin(), vertices_.rbegin() + 1, vertices_.rend());}

void Cell::addEdge(Edge* e, int i) { edges_.insert(edges_.begin()+i, e); }
int Cell::removeEdge(Edge* e) 
{ 
	std::vector<Edge*>::const_iterator it = std::find(edges_.begin(), edges_.end(), e);
	if (it != edges_.end())
	{
		edges_.erase(it); 
		e->removeCellJunction(this);
	}
	return std::distance(edges_.cbegin(), it);
}


const bool Cell::hasEdge(Edge* e) const
{
	std::vector<Edge*>::const_iterator it = std::find(edges_.begin(), edges_.end(), e);
	return it != edges_.end();
}
const bool Cell::onBoundary() const
{
	for (Edge* e : edges_) { if (e->cellJunctions().size() < 2) return true; }
	return false;
}

void Cell::findNeighbours()
{
	neighbours_ = {};
	std::unordered_set<Cell*> seen_cells;
	std::vector<std::pair<Cell*, double>> neighbour_cells_vec;

	for (Vertex* v : vertices_) 
	{
		for (Cell* c : v->cellContacts()) 
		{
			if (c == this || !seen_cells.insert(c).second) continue;
			c->calcR_0();
			Vec vec = c->r_0() - r_0_;
			double theta = std::atan2(vec.y(), vec.x());
			neighbour_cells_vec.push_back({c, theta});
		}
	}

	std::sort(neighbour_cells_vec.begin(), neighbour_cells_vec.end(), 
		[](const std::pair<Cell*, double>& c1, const std::pair<Cell*, double>& c2) { return c1.second < c2.second; });

	for (const std::pair<Cell*, double>& ce : neighbour_cells_vec) neighbours_.push_back(ce.first);
}

const int Cell::longestEdge_i() const
{
	double longest_l = 0; int i_l;
	for (int i= 0; i < edges_.size(); i++) 
	{
		if (edges_[i]->l() > longest_l) 
		{
			longest_l = edges_[i]->l();
			i_l = i;
		}
	} return i_l;
}


void Cell::extrude()
{
	
	//simply destroy cell if it is on a boundary
	if (onBoundary()) 
	{ 
		std::vector<Cell*> neighbours_copy = neighbours_;
		std::vector<Vertex*> vertices_copy = vertices_;
		T->destroyCell(this);
		//update neighbours and contacts orders
		for (Cell* c : neighbours_copy) c->findNeighbours();
		for (Vertex* v : vertices_copy) if (T->v_alive(v)) v->orderCellContacts();
		return; 
	}
	for (Vertex* v : vertices_) { if (v->edgeContacts().size() > 3) return; }
	
	calcR_0(); 													//calculate centroid, create vertex at centroid, and detatch cell vertices and edges from cell
	Vertex* v_new = T->createVertex(r_0_);
	for (Vertex* v : vertices_) { v->removeCellContact(this); }
	for (Edge* e : edges_) { e->removeCellJunction(this); }

	for (Edge* e : edges_)
	{	
		Cell* c = *(e->cellJunctions().begin()); 				//cell outside edge
		
		//tell vertices that it is no longer in contact with the removed edges
		e->v1()->removeEdgeContact(e);
		e->v2()->removeEdgeContact(e);
	
		c->removeEdge(e);										//disconnect edge from outside cell, edge will delete itself because it has no cell junctions	
	}
	for (Vertex* v : vertices_)
	{
		Edge* e = *(v->edgeContacts().begin());				//only element left in the vertex edge_contacts (incident edge)
		e->swapVertex(v, v_new);							//reconnect incident edges so that they meet at r_0, this deletes the old vertex and tells cells it no longer has this vertex
	}
	for (Cell* c : v_new->cellContacts()) while (!(c->valid())) c->rotateVertices();
	
	
	std::vector<Cell*> neighbours_copy = neighbours_;
	vertices_ = {}; edges_ = {};
	T->destroyCell(this);
	v_new->orderCellContacts();
	for (Cell* c : neighbours_copy) c->findNeighbours();
	
	std::cout << "cell extruded\n";
}

void Cell::divide()
{
	if (vertices_.size() <= 3 || onBoundary()) { return; }
	
	int i_va = longestEdge_i();
	Edge* e_a = edges_[i_va];
	int i_vb = (i_va+edges_.size()/2)%edges_.size();
	Edge* e_b = edges_[i_vb]; //edge opposite longest edge
	
	//midpoints of above edges
	Point a = CGAL::midpoint(e_a->v1()->r(), e_a->v2()->r());
	Point b = CGAL::midpoint(e_b->v1()->r(), e_b->v2()->r());
	Vertex* v_a = T->createVertex(a); Vertex* v_b = T->createVertex(b); 
	
	//add newly created vertices to relevent cells at correct index
	Cell* c_a = nullptr; Cell* c_b = nullptr; 
	int i_ea; int i_eb;
	auto addNewVertices = [this, c_a, c_b](Cell*& c_x, Edge* e_x, Vertex* v_x, int& i_ex)
	{
		if (e_x->cellJunctions().size() == 2)
		{
			auto it = e_x->cellJunctions().begin();
			c_x = (*it == this) ? *std::next(it) : *it;
			
			//place vertex between e_xv1 and e_xv2
			const std::vector<Vertex*>& c_x_vertices = c_x->vertices();			
			std::vector<Vertex*>::const_iterator it1 = std::find(c_x_vertices.begin(), c_x_vertices.end(), e_x->v1());
			std::vector<Vertex*>::const_iterator it2 = std::find(c_x_vertices.begin(), c_x_vertices.end(), e_x->v2());
			if (it1 == c_x_vertices.begin() && it2 != (it1+1)) i_ex = 0;
			else if (it2 == c_x_vertices.begin() && it1 != (it2+1)) i_ex = 0;
			else i_ex = (it1 < it2) ? std::distance(c_x_vertices.begin(), it2) : std::distance(c_x_vertices.begin(), it1);
			T->cellNewVertex(c_x, v_x, i_ex);
		}
	};
	addNewVertices(c_a, e_a, v_a, i_ea); addNewVertices(c_b, e_b, v_b, i_eb);

	Edge* e_new = T->createEdge(v_a, v_b); //edge dividng cell
	Edge* e_a1 = nullptr; 
	Edge* e_a2 = nullptr;
	Edge* e_b1 = nullptr;
	Edge* e_b2 = nullptr;
	
	//create edges so that e_x1 is always before e_x2 in vertex/edge ordering	
	auto updateAB = [this](Cell* const c_x, Edge* const e_x, Vertex* const v_x, Edge*& e_x1, Edge*& e_x2)
	{
		const std::vector<Vertex*>& c_x_vertices = c_x->vertices();
		std::vector<Vertex*>::const_iterator it_v1 = std::find(c_x_vertices.begin(), c_x_vertices.end(), e_x->v1());
		std::vector<Vertex*>::const_iterator it_v2 = std::find(c_x_vertices.begin(), c_x_vertices.end(), e_x->v2());
		
		int i_vx1 = std::distance(c_x_vertices.begin(), it_v1); int i_vx2 = std::distance(c_x_vertices.begin(), it_v2);
		if (((i_vx2+2) % c_x_vertices.size()) == i_vx1) { std::swap(it_v1, it_v2); std::swap(i_vx1, i_vx2); }

		e_x1 = T->createEdge(*it_v1, v_x);
		e_x2 = T->createEdge(*it_v2, v_x);
		if (c_x != nullptr)
		{
			const std::vector<Edge*>& c_x_edges = c_x->edges();
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
	
	
	auto isEdge = [this](Edge* e, Vertex* v_1, Vertex* v_2) { return ( (v_1 == e->v1() && v_2 == e->v2()) || (v_1 == e->v2() && v_2 == e->v1()) ); };	
	//find vertices of the two new cells in rotaional order
	std::vector<Vertex*> c_p_vertices, c_q_vertices;
	std::vector<Edge*> c_p_edges, c_q_edges;
	bool f_q = !(i_va < i_vb); size_t n = vertices_.size();
	for (int i = 0; i < n; i++) //neeed to fix
	{
		if (!f_q)
		{ 
			c_p_vertices.push_back(vertices_[i]);
			if (isEdge(e_a, vertices_[i], vertices_[(i+1)%n]))
			{
				c_p_vertices.push_back(v_a);
				c_p_vertices.push_back(v_b);
				f_q = true;
			}
		}	
		else
		{
			c_q_vertices.push_back(vertices_[i]);
			if (isEdge(e_b, vertices_[i], vertices_[(i+1)%n]))
			{
				c_q_vertices.push_back(v_b);
				c_q_vertices.push_back(v_a);		
				f_q = false;
			}
		}
	}
	
	//update edges of two new cells in rotational order
	size_t np = c_p_vertices.size();
	for (int i = 0; i < np; i++)
	{
		if (isEdge(e_a2, c_p_vertices[i], c_p_vertices[(i+1)%np])) c_p_edges.push_back(e_a2);
		else if (isEdge(e_b1, c_p_vertices[i], c_p_vertices[(i+1)%np])) c_p_edges.push_back(e_b1);
		else if (isEdge(e_new, c_p_vertices[i], c_p_vertices[(i+1)%np])) c_p_edges.push_back(e_new);
		else for (Edge* e : edges_) if (isEdge(e, c_p_vertices[i], c_p_vertices[(i+1)%np])) c_p_edges.push_back(e);
	}
	size_t nq = c_q_vertices.size();
	for (int i = 0; i < nq; i++)
	{
		if (isEdge(e_b2, c_q_vertices[i], c_q_vertices[(i+1)%nq])) c_q_edges.push_back(e_b2);
		else if (isEdge(e_a1, c_q_vertices[i], c_q_vertices[(i+1)%nq])) c_q_edges.push_back(e_a1);
		else if (isEdge(e_new, c_q_vertices[i], c_q_vertices[(i+1)%nq])) c_q_edges.push_back(e_new);
		else for (Edge* e : edges_) if (isEdge(e, c_q_vertices[i], c_q_vertices[(i+1)%nq])) c_q_edges.push_back(e);
	}
	
	Cell* c_p = T->createCell(c_p_vertices, c_p_edges);
	Cell* c_q = T->createCell(c_q_vertices, c_q_edges);
	
	std::vector<Cell*> neighbour_copy = neighbours_;
	T->destroyCell(this);
	c_p->findNeighbours(); c_q->findNeighbours();
	for (Cell* c : neighbour_copy) c->findNeighbours();
	
	for (Vertex* v : c_p_vertices) v->orderCellContacts();
	for (Vertex* v : c_q_vertices) v->orderCellContacts();
	std::cout << "cell divided\n";
}


void Cell::calcR_0()
{
    double x_sum = 0;
    double y_sum = 0;
    for (Vertex* v : vertices_) 
    {
        x_sum += v->r().x();
        y_sum += v->r().y();
    }
    r_0_ = Point(x_sum/vertices_.size(), y_sum/vertices_.size());
}

void Cell::calcA()
{
	A_ = 0;
	size_t n = vertices_.size();
	for (int i = 0; i < n; i++)
	{
		Point r_i = vertices_[i]->r();
		Point r_j = vertices_[(i+1)%n]->r();
		A_ += r_i.x()*r_j.y() - r_j.x()*r_i.y();
	} A_*= 0.5;
}

void Cell::calcL()
{
	L_ = 0; 
	for (Edge* e : edges_) L_ += e->l(); //edge lengths must already be calculated
}

void Cell::calcT_A() { T_A_ = param::K_a*(A_-param::A_0); }

void Cell::calcG()
{
	G[0]=0;	G[1]=0; G[2]=0;
	calcR_0();
	double x_0 = r_0_.x(); double y_0 = r_0_.y();
	for (Vertex* v : vertices_)
	{
		double x_v = v->r().x();
		double y_v = v->r().y();
		G[0]+=(x_v-x_0)*(x_v-x_0);	
		G[1]+=(x_v-x_0)*(y_v-y_0);
		G[2]+=(y_v-y_0)*(y_v-y_0);
	}
	
	double f = 1.0/vertices_.size();
	G[0]*=f; G[1]*=f; G[2]*=f;
	
	lambda = 0.5*( G[0]+G[2] + std::sqrt( (G[0]+G[2])*(G[0]+G[2]) - 4*(G[0]*G[2]-G[1]*G[1]) ) );
	n_ = Vec(1, (lambda-G[0])/G[1]) / std::sqrt( 1 + ((lambda-G[0])/G[1])*((lambda-G[0])/G[1]) );
	
	Z_ = 0.5*(G[0]-G[2]);
	X_ = G[1];
}

void Cell::calcm()
{
	size_t n = neighbours_.size();
	double w = 0;
	for (int i = 0; i < n; i++) w += T->D_angle(neighbours_[i], neighbours_[(i+1)%n]);
	m_ = w*boost::math::constants::one_div_two_pi <double>();
}



bool Cell::valid()
{
	if (vertices_.size() != edges_.size()) return false;
	size_t n = vertices_.size();
	for (int i = 0; i < n; i++) if ((edges_[i]->hasVertex(vertices_[i]) && edges_[i]->hasVertex(vertices_[(i+1)%n])) == false) return false;
	return true;
}


/*void Cell::outputVertices() const 
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
}*/
