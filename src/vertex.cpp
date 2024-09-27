#include "vertex.h"
#include "tissue.h"


Vertex::Vertex(Tissue* T, Point r) : T(T), r_(r), force_(Vec(0,0)) 
{ 
	not_boundary_cell = 1; 
	cell_contacts_ordered.reserve(8);
}
Vertex::Vertex() = default;

bool Vertex::operator==(const Vertex& other) const { return r_ == other.r_; }

bool Vertex::onBoundaryCell()
{
	for (Cell* c : cell_contacts_) 
	{
		if (c->onBoundary()) { not_boundary_cell = 0; return true; }
	}
	not_boundary_cell = 1;
	return false;
}

const Point& Vertex::r() const { return r_; }
const double Vertex::m() const { return m_; }
const std::unordered_set<Cell*>& Vertex::cellContacts() const { return cell_contacts_; }
const std::unordered_set<Edge*>& Vertex::edgeContacts() const { return edge_contacts_; }

void Vertex::addCellContact(Cell* c) { cell_contacts_.insert(c); }
void Vertex::removeCellContact(Cell* c) { cell_contacts_.erase(c); }

void Vertex::addEdgeContact(Edge* e) { edge_contacts_.insert(e); }
void Vertex::removeEdgeContact(Edge* e) 
{
	edge_contacts_.erase(e);
	if (edge_contacts_.size() == 0) T->destroyVertex(this);
}


Vec Vertex::calcSurfaceForce()
{
	Vec f_A(0,0);
	for (Cell* c : cell_contacts_) 
	{
		const std::vector<Vertex*>& c_vertices = c->vertices();
		int n = c_vertices.size(); int S = c->S();
		std::vector<Vertex*>::const_iterator it = std::find(c_vertices.begin(), c_vertices.end(), this);
		int j = std::distance(c_vertices.begin(), it);
		
		double dAdx = S*0.5*(c_vertices[(j+1)%n]->r().y() - c_vertices[(j-1+n)%n]->r().y());
		double dAdy = S*0.5*(c_vertices[(j-1+n)%n]->r().x() - c_vertices[(j+1)%n]->r().x());
		f_A -= c->T_A()*Vec(dAdx, dAdy);
	}
	return f_A;
}

Vec Vertex::calcLineForce()
{
	Vec f_L(0,0);
	for (Edge* e : edge_contacts_)
	{
		Vertex* v; //other vertex in edge
		(this == e->v1()) ? v = e->v2() : v = e->v1();
		
		double x_diff = r_.x() - v->r().x();
		double y_diff = r_.y() - v->r().y();
		double dldx = x_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		double dldy = y_diff/std::sqrt(x_diff*x_diff+y_diff*y_diff);
		f_L -= e->T_l()*Vec(dldx, dldy);
	}
	return f_L;
}

void Vertex::calcForce() { force_ += calcSurfaceForce()+calcLineForce(); }
void Vertex::applyForce() { r_ += not_boundary_cell*param::a*param::dt*force_ + 100*(1-not_boundary_cell)*param::a*param::dt*Vec(-r_.y(),r_.x()); }
void Vertex::shearForce() { force_ = Vec(-r_.y(),r_.x()); } //anticlockwise shear


void Vertex::orderCellContacts()
{
	std::vector<std::pair<Cell*, double>> contacts;
	for (Cell* c : cell_contacts_) 
	{
		c->calcR_0();
		Vec vec = c->r_0()-r_;
		double theta = std::atan2(vec.y(), vec.x());
		contacts.push_back({c,theta});
	}
	
	std::sort(contacts.begin(), contacts.end(), 
		[](const std::pair<Cell*, double>& c1, const std::pair<Cell*, double>& c2) { return c1.second < c2.second; });
	cell_contacts_ordered = contacts;
}

void Vertex::T1split()
{
	if (cell_contacts_.size() != 4 || edge_contacts_.size() != 4) return;
	for (Cell* c : cell_contacts_) if (c->onBoundary()) return; 
	
	//affected cells, a,b change vertex p,q gets new edge
	orderCellContacts();
	Cell* const c_a = cell_contacts_ordered[0].first; Cell* const c_b = cell_contacts_ordered[2].first;
	Cell* const c_p = cell_contacts_ordered[1].first; Cell* const c_q = cell_contacts_ordered[3].first;

	const std::vector<Vertex*>& c_p_vertices = c_p->vertices();
	auto it_id = std::find(c_p_vertices.begin(), c_p_vertices.end(), this);
	int i = std::distance(c_p_vertices.begin(), it_id);
	
	//new vertices
	Vertex* v_a = nullptr; 
	Vertex* v_b = nullptr;
	auto updateAB = [this](Cell* const c_x, Vertex*& v_x)
	{
		//find position of and create vertex
		c_x->calcR_0();
		Vec vec_x = c_x->r_0() - r_;
		Point x = r_ + 0.1*vec_x;
		v_x = T->createVertex(x);

		//attatch relevent edges to vertex
		for (Edge* e : c_x->edges()) e->swapVertex(this, v_x);
		if (!(c_x->valid())) c_x->rotateVertices();
	};
	updateAB(c_a, v_a); updateAB(c_b, v_b);

	//edge that vertex is split into
	Edge* const e_new = T->createEdge(v_a, v_b);

	//update edges for cells p, q
	auto updateEdges = [this, c_a, c_b, v_a, v_b, e_new](Cell* const c_x)
	{
		const std::vector<Vertex*>& c_x_vertices = c_x->vertices();
		const std::vector<Edge*>& c_x_edges = c_x->edges();
		size_t n = c_x_edges.size();
		for (int i = 0; i < n; i++)
		{
			if (c_x_edges[i]->hasVertex(v_a) && c_x_edges[(i+1)%n]->hasVertex(v_b))
			{
				T->cellNewEdge(c_x, e_new, i+1);
				break;
			}
			else if (c_x_edges[i]->hasVertex(v_b) && c_x_edges[(i+1)%n]->hasVertex(v_a))
			{
				T->cellNewEdge(c_x, e_new, i+1);
				break;
			}
		}		
	};
	updateEdges(c_p); updateEdges(c_q);
	
	auto updateVertices = [this, c_a, c_b, v_a, v_b, e_new](Cell* const c_x, bool before)
	{
		const std::vector<Vertex*>& c_x_vertices = c_x->vertices();
		std::vector<Vertex*>::const_iterator it_va = std::find(c_x_vertices.begin(), c_x_vertices.end(), v_a);
		
		if (it_va != c_x_vertices.end())
		{
			int i = std::distance(c_x_vertices.begin(), it_va);
			if (before) T->cellNewVertex(c_x, v_b, i);
			else T->cellNewVertex(c_x, v_b, i+1);
			while (!(c_x->valid())) c_x->rotateVertices();
		}
	};
	updateVertices(c_p, true); updateVertices(c_q, false);

	T->destroyVertex(this);
	v_a->orderCellContacts(); v_b->orderCellContacts();
	for (Cell* c : v_a->cellContacts()) c->findNeighbours(); 
	for (Cell* c : v_b->cellContacts()) c->findNeighbours();
	//std::cout << "T1 split\n";
}

void Vertex::calcm()
{
	//cell order already known
	size_t n = cell_contacts_.size();
	double w = 0;
	for (int i = 0; i < n; i++) w += T->D_angle(cell_contacts_ordered[i].first, cell_contacts_ordered[(i+1)%n].first);
	m_ = w*boost::math::constants::one_div_two_pi <double>();
}
