#include "Vertex.h"
#include "Edge.h"

Vertex::Vertex() : edges()
{
	this->ID = -1;
	this->connection_degree = 0;
	this->radius = 0;
	// this->node_type = 1;//State repalced the node_type
	this->parent = nullptr;
	this->next = nullptr;
	this->next_in_region = nullptr;
	this->last = this;
	this->color = -1;
	this->is_visited = false;
	this->state = State::Activate;
	this->size = 1;
	this->direct = nullptr;
	this->ptr_features = nullptr;
	this->positions = nullptr;
	this->position_dimension = -1;
	this->bit_sign = 0;
	this->bits_num = 0;
}

Vertex::Vertex(int id, int connection_degree, double radius) : Vertex()
{
	this->ID = id;
	this->connection_degree = connection_degree;
	this->radius = radius;
}

Vertex::Vertex(int ID, RootFeatures *features) : Vertex()
{
	this->ID = ID;
	// this->ptr_features = (MRSFeatures*)features;
	this->ptr_features = (RootFeatures *)features;
}

// ���ƹ��캯��ʵ��
Vertex::Vertex(const Vertex &other) : ID(other.ID),
									  position_dimension(other.position_dimension),
									  radius(other.radius),
									  connection_degree(other.connection_degree),
									  color(other.color),
									  size(other.size),
									  state(other.state),
									  parent(other.parent),
									  next(other.next),
									  next_in_region(other.next_in_region),
									  last(other.last),
									  is_visited(other.is_visited),
									  edges(other.edges),
									  direct(other.direct),
									  ptr_features(other.ptr_features),
									  bit_sign(other.bit_sign),
									  bits_num(other.bits_num)
{
	if (other.positions != nullptr && position_dimension > 0)
	{
		delete[] this->positions;
		this->positions = new double[position_dimension];
		std::copy(other.positions, other.positions + position_dimension, positions);
	}
	else
	{
		positions = nullptr;
	}
}

// ��ֵ�����ʵ��
Vertex &Vertex::operator=(const Vertex &other)
{
	if (this != &other)
	{
		// �ͷ�������Դ
		delete[] positions;

		// ���Ƴ�Ա����
		ID = other.ID;
		position_dimension = other.position_dimension;
		radius = other.radius;
		connection_degree = other.connection_degree;
		color = other.color;
		size = other.size;
		state = other.state;
		parent = other.parent;
		next = other.next;
		next_in_region = other.next_in_region;
		last = other.last;
		is_visited = other.is_visited;
		edges = other.edges;
		direct = other.direct;
		ptr_features = other.ptr_features;
		bit_sign = other.bit_sign;
		bits_num = other.bits_num;

		// ��� positions ����
		if (other.positions != nullptr && position_dimension > 0)
		{
			delete[] positions;
			positions = new double[position_dimension];
			std::copy(other.positions, other.positions + position_dimension, positions);
		}
		else
		{
			positions = nullptr;
		}
	}
	return *this;
}

Vertex::~Vertex()
{
	this->edges.clear();
	this->edges.~vector();
	// this->ID;
	// if(ptr_features !=nullptr) delete ptr_features;
	if (this->positions != nullptr)
		delete[] this->positions;
}

void Vertex::push(Edge *ptr_e)
{
	this->edges.push_back(ptr_e->ID);
}

void Vertex::push_range(int *first, int *last)
{
	this->edges.insert(this->edges.end(), first, last);
}

void Vertex::add_next_vertex(Vertex *ptr_next)
{
	this->last->next = ptr_next;
	ptr_next->parent = this->last;
	this->last = ptr_next->last;
	this->size += ptr_next->size;
	if (this->is_built_tree)
	{
		this->merged_vertice.push_back(ptr_next->ID);
	}
}

void Vertex::remove_edge(int edge_id)
{
	this->edges.erase(std::remove(this->edges.begin(), this->edges.end(), edge_id), this->edges.end());
	this->edges.shrink_to_fit();
}

void Vertex::change_direct(Edge *new_direct, NNG *nng)
{
	// int ttt = 0;
	if (new_direct == nullptr)
	{
		throw std::runtime_error("Cannot change the direction of vertex into a null ptr.)");
	}

	if (this->direct == nullptr)
	{
		this->direct = new_direct;
		this->direct->shortest_arc_num++;
		if (this->direct->shortest_arc_num == 2)
		{

			nng->push((Edge *)this->direct);
		}
		return;
	}
	if (this->direct == new_direct)
	{
		return;
	}

	if (this->direct->shortest_arc_num == 2)
		nng->remove((Edge *)this->direct);
	this->direct->shortest_arc_num--;
	this->direct = new_direct;
	this->direct->shortest_arc_num++;
	if (this->direct->shortest_arc_num == 2)
	{
		nng->push((Edge *)this->direct);
	}
}

void Vertex::set_position_dimention(int dimention)
{
	if (this->position_dimension > 0)
	{
		delete[] this->positions; // �ͷ�֮ǰ���ڴ�
		this->positions = new double[dimention];
		this->position_dimension = dimention;
	}
	else
	{
		this->position_dimension = dimention;
		this->positions = new double[dimention];
	}
}

void Vertex::append_sign(bool sign)
{
	// ͨ������������
	this->bit_sign = this->bit_sign << 1;
	this->bit_sign += sign;
	this->bits_num++;
}

int Vertex::get_sign()
{
	return this->bit_sign;
}

void Vertex::reset()
{
	this->bits_num = 0;
	this->bit_sign = 0;
	this->direct = nullptr;
	this->is_visited = false;
	this->state = State::Activate;
}

bool Vertex::operator<(const Vertex &other) const
{
	return this->edges.size() < other.edges.size();
}