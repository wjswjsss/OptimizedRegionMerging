#include "Edge.h"
#include"Vertex.h"
Edge::Edge()
{
	this->ID = -1;
	this->shortest_arc_num = 0;
	this->length = 1;
	this->weight = -1.0;
	//this->type = 1;
	this->left_vertex = nullptr;
	this->right_vertex = nullptr;
	this->state = State::Activate;
	this->next_edge_in_region = nullptr;
	this->is_checked = false;
	//this->direct = nullptr;
}

Edge::Edge(int ID, Vertex* left_vertex, Vertex* right_vertex):Edge()
{
	this->ID = ID;

	this->left_vertex = left_vertex;
	this->right_vertex = right_vertex;
	//this->state = State::Activate;
}

Edge::Edge(int ID, double weight, Vertex* left_vertex, Vertex* right_vertex): Edge(ID, left_vertex, right_vertex)
{
	this->weight = weight;
}

Edge::~Edge()
{
	//if( left_vertex !=nullptr) delete left_vertex;
	//if(right_vertex != nullptr) delete right_vertex;
	//if(next_edge_in_region != nullptr) delete next_edge_in_region;
}

Vertex* Edge::get_other_vertex(Vertex* aim_v)
{
	if (this->left_vertex == aim_v)
	{
		return this->right_vertex;
	}
	else if (this->right_vertex == aim_v)
	{
		return this->left_vertex;
	}
	else
	{
		throw std::runtime_error("Error: edge is not connected to the vertex!");
	}
}

bool Edge::compareByValue(const Edge& other) const
{
	return this->weight == other.weight;
}

bool Edge::compareByAddress(const Edge& other) const
{
	return this == &other;
}

bool Edge::compareByEnds(const Edge& other) const
{
	bool ends_match = this->left_vertex == other.left_vertex && this->right_vertex == other.right_vertex ||
		this->left_vertex == other.right_vertex && this->right_vertex == other.left_vertex;
	return ends_match;
}

bool Edge::operator<(const Edge& other) const {
	return this->weight < other.weight;
}

bool Edge::operator>(const Edge& other) const {
	return this->weight > other.weight;
}

bool Edge::operator==(const Edge& other) const {
	return this == &other;
}

bool Edge::operator!=(const Edge& other) const {
	return this != &other;;
}