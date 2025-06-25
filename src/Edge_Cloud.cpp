#include "Edge_Cloud.h"
#include "Vertex_Cloud.h"
Edge_Cloud::Edge_Cloud():Edge()
{
	this->E_Sij = 0.0;
}

Edge_Cloud::Edge_Cloud(int ID, Vertex_Cloud* left_vertex, Vertex_Cloud* right_vertex)
	: Edge(ID, left_vertex, right_vertex)
{
	this->E_Sij = 0.0;
}

Edge_Cloud::Edge_Cloud(int ID, double weight, Vertex_Cloud* left_vertex, Vertex_Cloud* right_vertex)
	: Edge(ID, weight, left_vertex, right_vertex)
{
	this->E_Sij = 0.0;
}

Edge_Cloud::~Edge_Cloud()
{
}

double Edge_Cloud::GetWeight()
{
	return this->weight;
}
