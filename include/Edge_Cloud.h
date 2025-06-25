#pragma once
#include"Edge.h"

class Vertex_Cloud;

class Edge_Cloud : public Edge
{
public:
	Edge_Cloud();
	Edge_Cloud(int ID, Vertex_Cloud* left_vertex, Vertex_Cloud* right_vertex);
	Edge_Cloud(int ID, double weight, Vertex_Cloud* left_vertex, Vertex_Cloud* right_vertex);
	~Edge_Cloud();
	double GetWeight() override;

public:
	double E_Sij = 0.0;

};

