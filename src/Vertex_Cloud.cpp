#include "Vertex_Cloud.h"

Vertex_Cloud::Vertex_Cloud():Vertex()
{
	this->connection_degree = 0;
	this->Org_ID = -1;
}

Vertex_Cloud::Vertex_Cloud(int id):Vertex()
{
	this->ID = id;
}

Vertex_Cloud::Vertex_Cloud(int id,int k): Vertex_Cloud(id)
{
	this->connection_degree = k;
}

Vertex_Cloud::Vertex_Cloud(int id,int k, PointCloudFeatures* features): 
	Vertex(id, (RootFeatures*)features)
{
	//this->highest_v = this;
	this->connection_degree = k;
	this->Org_ID = -1;
}

Vertex_Cloud::~Vertex_Cloud()
{
}

void Vertex_Cloud::show_vertice()
{
	Vertex_Cloud* SIL = this;
	while (SIL != nullptr)
	{
		std::cout << SIL->ID << " ";
		SIL = (Vertex_Cloud*)SIL->next;
	}
	std::cout << std::endl;
}

int Vertex_Cloud::GetID()
{
	return this->ID;
}
