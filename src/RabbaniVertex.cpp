#include "RabbaniVertex.h"
#include "BaseFeatures.h"
#include <iostream>

class RabbaniFeatures;

RabbaniVertex::RabbaniVertex()
    : Vertex(), normal(Eigen::Vector3d::Zero()), Org_ID(-1)
{
    this->connection_degree = 0;
}

RabbaniVertex::RabbaniVertex(int id)
    : Vertex(), normal(Eigen::Vector3d::Zero())
{
    this->ID = id;
    this->connection_degree = 0;
    this->Org_ID = -1;
}

RabbaniVertex::RabbaniVertex(int id, int k)
    : RabbaniVertex(id)
{
    this->connection_degree = k;
}

RabbaniVertex::RabbaniVertex(int id, int k, RabbaniFeatures *features)
    : Vertex(id, (RootFeatures *)features), normal(Eigen::Vector3d::Zero())
{
    this->connection_degree = k;
    this->Org_ID = -1;
}

RabbaniVertex::~RabbaniVertex() {}

void RabbaniVertex::show_vertices()
{
    RabbaniVertex *cur = this;
    while (cur)
    {
        std::cout << cur->ID << " ";
        cur = static_cast<RabbaniVertex *>(cur->next);
    }
    std::cout << "\n";
}

int RabbaniVertex::GetID()
{
    return this->ID;
}
