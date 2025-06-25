#pragma once

#include "Vertex.h"
#include "BaseFeatures.h"
#include <Eigen/Core> // for Eigen::Vector3d

class Edge_Cloud;
class NNG;
class RabbaniFeatures;

class RabbaniVertex : public Vertex
{
public:
    // --- our new attribute
    Eigen::Vector3d normal = Eigen::Vector3d::Zero();

    // constructors
    RabbaniVertex();
    RabbaniVertex(int id);
    RabbaniVertex(int id, int k);
    RabbaniVertex(int id, int k, RabbaniFeatures *features);

    ~RabbaniVertex();

    // show the linked vertices by ID
    void show_vertices();
    virtual int GetID() override;

    // if you still need Org_ID semantics:
    int Org_ID;
};
