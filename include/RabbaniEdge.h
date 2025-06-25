// RabbaniEdge.h
#pragma once

#include "Edge.h"
#include "RabbaniVertex.h" // custom vertex with normal
#include <Eigen/Core>      // for Eigen::Vector3d, std::clamp

class RabbaniEdge : public Edge
{
public:
    // constructors
    RabbaniEdge();
    RabbaniEdge(int ID, RabbaniVertex *left_vertex, RabbaniVertex *right_vertex);
    ~RabbaniEdge();

    // override to return stored angle difference
    double GetWeight() override;

    // store the computed angle (in radians)
    double angle_diff = 0.0;
};