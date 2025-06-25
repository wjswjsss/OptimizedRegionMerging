// RabbaniEdge.cpp
#include "RabbaniEdge.h"
#include <algorithm> // for std::clamp
#include <cmath>     // for std::acos

RabbaniEdge::RabbaniEdge()
    : Edge(), angle_diff(0.0)
{
}

RabbaniEdge::RabbaniEdge(int ID, RabbaniVertex *left_vertex, RabbaniVertex *right_vertex)
    : Edge(ID, left_vertex, right_vertex)
{
    // compute angle between normals once and store it
    // const auto &n1 = static_cast<RabbaniVertex *>(left_vertex)->normal;
    // const auto &n2 = static_cast<RabbaniVertex *>(right_vertex)->normal;
    // double cosang = n1.dot(n2) / (n1.norm() * n2.norm());
    // cosang = std::clamp(cosang, -1.0, 1.0);
    // angle_diff = std::acos(cosang);
}

RabbaniEdge::~RabbaniEdge()
{
}

double RabbaniEdge::GetWeight()
{
    return angle_diff;
}