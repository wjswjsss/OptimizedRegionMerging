#pragma once
#include "BaseFeatures.h"
#include <iostream>
#include <map>
#include <vector>
#include "BaseFeatures.h"

#include "Vertex.h"
class Edge_Cloud;
class NNG;

//
#include <Eigen/Core> // or PCL’s Normal type
//

class Vertex_Cloud : public Vertex
{
public:
	Vertex_Cloud();
	Vertex_Cloud(int id);
	Vertex_Cloud(int id, int k);
	Vertex_Cloud(int id, int k, PointCloudFeatures *features);
	~Vertex_Cloud();

public:
	//
	// Eigen::Vector3d normal = Eigen::Vector3d::Zero(); // <<< REGION normal here
	//

	void show_vertice();
	virtual int GetID() override;
	// Vertex_Cloud* highest_v = nullptr;
	int Org_ID;
};
