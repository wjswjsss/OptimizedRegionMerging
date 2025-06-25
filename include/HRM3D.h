#pragma once
#include "Point_Cloud.h"
#include "RAG_Cloud.h"
#include "RAGCallbacks.h"

using H_GetFeat = RAGCallbacks::GetFeatures<PointCloudFeatures, Edge_Cloud, Vertex_Cloud>;
using H_Update = RAGCallbacks::UpdateFeatures<PointCloudFeatures, Edge_Cloud, Vertex_Cloud>;
using H_Merge = RAGCallbacks::MergingCriteria<PointCloudFeatures, Edge_Cloud, Vertex_Cloud>;

class HRM3D
{
public:
	H_GetFeat get_features_func;
	H_Update update_features_func;
	H_Merge merging_criteria;
	HRM3D(const char *);
	~HRM3D();

private:
	Point_Cloud *cloud;
	RAG_Cloud *rag;

public:
	void run(double scale_parameter, double alpha = 1.0, double search_radius = 0.1, int region_size = 2500, RunMode mode = RunMode::None);
	void assignment_org_ID(Vertex_Cloud *root);
	void assignment_value_to_org_vertex(Vertex_Cloud *root);
	int check_num(Vertex_Cloud *root);
	void output(const char *);
};
