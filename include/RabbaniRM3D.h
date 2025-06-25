// RabbaniRM3D.h
#pragma once

#include "Point_Cloud.h"
// #include "RAG_Cloud.h"
#include "RabbaniVertex.h"
#include "RabbaniEdge.h"
#include "RabbaniRAG.h"
#include "RAGCallbacks.h"

using R_GetFeat = RAGCallbacks::GetFeatures<RabbaniFeatures, RabbaniEdge, RabbaniVertex>;
using R_Update = RAGCallbacks::UpdateFeatures<RabbaniFeatures, RabbaniEdge, RabbaniVertex>;
using R_Merge = RAGCallbacks::MergingCriteria<RabbaniFeatures, RabbaniEdge, RabbaniVertex>;

class RabbaniRM3D
{
public:
    // feature and merging function types expect RabbaniEdge*
    // PointCloudNamespace::GetFeatures get_features_func;
    // PointCloudNamespace::UpdateFeatures update_features_func;
    // PointCloudNamespace::MergingCriteria merging_criteria;
    R_Update update_features_func;
    R_Merge merging_criteria;

    RabbaniRM3D(const char *cloud_file);
    ~RabbaniRM3D();

private:
    Point_Cloud *cloud;
    RabbaniRAG *rag;

public:
    void run(double scale_parameter,
             double search_radius = 0.1,
             int region_size = 2500,
             RunMode mode = RunMode::None);

    void assignment_org_ID(RabbaniVertex *root);
    void assignment_value_to_org_vertex(RabbaniVertex *root);
    int check_num(RabbaniVertex *root);
    void output(const char *out_file);
};