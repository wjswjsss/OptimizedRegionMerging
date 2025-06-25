// RabbaniRAG.h
#pragma once

#include "Point_Cloud.h"
#include "RabbaniVertex.h"
#include "RabbaniEdge.h"
#include "BaseFeatures.h"
#include "RegionAdjacencyGraph.h"
#include <pcl/point_cloud.h>
#include <pcl/search/kdtree.h>
#include "RAGCallbacks.h"

using R_GetFeat = RAGCallbacks::GetFeatures<RabbaniFeatures, RabbaniEdge, RabbaniVertex>;
using R_Update = RAGCallbacks::UpdateFeatures<RabbaniFeatures, RabbaniEdge, RabbaniVertex>;
using R_Merge = RAGCallbacks::MergingCriteria<RabbaniFeatures, RabbaniEdge, RabbaniVertex>;

class RabbaniRAG : public RegionAdjacencyGraph<RabbaniVertex, RabbaniEdge, RabbaniFeatures>
{
public:
    RabbaniRAG(Point_Cloud *pc,
               double search_radius,
               R_Update update_features,
               R_Merge mc);

    RabbaniRAG(RabbaniVertex *root,
               int size,
               double search_radius,
               R_Update update_features,
               R_Merge mc);

    ~RabbaniRAG();

    // perform a single round of region merging
    // RabbaniVertex *run(double scale_parameter,
    //                    RunMode mode,
    //                    int region_size);

    void output(const char *output_path) override;

private:
    Point_Cloud *pc;
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud;
    pcl::PointCloud<pcl::Normal>::Ptr normal_cloud;
    RabbaniFeatures *features;
    double search_radius;

    // + -----------------------------------------------+
    // | for initial the vertices, features, and edges. |
    // + -----------------------------------------------+
    void connect_vertex_and_edge() override;   // high level call.
    void initial_vertice();                    // Phase 1
    int kd_tree_radius_consturct_topo(double); // Phase 2
    void initial_edges(int);                   // Phase 3
    // ------------------------------------------------ +
};
