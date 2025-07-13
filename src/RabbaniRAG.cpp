// RabbaniRAG.cpp
#include "RabbaniRAG.h"
#include "RabbaniVertex.h"
#include "RAGCallbacks.h"

#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d.h>

#include <iostream>
#include <vector>

using R_GetFeat = RAGCallbacks::GetFeatures<RabbaniFeatures, RabbaniEdge, RabbaniVertex>;
using R_Update = RAGCallbacks::UpdateFeatures<RabbaniFeatures, RabbaniEdge, RabbaniVertex>;
using R_Merge = RAGCallbacks::MergingCriteria<RabbaniFeatures, RabbaniEdge, RabbaniVertex>;

RabbaniRAG::RabbaniRAG(Point_Cloud *pc,
                       double search_radius,
                       R_Update update_features,
                       R_Merge mc)
    : RegionAdjacencyGraph(pc->count, -1, update_features, mc),
      pc(pc), search_radius(search_radius)
{
    this->cloud = pc->pcl_cloud;
    this->vertice_num = cloud->size();
    this->features = new RabbaniFeatures[vertice_num];
    this->vertice = new RabbaniVertex[vertice_num];
    this->normal_cloud = nullptr;

    // 1) Estimate normals once in constructor:
    pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
    ne.setInputCloud(cloud);
    auto tree = std::make_shared<pcl::search::KdTree<pcl::PointXYZ>>();
    ne.setSearchMethod(tree);
    ne.setRadiusSearch(search_radius);
    normal_cloud.reset(new pcl::PointCloud<pcl::Normal>);
    ne.compute(*normal_cloud);

    // // TODO: wait for implementation.
    this->connect_vertex_and_edge();

    // build_adjacency(nullptr);
}

RabbaniRAG::RabbaniRAG(RabbaniVertex *root,
                       int size,
                       double search_radius,
                       R_Update update_features,
                       R_Merge mc)
    : RegionAdjacencyGraph(size, -1, update_features, mc),
      search_radius(search_radius)
{
    vertice_num = size;
    features = new RabbaniFeatures[size];
    vertice = new RabbaniVertex[size];
    // ? Should come back later. (when running multiple rounds)
    // build_adjacency(root);
}

RabbaniRAG::~RabbaniRAG()
{
    delete[] vertice;
    delete[] edges;
    delete[] features;
}

// ! Calling parent's method
// RabbaniVertex *RabbaniRAG::run(double scale_parameter,
//                                RunMode mode,
//                                int region_size)
// {
//     // single-round merge: build adjacency, then invoke base run
//     // build_adjacency(nullptr);
//     return static_cast<RabbaniVertex *>(RegionAdjacencyGraph::run(scale_parameter, mode, region_size));
// }

void RabbaniRAG::output(const char *output_path)
{
    std::cout << "NOT IMPLEMENTED: output to " << output_path << std::endl;
    // // TODO: implement output to file (e.g., color labels)
}

void RabbaniRAG::connect_vertex_and_edge()
{
    this->initial_vertice();
    this->edge_num = this->kd_tree_radius_consturct_topo(this->search_radius);
    this->initial_edges(this->edge_num);
}

void RabbaniRAG::initial_vertice()
{
    for (int i = 0; i < vertice_num; ++i)
    {
        // extract point and normal
        const auto &pt = cloud->points[i];
        const auto &pn = normal_cloud->points[i];

        // build Eigen normal, normalize with fallback
        Eigen::Vector3d n(pn.normal_x, pn.normal_y, pn.normal_z);
        if (n.norm() > 1e-6)
            n.normalize();
        else
            n = Eigen::Vector3d::UnitZ();

        // initialize vertex and feature
        RabbaniVertex *v = &vertice[i];
        v->ID = i;
        v->set_position_dimention(3);
        v->positions[0] = pt.x;
        v->positions[1] = pt.y;
        v->positions[2] = pt.z;
        v->normal = n;

        features[i].set_feature(v, pt.x, pt.y, pt.z, n);
        v->ptr_features = reinterpret_cast<RootFeatures *>(&features[i]);

        // chain the region list
        if (i + 1 < vertice_num)
            v->next_in_region = &vertice[i + 1];
        else
            v->next_in_region = nullptr;
    }
}

int RabbaniRAG::kd_tree_radius_consturct_topo(double r)
{
    // 1) build the FLANN tree on the point cloud
    pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
    kdtree.setInputCloud(cloud);

    int edge_idx = 0;
    // 2) for each point i
    for (int i = 0; i < vertice_num; ++i)
    {
        std::vector<int> nbrs;
        std::vector<float> dists;
        const auto &pt = cloud->points[i];
        RabbaniVertex *v = &vertice[i];

        // 3) radius search around pt
        if (kdtree.radiusSearch(pt, r, nbrs, dists) > 0)
        {
            // skip nbrs[0] (it’s always i itself)
            for (size_t j = 1; j < nbrs.size(); ++j)
            {
                int ni = nbrs[j];
                // only add each pair once
                if (ni > i)
                {
                    RabbaniVertex *w = &vertice[ni];
                    v->edges.push_back(edge_idx);
                    w->edges.push_back(edge_idx);
                    ++edge_idx;
                }
            }

            // 4) if truly isolated (only itself in nbrs), connect to 2-nd NN
            if (nbrs.size() == 1)
            {
                std::vector<int> knn_i(2);
                std::vector<float> knn_d(2);
                if (kdtree.nearestKSearch(pt, 2, knn_i, knn_d) > 1)
                {
                    int ni = knn_i[1]; // the closest other point
                    RabbaniVertex *w = &vertice[ni];
                    v->edges.push_back(edge_idx);
                    w->edges.push_back(edge_idx);
                    ++edge_idx;
                }
            }
        }
    }

    // 5) return how many edges you inserted (used next by init_edges)
    return edge_idx;
}

void RabbaniRAG::initial_edges(int edge_count)
{
    // 1) allocate the edge array and record count
    this->edges = new RabbaniEdge[edge_count];
    this->edge_num = edge_count;

    // 2) walk each vertex’s adjacency list to hook up left/right pointers
    for (int i = 0; i < vertice_num; ++i)
    {
        RabbaniVertex *v = &vertice[i];
        for (int eid : v->edges)
        {
            RabbaniEdge *e = &edges[eid];
            e->ID = eid;
            if (e->left_vertex == nullptr)
            {
                e->left_vertex = v;
            }
            else if (e->right_vertex == nullptr)
            {
                e->right_vertex = v;
            }
            else
            {
                throw std::runtime_error("Repeated edge error in init_edges");
            }
        }
    }

    // 3) compute each edge’s weight (angle difference) via the merging_criteria callback
    for (int eid = 0; eid < edge_num; ++eid)
    {
        RabbaniEdge *e = &edges[eid];
        auto *A = static_cast<RabbaniFeatures *>(e->left_vertex->ptr_features);
        auto *B = static_cast<RabbaniFeatures *>(e->right_vertex->ptr_features);

        if (e->left_vertex && e->right_vertex)
        {
            double angle = merging_criteria(A, B, e);
            e->weight = angle; // store the angle as weight
            // store it in your edge’s field
            e->angle_diff = angle;
        }
        else
        {
            throw std::runtime_error("Edge endpoints not set before computing weight");
        }
    }
}
