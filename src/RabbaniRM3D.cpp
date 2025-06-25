// RabbaniRM3D.cpp

#include <thread>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/point_types.h>
#include <pcl/PointIndices.h>
#include <pcl/console/print.h>

#include "RabbaniRM3D.h"
#include "RabbaniRAG.h"
#include <iostream>
#include <chrono>
#include "BaseFeatures.h"
void RabbaniUpdateFeatures(RabbaniFeatures *A,
                           RabbaniFeatures *B,
                           RabbaniEdge * /*e*/);

double RabbaniAngleCriterion(RabbaniFeatures *A,
                             RabbaniFeatures *B,
                             RabbaniEdge *e);

// Placeholder free functions (to be implemented):
// double normal_angle_difference(PointCloudFeatures* A, PointCloudFeatures* B, RabbaniEdge* e);
// void   merge_region_normals(PointCloudFeatures* A, PointCloudFeatures* B, RabbaniEdge* e);

// void RabbaniUpdateFeatures(RabbaniFeatures *, RabbaniFeatures *, RabbaniEdge *);
// double calculate_total_E_Sij(RabbaniFeatures *, RabbaniFeatures *);
// double RabbaniAngleCriterion(RabbaniFeatures *, RabbaniFeatures *, RabbaniEdge *);

RabbaniRM3D::RabbaniRM3D(const char *cloud_file)
{
    // this->get_features_func = /* your feature extractor */;
    this->update_features_func = RabbaniUpdateFeatures;
    this->merging_criteria = RabbaniAngleCriterion;
    this->cloud = new Point_Cloud(cloud_file);
}

RabbaniRM3D::~RabbaniRM3D()
{
    delete this->cloud;
}

void RabbaniRM3D::run(double scale_parameter, double search_radius, int region_size, RunMode mode)
{
    // single-round region merging (multi-round to be added later)
    this->rag = new RabbaniRAG(
        this->cloud,
        search_radius,
        this->update_features_func,
        this->merging_criteria);

    std::cout << RED_TEXT << "Running single-round RabbaniRM3D merge..." << RESET_COLOR << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();

    // perform one merge round
    RabbaniVertex *root = static_cast<RabbaniVertex *>(
        this->rag->run(scale_parameter, mode, region_size));

    auto t1 = std::chrono::high_resolution_clock::now();
    long long duration = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();

    std::cout << "Merge time: " << duration << " ms" << std::endl;

    // initialize original IDs
    this->assignment_org_ID(root);

    // optionally, you could update values against pre-existing root here
    // this->assignment_value_to_org_vertex(root);

    int sz = this->check_num(root); // check how many clusters Survived.
    std::cout << "Resulting region count: " << sz << std::endl;

    std::cout << "Total run time: " << duration << " ms" << std::endl;

    // Start the visualization
    // 1. Grab the head of the region list and the original PCL cloud pointer:
    RabbaniVertex *region_root = root;
    auto originalCloud = this->cloud->pcl_cloud; // assume pcl_cloud is a pcl::PointCloud<PointXYZ>::Ptr

    // 2. Create the PCLVisualizer in‐line (no thread):
    pcl::visualization::PCLVisualizer::Ptr viewer{
        new pcl::visualization::PCLVisualizer("Region Visualization")};
    viewer->setBackgroundColor(0, 0, 0);
    viewer->initCameraParameters();
    viewer->addCoordinateSystem(1.0);

    // 3. Define a simple color palette (expand if you expect >10 regions):
    std::vector<std::tuple<uint8_t, uint8_t, uint8_t>> colors = {
        {255, 0, 0},   // red
        {0, 255, 0},   // green
        {0, 0, 255},   // blue
        {255, 255, 0}, // yellow
        {255, 0, 255}, // magenta
        {0, 255, 255}, // cyan
        {128, 128, 0}, // olive
        {128, 0, 128}, // purple
        {0, 128, 128}, // teal
        {128, 64, 0}   // brown
    };

    size_t totalPoints = 0;
    int regionIdx = 0;
    for (RabbaniVertex *r = region_root;
         r != nullptr;
         r = static_cast<RabbaniVertex *>(r->next_in_region))
    {
        // Build a PointCloud<PointXYZRGB> for this region:
        pcl::PointCloud<pcl::PointXYZRGB>::Ptr regionCloud{new pcl::PointCloud<pcl::PointXYZRGB>()};
        regionCloud->header = originalCloud->header;
        regionCloud->is_dense = originalCloud->is_dense;

        auto [r_col, g_col, b_col] = colors[regionIdx % colors.size()];

        // Traverse all vertices (points) in this region via the 'next' link:
        for (RabbaniVertex *v = r; v != nullptr; v = static_cast<RabbaniVertex *>(v->next))
        {
            int ptIdx = v->ID;
            const auto &basePt = originalCloud->points[ptIdx];

            pcl::PointXYZRGB p;
            p.x = basePt.x;
            p.y = basePt.y;
            p.z = basePt.z;
            p.r = r_col;
            p.g = g_col;
            p.b = b_col;
            regionCloud->points.push_back(p);
        }

        regionCloud->width = static_cast<uint32_t>(regionCloud->points.size());
        regionCloud->height = 1;

        totalPoints += regionCloud->points.size();

        std::string cloudName = "region_" + std::to_string(regionIdx);
        viewer->addPointCloud<pcl::PointXYZRGB>(regionCloud, cloudName);
        viewer->setPointCloudRenderingProperties(
            pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, cloudName);

        ++regionIdx;
    }

    std::cout << "Total points to visualize: " << totalPoints << std::endl;

    // 4. Now spin in the same thread until the window closes:
    while (!viewer->wasStopped())
    {
        viewer->spinOnce(100);
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

void RabbaniRM3D::assignment_org_ID(RabbaniVertex *root)
{
    for (RabbaniVertex *cur = root; cur; cur = static_cast<RabbaniVertex *>(cur->next_in_region))
    {
        cur->Org_ID = cur->ID;
    }
}

void RabbaniRM3D::assignment_value_to_org_vertex(RabbaniVertex *root)
{
    for (RabbaniVertex *cur = root; cur; cur = static_cast<RabbaniVertex *>(cur->next_in_region))
    {
        for (RabbaniVertex *nxt = static_cast<RabbaniVertex *>(cur->next); nxt; nxt = static_cast<RabbaniVertex *>(nxt->next))
        {
            RabbaniVertex *org_v = &this->rag->vertice[cur->Org_ID];
            RabbaniVertex *org_n = &this->rag->vertice[nxt->Org_ID];
            org_v->add_next_vertex(org_n);
            org_n->state = State::Freeze;
        }
    }
}

int RabbaniRM3D::check_num(RabbaniVertex *root)
{
    int count = 0;
    for (RabbaniVertex *cur = root; cur; cur = static_cast<RabbaniVertex *>(cur->next_in_region))
    {
        ++count;
    }
    return count;
}

void RabbaniRM3D::output(const char *out_file)
{
    std::cout << "Outputting results to " << out_file << std::endl;
    this->rag->output(out_file);
}

double RabbaniAngleCriterion(RabbaniFeatures *A,
                             RabbaniFeatures *B,
                             RabbaniEdge *e)
{
    const Eigen::Vector3d &n1 = A->normal;
    const Eigen::Vector3d &n2 = B->normal;
    double cosang = n1.dot(n2) / (n1.norm() * n2.norm());
    cosang = std::clamp(cosang, -1.0, 1.0);
    double angle = std::acos(cosang);
    e->angle_diff = angle;
    return angle;
}

void RabbaniUpdateFeatures(RabbaniFeatures *A,
                           RabbaniFeatures *B,
                           RabbaniEdge * /*e*/)
{
    // 1) new total size
    int total = A->size + B->size;

    // 2) weighted normal fusion
    Eigen::Vector3d nm = A->normal * double(A->size) + B->normal * double(B->size);
    if (nm.norm() > 1e-6)
        nm.normalize();
    A->normal = nm;

    // 3) weighted centroid fusion
    A->x_bar = (A->x_bar * A->size + B->x_bar * B->size) / double(total);
    A->y_bar = (A->y_bar * A->size + B->y_bar * B->size) / double(total);
    A->z_bar = (A->z_bar * A->size + B->z_bar * B->size) / double(total);

    // 4) update size
    A->size = total;

    // 5) highest_v remains the same representative
}