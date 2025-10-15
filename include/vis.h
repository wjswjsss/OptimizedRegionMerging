// Function to visualize Rabbani regions as one colored point cloud
#pragma once

// your own vertex/region definition
#include "RabbaniVertex.h"

// PCL core and types
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

// PCL visualizer
#include <pcl/visualization/pcl_visualizer.h>

// C++ random number generation
#include <random>

// thread sleep
#include <thread>
#include <chrono>

void visualizeRabbaniRegions(RabbaniVertex *region_root,
                             pcl::PointCloud<pcl::PointXYZ>::Ptr originalCloud)
{
    // 1. Prepare one big colored cloud
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr colored_cloud(new pcl::PointCloud<pcl::PointXYZRGB>);
    colored_cloud->header = originalCloud->header;
    colored_cloud->is_dense = originalCloud->is_dense;

    // 2. Walk through each region in the linked list
    std::mt19937 rng{std::random_device{}()};
    std::uniform_int_distribution<int> dist(0, 255);
    size_t regionCount = 0;
    size_t totalPoints = 0;

    for (RabbaniVertex *region = region_root; region; region = static_cast<RabbaniVertex *>(region->next_in_region))
    {
        // Pick a random color for this region
        uint8_t r = static_cast<uint8_t>(dist(rng));
        uint8_t g = static_cast<uint8_t>(dist(rng));
        uint8_t b = static_cast<uint8_t>(dist(rng));

        // Traverse all points in this region
        for (RabbaniVertex *v = region; v; v = static_cast<RabbaniVertex *>(v->next))
        {
            const auto &src = originalCloud->points[v->ID];
            pcl::PointXYZRGB p;
            p.x = src.x;
            p.y = src.y;
            p.z = src.z;
            p.r = r;
            p.g = g;
            p.b = b;
            colored_cloud->points.push_back(p);
            ++totalPoints;
        }

        ++regionCount;
    }

    colored_cloud->width = static_cast<uint32_t>(colored_cloud->points.size());
    colored_cloud->height = 1;

    std::cout << "Regions: " << regionCount
              << " | Total points: " << totalPoints << std::endl;

    // 3. Visualize just the colored cloud
    pcl::visualization::PCLVisualizer viewer("Region Viewer");
    viewer.setBackgroundColor(0, 0, 0);
    viewer.addPointCloud<pcl::PointXYZRGB>(colored_cloud, "colored_cloud");
    viewer.setPointCloudRenderingProperties(
        pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "colored_cloud");

    // 4. Spin until the window closes
    while (!viewer.wasStopped())
    {
        viewer.spinOnce(100);
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
    }
}