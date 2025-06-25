// VirtualSeedsCpp.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// #pragma warning(disable:4996)
// #pragma execution_character_set("utf-8")

#include <iostream>
#include <chrono>
#include <cmath>

// gdal
// #include"gdal_priv.h"

// my cpp head
// #include"Image.h"
// #include"RAG_IMG.h"
// #include"HSWO.h"
// #include"MRS.h"
#include "Point_Cloud.h"
#include "HRM3D.h"
#include "RabbaniRM3D.h"
// #include"GaiaDR3.h"
// #include"BloomingTree.h"

// laslib
#include "lasreader.hpp"
#include "laswriter.hpp"
// pcl
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/passthrough.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>
#include <pcl/segmentation/supervoxel_clustering.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/surface/convex_hull.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <chrono>
#include <thread>

#include "cpl_conv.h"
#define RED_TEXT "\033[31m"
#define GREEN_TEXT "\033[32m"
#define RESET_COLOR "\033[0m"

using namespace std;
void OpenLAS();
void OpenPCD();
void OpenLASRabbani();

int main()
{

	OpenLASRabbani();
	// OpenPCD();

	std::cout << "END" << std::endl;
}

// laslib
void OpenLAS()
{
	// std::string path =      "E:\\01paper\\mypaper7_clustering_acceleation\\data\\pointcloud\\tree\\plot_1_annotated.las";
	std::string path = "E:\\01paper\\mypaper7_clustering_acceleation\\data\\pointcloud\\tree\\plot_1_annotated_tree.las"; // plot_1_annotated_tree_160k
	// std::string path = "E:\\01paper\\mypaper7_clustering_acceleation\\data\\pointcloud\\tree\\plot_1_annotated_tree_160k.las";//

	std::string save_path = "E:\\01paper\\mypaper7_clustering_acceleation\\data\\pointcloud\\tree\\plot_1_annotated_tree_Seg.las";

	auto start_time = std::chrono::high_resolution_clock::now();

	// Point_Cloud* pr_C = new Point_Cloud(path.c_str());
	// pr_C->save_as_las_file(save_path.c_str(), "treeID", 0.0);

	int region_size = 10000;
	HRM3D *hrm3d = new HRM3D(path.c_str());
	double scale_parameter = 4.0;
	double alpha = 0.3;
	double search_radius = 0.1;
	hrm3d->run(scale_parameter, alpha, search_radius, region_size, RunMode::Region_Pruning_Cycle_MultiCores);
	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
	long long t = duration.count();
	// hrm3d->output(save_path.c_str());
	delete hrm3d;
	if (t < 10)
		std::cout << "Total time:" << GREEN_TEXT << t << RESET_COLOR << "ms" << std::endl;
	else
		std::cout << "Total time:" << RED_TEXT << t << RESET_COLOR << "ms" << std::endl;
}

// + -------------------------------------- +
// | Rabbani criteria driven Region Merging |
// + -------------------------------------- +
void OpenLASRabbani()
{
	std::cout << RED_TEXT << "U r running the Rabbani Criteria based Region Merging Algo." << RESET_COLOR
			  << std::endl;
	// std::string path =      "E:\\01paper\\mypaper7_clustering_acceleation\\data\\pointcloud\\tree\\plot_1_annotated.las";
	std::string path = "../../../pointCloudData/output.las"; // plot_1_annotated_tree_160k
	// std::string path = "E:\\01paper\\mypaper7_clustering_acceleation\\data\\pointcloud\\tree\\plot_1_annotated_tree_160k.las";//

	std::string save_path = "./";

	auto start_time = std::chrono::high_resolution_clock::now();

	// Point_Cloud* pr_C = new Point_Cloud(path.c_str());
	// pr_C->save_as_las_file(save_path.c_str(), "treeID", 0.0);

	int region_size = 200;

	// + -------------------------------------------------------- +
	// | Instantiate Rabbani criteria driven Region Merging Class |
	// + -------------------------------------------------------- +
	RabbaniRM3D *rabbani = new RabbaniRM3D(path.c_str());

	double scale_parameter = 44.9;
	double search_radius = 100;

	// + --------------------------------------------------------------- +
	// | Calling Rabbani criteria driven Region Merging Class run method |
	// + --------------------------------------------------------------- +
	std::cout << RED_TEXT << "[Call] RabbaniRM3D::run" << RESET_COLOR << std::endl;

	rabbani->run(scale_parameter, search_radius, region_size, RunMode::Region_Pruning_Cycle_MultiCores);

	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
	long long t = duration.count();
	// hrm3d->output(save_path.c_str());
	delete rabbani;
	if (t < 10)
		std::cout << "Total time:" << GREEN_TEXT << t << RESET_COLOR << "ms" << std::endl;
	else
		std::cout << "Total time:" << RED_TEXT << t << RESET_COLOR << "ms" << std::endl;
}

/// test pcl
void OpenPCD()
{
	// string pcd_file = "E:\\01paper\\mypaper7_clustering_acceleation\\data\\pointcloud\\table_scene_lms400.pcd";
	// string pcd_file = "E:\\01paper\\mypaper7_clustering_acceleation\\data\\pointcloud\\Case1_flip0_ss0.pcd";
	string pcd_file = "E:\\01paper\\mypaper7_clustering_acceleation\\data\\pointcloud\\FORinstance_dataset\\NIBIO_pcd\\plot_1_annotated_tree.pcd";
	const char *pcd_file_char = pcd_file.c_str();
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	if (pcl::io::loadPCDFile<pcl::PointXYZ>(pcd_file_char, *cloud) == -1)
	{
		PCL_ERROR("Couldn't read file your_file.pcd \n");
	}

	pcl::PCLPointCloud2::Ptr cloud_blob(new pcl::PCLPointCloud2);

	pcl::toPCLPointCloud2(*cloud, *cloud_blob);
	int intensity_field_index = pcl::getFieldIndex(*cloud_blob, "Intensity");
	int treeID = pcl::getFieldIndex(*cloud_blob, "treeID");
	// std::cout << "Fields: " << fields_list << std::endl;
	std::string fields_list = pcl::getFieldsList(*cloud);

	// int num_fields = cloud_blob->fields.size();

	std::cout << "fileds: " << fields_list << std::endl;

	// std::cout << "Loaded "
	//     << cloud->width * cloud->height
	//     << " data points from your_file.pcd"
	//     << std::endl;
	// std::cout << cloud->points.size() << std::endl;
	for (size_t i = 0; i < 10; ++i)
	{
		std::cout << "    " << cloud->points[i].x
				  << " " << cloud->points[i].y
				  << " " << cloud->points[i].z << std::endl;
	}
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧:
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
