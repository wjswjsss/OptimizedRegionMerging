// VirtualSeedsCpp.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
// #pragma warning(disable:4996)
// #pragma execution_character_set("utf-8")

#include <iostream>
#include <chrono>
#include <cmath>

// #include <liblas/liblas.hpp>
#include <fstream>
#include <iostream>
#include <string>

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

#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdint>
#include "cpl_conv.h"
#define RED_TEXT "\033[31m"
#define GREEN_TEXT "\033[32m"
#define RESET_COLOR "\033[0m"

using namespace std;
void OpenLAS();
void OpenPCD();
void OpenLASRabbani();
std::string centerPointCloud(const std::string &inputPath);

int main()
{
	OpenLASRabbani();
	// OpenPCD();

	std::cout << "END" << std::endl;
}
void printMode(RunMode runMode)
{
	switch (runMode)
	{
	case RunMode::None:
		std::cout << "RunMode: None" << std::endl;
		break;
	case RunMode::Region:
		std::cout << "RunMode: Region" << std::endl;
		break;
	case RunMode::Pruning:
		std::cout << "RunMode: Pruning" << std::endl;
		break;
	case RunMode::Cycle:
		std::cout << "RunMode: Cycle" << std::endl;
		break;
	case RunMode::Region_Pruning:
		std::cout << "RunMode: Region_Pruning" << std::endl;
		break;
	case RunMode::Region_Cycle:
		std::cout << "RunMode: Region_Cycle" << std::endl;
		break;
	case RunMode::Pruning_Cycle:
		std::cout << "RunMode: Pruning_Cycle" << std::endl;
		break;
	case RunMode::Region_Pruning_Cycle:
		std::cout << "RunMode: Region_Pruning_Cycle" << std::endl;
		break;
	case RunMode::Region_Pruning_Cycle_MultiCores:
		std::cout << "RunMode: Region_Pruning_Cycle_MultiCores" << std::endl;
		break;
	case RunMode::Region_Pruning_Cycle_GPU:
		std::cout << "RunMode: Region_Pruning_Cycle_GPU" << std::endl;
		break;
	default:
		std::cout << "RunMode: Unknown" << std::endl;
		break;
	}
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
	std::cout << RED_TEXT
			  << "U r running the Rabbani Criteria based Region Merging Algo."
			  << RESET_COLOR << std::endl;

	const char *rawPath = "../../../pointCloudData/output.las";
	std::string path = rawPath;
	path = centerPointCloud(path); // Center the point cloud before processing

	// now path either points to the original or the new centered file
	std::cout << "[Using point cloud file] " << path << std::endl;

	// rest of your code unchanged, but using 'path'
	auto start_time = std::chrono::high_resolution_clock::now();
	int region_size = 2500;
	double scale_parameter = 4;
	double search_radius = 20;
	RunMode run_mode = RunMode::Region_Pruning_Cycle_MultiCores;

	std::cout << RED_TEXT << "[Call] RabbaniRM3D::run" << RESET_COLOR << std::endl;
	printMode(run_mode);

	RabbaniRM3D *rabbani = new RabbaniRM3D(path.c_str());
	rabbani->run(scale_parameter, search_radius, region_size, run_mode);
	delete rabbani;

	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
	if (duration < 10)
		std::cout << "Total time: " << GREEN_TEXT << duration << RESET_COLOR << "ms\n";
	else
		std::cout << "Total time: " << RED_TEXT << duration << RESET_COLOR << "ms\n";
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

// new Functions
//------------------------------------------------------------------------------
// Compute centroid, shift all points by subtracting the mean, write out new file
//------------------------------------------------------------------------------
///
/// Compute the centroid of the input .las, subtract it from every point,
/// write out to "<inputBase>_centered.las", log the mean, and return the new path.
///
std::string centerPointCloud(const std::string &inputPath)
{
	std::cout << "[INFO] centerPointCloud() called with inputPath = “"
			  << inputPath << "”\n";

	// 1) OPEN READER
	LASreadOpener ropr;
	ropr.set_file_name(inputPath.c_str());
	std::cout << "[INFO] Attempting to open LASreader…\n";
	LASreader *reader = ropr.open();
	if (!reader)
	{
		std::cerr << "[ERROR] LASreadOpener::open() returned nullptr for “"
				  << inputPath << "”\n";
		throw std::runtime_error("Failed to open LAS reader for " + inputPath);
	}
	std::cout << "[INFO] LASreader opened at ptr=" << reader << "\n";

	// 1a) LOG HEADER METADATA (with correct field names)
	{
		auto &H = reader->header;
		std::cout << "[INFO] Header:\n"
				  << "       Version:             " << int(H.version_major)
				  << "." << int(H.version_minor) << "\n"
				  << "       Point format:        " << int(H.point_data_format) << "\n"
				  << "       Point count:         " << H.number_of_point_records << "\n"
				  << "       Point record length: " << H.point_data_record_length << "\n"
				  << "       Bounds X:            [" << H.min_x << ", " << H.max_x << "]\n"
				  << "       Bounds Y:            [" << H.min_y << ", " << H.max_y << "]\n"
				  << "       Bounds Z:            [" << H.min_z << ", " << H.max_z << "]\n";
	}

	// 2) ACCUMULATE SUM & COUNT
	std::uint64_t count = 0;
	double sumX = 0, sumY = 0, sumZ = 0;
	while (reader->read_point())
	{
		sumX += reader->get_x();
		sumY += reader->get_y();
		sumZ += reader->get_z();
		++count;
	}
	if (count == 0)
	{
		reader->close();
		delete reader;
		throw std::runtime_error("No points in " + inputPath);
	}

	// 3) COMPUTE & LOG MEANS
	double meanX = sumX / count;
	double meanY = sumY / count;
	double meanZ = sumZ / count;
	std::cout << "[INFO] Computed means over " << count << " points → "
			  << "meanX=" << meanX << "  meanY=" << meanY
			  << "  meanZ=" << meanZ << "\n";

	// 4) REOPEN READER FOR WRITING PASS
	std::cout << "[INFO] Reopening reader for write pass…\n";
	reader->close();
	delete reader;
	// reader = ropr.open();
	// if (!reader)
	// {
	// 	std::cerr << "[ERROR] Re-open failed!\n";
	// 	throw std::runtime_error("Failed to reopen LAS reader for write pass");
	// }
	LASreadOpener ropr2;
	ropr2.set_file_name(inputPath.c_str());
	std::cout << "[INFO] Re-opening with fresh opener…\n";
	reader = ropr2.open();
	if (!reader)
	{
		std::cerr << "[ERROR] Re-open still failed!\n";
		throw std::runtime_error("Failed to reopen LAS reader for write pass");
	}
	std::cout << "[INFO] Reader reopened at ptr=" << reader << "\n";

	// 5) PREPARE WRITER
	std::string outPath = inputPath.substr(0, inputPath.find_last_of('.')) + "_centered.las";
	std::cout << "[INFO] Preparing to write centered file to “"
			  << outPath << "”\n";

	// 5a) Quick file-creatability test
	{
		std::ofstream ofs_test(outPath, std::ios::binary);
		if (!ofs_test)
		{
			std::cerr << "[ERROR] Cannot create “" << outPath
					  << "” — check directory & permissions\n";
		}
		else
		{
			std::cout << "[INFO] File creation test passed\n";
		}
	}

	// 5b) Actual LASwriter open
	std::cout << "[INFO] Calling LASwriteOpener::open(&reader->header)…\n";
	LASwriteOpener wopr;
	wopr.set_file_name(outPath.c_str());
	LASwriter *writer = nullptr;
	try
	{
		writer = wopr.open(&reader->header);
	}
	catch (const std::exception &e)
	{
		std::cerr << "[EXCEPTION] wopr.open() threw: " << e.what() << "\n";
	}
	if (!writer)
	{
		std::cerr << "[ERROR] wopr.open() returned nullptr\n";
		reader->close();
		delete reader;
		throw std::runtime_error("Failed to open LAS writer for " + outPath);
	}
	std::cout << "[INFO] LASwriter opened at ptr=" << writer << "\n";

	// 6) WRITE SHIFTED POINTS
	std::cout << "[INFO] Writing centered points…\n";
	while (reader->read_point())
	{
		double x = reader->get_x(),
			   y = reader->get_y(),
			   z = reader->get_z();
		reader->point.set_x(x - meanX);
		reader->point.set_y(y - meanY);
		reader->point.set_z(z - meanZ);
		writer->write_point(&reader->point);
	}
	std::cout << "[INFO] Finished writing.  Total points written: "
			  << writer->npoints << "\n";

	// 7) CLEAN UP
	std::cout << "[INFO] Cleaning up reader & writer…\n";
	writer->close();
	delete writer;
	reader->close();
	delete reader;
	std::cout << "[INFO] centerPointCloud() done.  Output = “"
			  << outPath << "”\n";

	return outPath;
}