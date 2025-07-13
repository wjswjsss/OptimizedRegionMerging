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
void analyzeLAS();
void centerLAS(const char *inFile, const char *outFile);

int main()
{
	analyzeLAS();
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
	std::cout << RED_TEXT << "U r running the Rabbani Criteria based Region Merging Algo." << RESET_COLOR
			  << std::endl;

	const char *inputPath = "../../../pointCloudData/output.las";
	const char *outputPath = "../../../pointCloudData/outputCentered.las";

	try
	{
		std::cout << "Centering point cloud...\n";
		centerLAS(inputPath, outputPath);
		std::cout << "Point cloud centered and written to " << outputPath << "\n";
	}
	catch (std::exception &e)
	{
		std::cerr << "Error: " << e.what() << "\n";
		exit(1);
	}

	// std::string path =      "E:\\01paper\\mypaper7_clustering_acceleation\\data\\pointcloud\\tree\\plot_1_annotated.las";
	std::string path = "../../../pointCloudData/outputCentered.las"; // plot_1_annotated_tree_160k
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
	RunMode run_mode = RunMode::Pruning_Cycle;
	printMode(run_mode);
	rabbani->run(scale_parameter, search_radius, region_size, run_mode);

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

// Helper to compute percentile from a sorted vector
double percentile(const std::vector<double> &data, double p)
{
	if (data.empty())
		return 0.0;
	double idx = p / 100.0 * (data.size() - 1);
	size_t lo = static_cast<size_t>(std::floor(idx));
	size_t hi = static_cast<size_t>(std::ceil(idx));
	if (hi == lo)
		return data[lo];
	double weight = idx - lo;
	return data[lo] * (1.0 - weight) + data[hi] * weight;
}

void analyzeLAS()
{
	// 1. hard-coded paths
	const char *lasFilePath = "../../../pointCloudData/output.las";

	LASreadOpener lasreadopener;
	lasreadopener.set_file_name(lasFilePath);
	lasreadopener.parse(0, nullptr);
	LASreader *reader = lasreadopener.open();
	if (!reader)
	{
		std::cerr << "ERROR: could not open LAS file: "
				  << lasFilePath << std::endl;
		return;
	}

	std::uint64_t count = 0;
	double sumX = 0, sumY = 0, sumZ = 0;
	double sumX2 = 0, sumY2 = 0, sumZ2 = 0;

	std::vector<double> xs, ys, zs;
	xs.reserve(1000000);
	ys.reserve(1000000);
	zs.reserve(1000000);

	while (reader->read_point())
	{
		double x = reader->get_x();
		double y = reader->get_y();
		double z = reader->get_z();

		sumX += x;
		sumY += y;
		sumZ += z;
		sumX2 += x * x;
		sumY2 += y * y;
		sumZ2 += z * z;
		xs.push_back(x);
		ys.push_back(y);
		zs.push_back(z);
		++count;
	}

	reader->close();
	delete reader;

	if (count == 0)
	{
		std::cerr << "No points found in the LAS file.\n";
		return;
	}

	// Sort for min/max/percentiles
	std::sort(xs.begin(), xs.end());
	std::sort(ys.begin(), ys.end());
	std::sort(zs.begin(), zs.end());

	double meanX = sumX / count;
	double meanY = sumY / count;
	double meanZ = sumZ / count;

	double varX = sumX2 / count - meanX * meanX;
	double varY = sumY2 / count - meanY * meanY;
	double varZ = sumZ2 / count - meanZ * meanZ;

	double stdX = std::sqrt(varX);
	double stdY = std::sqrt(varY);
	double stdZ = std::sqrt(varZ);

	// Min / Max
	double minX = xs.front(), maxX = xs.back();
	double minY = ys.front(), maxY = ys.back();
	double minZ = zs.front(), maxZ = zs.back();

	// Percentiles to compute
	const double percentiles[] = {1, 5, 25, 50, 75, 95, 99};

	// Output
	std::cout << "Point count: " << count << "\n\n";

	std::cout << "Axis Statistics:\n";
	auto printAxis = [&](const char axis,
						 double minv, double maxv,
						 double meanv, double stdv,
						 const std::vector<double> &data)
	{
		std::cout << axis << "-axis:\n"
				  << "  Min = " << minv << "\n"
				  << "  Max = " << maxv << "\n"
				  << "  Mean = " << meanv << "\n"
				  << "  Std  = " << stdv << "\n";

		std::cout << "  Percentiles:\n";
		for (double p : percentiles)
		{
			double val = percentile(data, p);
			std::cout << "    " << p << "% → " << val << "\n";
		}
		std::cout << "\n";
	};

	printAxis('X', minX, maxX, meanX, stdX, xs);
	printAxis('Y', minY, maxY, meanY, stdY, ys);
	printAxis('Z', minZ, maxZ, meanZ, stdZ, zs);
}

/// Centers all points in `inFile` around the median of each axis
/// and writes the result to `outFile`.
void centerLAS(const char *inFile, const char *outFile)
{
	// --- FIRST PASS: collect coordinates and header ---
	LASreadOpener readerOpener;
	readerOpener.set_file_name(inFile);
	readerOpener.parse(0, nullptr);
	LASreader *reader1 = readerOpener.open();
	if (!reader1)
	{
		throw std::runtime_error(std::string("Cannot open input LAS: ") + inFile);
	}

	std::vector<double> xs, ys, zs;
	xs.reserve(1000000);
	ys.reserve(1000000);
	zs.reserve(1000000);

	while (reader1->read_point())
	{
		xs.push_back(reader1->get_x());
		ys.push_back(reader1->get_y());
		zs.push_back(reader1->get_z());
	}

	LASheader header = reader1->header; // copy for writer
	reader1->close();
	delete reader1;

	if (xs.empty())
	{
		throw std::runtime_error("Input LAS contains no points.");
	}

	// sort to find median
	auto find_median = [&](std::vector<double> &v)
	{
		std::sort(v.begin(), v.end());
		return v[v.size() / 2];
	};
	double medX = find_median(xs);
	double medY = find_median(ys);
	double medZ = find_median(zs);

	// --- SECOND PASS: subtract medians and write out ---
	LASreader *reader2 = readerOpener.open(); // reopen same file
	LASwriteOpener writerOpener;
	writerOpener.set_file_name(outFile);
	LASwriter *writer = writerOpener.open(&header);
	if (!writer)
	{
		reader2->close();
		delete reader2;
		throw std::runtime_error(std::string("Cannot open output LAS: ") + outFile);
	}

	while (reader2->read_point())
	{
		LASpoint p = reader2->point;
		double cx = p.get_x() - medX;
		double cy = p.get_y() - medY;
		double cz = p.get_z() - medZ;

		p.set_X(cx);
		p.set_Y(cy);
		p.set_Z(cz);
		writer->write_point(&p);
	}

	writer->close();
	delete writer;
	reader2->close();
	delete reader2;
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
