#pragma once
#include <iostream>
#include <vector>

// laslib
#include "lasreader.hpp"
#include "laswriter.hpp"
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h> // ���磬ʹ�������˲�

using namespace std;
class Point_Cloud
{
public: // ����
	Point_Cloud(const char *);

	~Point_Cloud();
	LASreader *get_las_reader();

public: // ��Ա
	std::string pro_info;
	uint32_t count;
	vector<int> labels;
	pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_cloud;
	// LASheader header;
	LASheader get_header();
	void save_as_las_file(const char *, const char *, double);

private:
	pcl::PointCloud<pcl::PointXYZ>::Ptr get_pcl_points();

private:
	LASreadOpener lasLoad;
	const char *file_path = nullptr;
	vector<int> get_unique_labels();
};
