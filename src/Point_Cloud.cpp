#include "Point_Cloud.h"
#include <vector>
#include <chrono>
#include "lasreader.hpp"
#include "laswriter.hpp"
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h> // ���磬ʹ�������˲�
#define RED_TEXT "\033[31m"
#define GREEN_TEXT "\033[32m"
#define RESET_COLOR "\033[0m"
Point_Cloud::Point_Cloud(const char *path)
{
	// auto start_time = std::chrono::high_resolution_clock::now();
	this->file_path = path;
	this->lasLoad.set_file_name(this->file_path);

	this->pcl_cloud = this->get_pcl_points();
	// this->labels = this->get_unique_labels();
	this->count = this->pcl_cloud->size();
	std::cout << "point num: " << this->count << std::endl;
}

Point_Cloud::~Point_Cloud() {}

LASheader Point_Cloud::get_header()
{
	this->lasLoad.reset();
	LASreader *reader = this->lasLoad.open();
	return reader->header;
}

void Point_Cloud::save_as_las_file(const char *output_path, const char *attribute_name, double attribute_value)
{

	this->lasLoad.reset();
	LASreader *reader = this->lasLoad.open();
	int attribute_index = reader->header.get_attribute_index(attribute_name);

	LASwriteOpener laswriteopener;
	laswriteopener.set_file_name(output_path); // ��������ļ���
	LASwriter *writer = laswriteopener.open(&reader->header);
	if (writer == 0)
	{
		fprintf(stderr, "ERROR: could not open LAS writer.\n");
	}

	while (reader->read_point())
	{
		LASpoint *point = &(reader->point);
		float value = point->get_attribute_as_float(attribute_index);
		// int int_value = int(value);
		if (attribute_value != value)
		{
			writer->write_point(point);
		}
	}
	reader->close();
	writer->close();
}

pcl::PointCloud<pcl::PointXYZ>::Ptr Point_Cloud::get_pcl_points()
{
	this->lasLoad.reset();
	LASreader *reader = this->lasLoad.open();
	// this->header = &reader->header;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	if (reader == nullptr)
		return nullptr;
	int num = 0;
	while (reader->read_point())
	{
		pcl::PointXYZ point;
		point.x = reader->point.get_x() - 599066.00;
		point.y = reader->point.get_y() - 6655080.00;
		point.z = reader->point.get_z();
		cloud->points.push_back(point);

		num++;
		// if (num == 100000)
		//{
		//	break;
		// }
	}
	// std::cout << "PCL Num: "<< num << std::endl;
	reader->close();
	return cloud;
}

vector<int> Point_Cloud::get_unique_labels()
{
	this->lasLoad.reset();
	LASreader *las_reader = this->lasLoad.open();
	vector<int> tree_labels;
	int Num = 0;
	int attribute_index = las_reader->header.get_attribute_index("treeID");
	while (las_reader->read_point())
	{
		float tree_id = las_reader->point.get_attribute_as_float(attribute_index);
		int tree_label = int(tree_id);
		// is 'label' existed
		if (std::find(tree_labels.begin(), tree_labels.end(), tree_label) == tree_labels.end())
		{
			// insert label
			tree_labels.push_back(tree_label);
		}
		Num++;
	}
	// std::cout << "tree_labels: ";
	// for (int label : tree_labels) {
	//	std::cout << label << " ";
	// }
	// std::cout << std::endl;
	las_reader->close();
	return tree_labels;
}

LASreader *Point_Cloud::get_las_reader()
{
	this->lasLoad.reset();
	return this->lasLoad.open();
}