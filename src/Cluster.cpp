#include "Cluster.h"
#include <cmath>

#include <iostream>
#include <vector>
#include <thread>
#include <algorithm>
#include <limits>
#include <future>

//#include "Vertex.h"

Cluster::Cluster()
{
	this->bit_sign = -1;
	this->dataset_size = 0;
	this->interval = 1.0;
	//this->vertice_id = nullptr;
	//this->vertice = nullptr;
	this->root_vertex = nullptr;
	this->minimum_size = 1;
	this->position_dimension = 2;
	this->left_child = nullptr;
	this->right_child = nullptr;
}

Cluster::Cluster(Vertex* root_vertex) : Cluster()
{
	//this->vertice = vertice;
	//this->vertice_id = new int[this->dataset_size];//save the ids
	//for (int i = 0; i < this->dataset_size; i++)
	//{
	//	this->vertice_id[i] = ids[i];
	//}
	this->root_vertex = root_vertex;
	//this->interval = interval;
	//if (this->interval <= 0) throw std::runtime_error("The invalid interval,whcih cannot be smaller or equal to zero!");
}

Cluster::Cluster(Vertex* root_vertex, int minimum_size, int pisition_dimension, int dataset_size) : Cluster(root_vertex)
{
	this->minimum_size = minimum_size;
	this->position_dimension = position_dimension;
	this->bit_sign = root_vertex->bit_sign;
	this->dataset_size = dataset_size;
}

Cluster::~Cluster()
{
	//if (this->vertice_id != nullptr) delete[] this->vertice_id;
	if (this->left_child != nullptr) delete this->left_child;
	if (this->right_child != nullptr) delete this->right_child;
}

std::pair<Cluster*, Cluster*> Cluster::divide_cluster()
{
	std::vector<double> minimum_coordinates(position_dimension, std::numeric_limits<double>::max());
	std::vector<double> maximum_coordinates(position_dimension, std::numeric_limits<double>::min());
	std::vector<double> range_coordinates(position_dimension, 0.0);

	//search the minimun and maximun coordinates at each axis
	Vertex* v1 = this->root_vertex;
	int w = 0;
	while (v1 != nullptr)
	{
		for (int coord = 0; coord < position_dimension; coord++)
		{
			minimum_coordinates[coord] = (v1->positions[coord] < minimum_coordinates[coord]) ?
				v1->positions[coord] : minimum_coordinates[coord];
			maximum_coordinates[coord] = (v1->positions[coord] > maximum_coordinates[coord]) ?
				v1->positions[coord] : maximum_coordinates[coord];
		}
		v1 = v1->next_in_region;
		w++;
	}

	if (w != this->dataset_size)
	{
		throw std::runtime_error("Wrong dataset size!");
	}

	//calculate the range of axis and search the longest range
	double longest_range = std::numeric_limits<double>::min();
	int longest_coord = -1;
	for (int coord = 0; coord < position_dimension; coord++)
	{
		range_coordinates[coord] = maximum_coordinates[coord] - minimum_coordinates[coord];
		if (range_coordinates[coord] > longest_range)
		{
			longest_range = range_coordinates[coord];
			longest_coord = coord;
		}
	}

	//cannot conduct the operation
	if (longest_range <= 0) throw std::runtime_error("The invalid longest range!");
	//if (this->interval == 0) throw std::runtime_error("The invalid interval!");

	//statistics for the longest axis histogram
	this->interval = longest_range / 100.0;
	int num_interval = int(longest_range / this->interval + 1.0);//向上取整
	std::vector<int> histogram(num_interval, 0);

	Vertex* v2 = this->root_vertex;
	while (v2 != nullptr)
	{
		double pos = v2->positions[longest_coord];
		int pos_his = int((pos - minimum_coordinates[longest_coord]) / this->interval);//向下取整
		histogram[pos_his]++;
		v2 = v2->next_in_region;
	}

	//search for the median of the histogram
	int half_size = int(this->dataset_size * 0.5);
	int current_size = 0;
	int index = 0;
	while (current_size < half_size && index < num_interval)
	{
		current_size += histogram[index];
		index++;
	}
	double median = (index - 0.5) * this->interval + minimum_coordinates[longest_coord];

	//divide the dataset by the median
	Vertex* left_root_v = nullptr;
	Vertex* left_next_v_in_region = nullptr;
	int left_size = 0;
	Vertex* right_root_v = nullptr;
	Vertex* right_next_v_in_region = nullptr;
	int right_size = 0;

	Vertex* v3 = this->root_vertex;
	while (v3 != nullptr)
	{
		if (v3->positions[longest_coord] < median)
		{
			v3->append_sign(false);
			if (left_root_v == nullptr)
			{
				left_root_v = v3;
				left_next_v_in_region = v3;
				//left_next_v_in_region->next_in_region = nullptr;
			}
			else
			{
				left_next_v_in_region->next_in_region = v3;
				left_next_v_in_region = v3;
				//left_next_v_in_region->next_in_region = nullptr;
			}
			left_size++;
			//left_ids.push_back(v->ID);
		}
		else
		{
			v3->append_sign(true);
			if (right_root_v == nullptr)
			{
				right_root_v = v3;
				right_next_v_in_region = v3;
				//right_next_v_in_region->next_in_region = nullptr;
			}
			else
			{
				right_next_v_in_region->next_in_region = v3;
				right_next_v_in_region = v3;
				//right_next_v_in_region->next_in_region = nullptr;
			}
			right_size++;
			//right_ids.push_back(v->ID);
		}
		v3 = v3->next_in_region;
	}
	if (left_next_v_in_region != nullptr) left_next_v_in_region->next_in_region = nullptr;
	if (right_next_v_in_region != nullptr) right_next_v_in_region->next_in_region = nullptr;

	if(left_root_v !=nullptr) this->left_child = new Cluster(left_root_v, this->minimum_size, this->position_dimension,left_size);
	if(right_root_v != nullptr) this->right_child = new Cluster(right_root_v,this->minimum_size, this->position_dimension, right_size);

	return std::make_pair(this->left_child, this->right_child);
}

MinMax parallel_find_min_max(Vertex* start_vertex, int num_vertices, int position_dimension) {
	MinMax result;
	result.minimum_coordinates.resize(position_dimension, std::numeric_limits<double>::max());
	result.maximum_coordinates.resize(position_dimension, std::numeric_limits<double>::min());

	Vertex* v = start_vertex;
	for (int i = 0; i < num_vertices && v != nullptr; ++i) {
		for (int coord = 0; coord < position_dimension; coord++) {
			result.minimum_coordinates[coord] = std::min(v->positions[coord], result.minimum_coordinates[coord]);
			result.maximum_coordinates[coord] = std::max(v->positions[coord], result.maximum_coordinates[coord]);
		}
		v = v->next_in_region;
	}

	return result;
}

void parallel_histogram(Vertex* start_vertex, int num_vertices, int longest_coord, double minimum_coordinate, double interval, std::vector<int>& histogram) {
	Vertex* v = start_vertex;
	for (int i = 0; i < num_vertices && v != nullptr; ++i) {
		double pos = v->positions[longest_coord];
		int pos_his = int((pos - minimum_coordinate) / interval);
		histogram[pos_his]++;
		v = v->next_in_region;
	}
}

void parallel_divide_vertices(Vertex* start_vertex, int num_vertices, int longest_coord, double median, Vertex** left_root, Vertex** left_tail, int* left_size, Vertex** right_root, Vertex** right_tail, int* right_size) {
	Vertex* v = start_vertex;
	for (int i = 0; i < num_vertices && v != nullptr; ++i) {
		if (v->positions[longest_coord] < median) {
			v->append_sign(false);
			if (*left_root == nullptr) {
				*left_root = v;
				*left_tail = v;
			}
			else {
				(*left_tail)->next_in_region = v;
				*left_tail = v;
			}
			(*left_size)++;
		}
		else {
			v->append_sign(true);
			if (*right_root == nullptr) {
				*right_root = v;
				*right_tail = v;
			}
			else {
				(*right_tail)->next_in_region = v;
				*right_tail = v;
			}
			(*right_size)++;
		}
		v = v->next_in_region;
	}
}

std::pair<Cluster*, Cluster*> Cluster::divide_cluster_parallel() {
	int num_threads = std::thread::hardware_concurrency();

	// 并行查找最小和最大坐标
	std::vector<std::future<MinMax>> min_max_futures;
	int chunk_size = this->dataset_size / num_threads;
	Vertex* current_vertex = this->root_vertex;
	for (int i = 0; i < num_threads; ++i) {
		int start_index = i * chunk_size;
		int end_index = (i == num_threads - 1) ? this->dataset_size : (i + 1) * chunk_size;
		int num_vertices = end_index - start_index;
		min_max_futures.push_back(std::async(std::launch::async, parallel_find_min_max, current_vertex, num_vertices, this->position_dimension));
		for (int j = 0; j < num_vertices && current_vertex != nullptr; ++j) {
			current_vertex = current_vertex->next_in_region;
		}
	}

	std::vector<double> minimum_coordinates(position_dimension, std::numeric_limits<double>::max());
	std::vector<double> maximum_coordinates(position_dimension, std::numeric_limits<double>::min());
	for (auto& future : min_max_futures) {
		MinMax result = future.get();
		for (int coord = 0; coord < position_dimension; coord++) {
			minimum_coordinates[coord] = std::min(result.minimum_coordinates[coord], minimum_coordinates[coord]);
			maximum_coordinates[coord] = std::max(result.maximum_coordinates[coord], maximum_coordinates[coord]);
		}
	}

	// 计算范围并找到最长范围
	double longest_range = std::numeric_limits<double>::min();
	int longest_coord = -1;
	for (int coord = 0; coord < position_dimension; coord++) {
		double range = maximum_coordinates[coord] - minimum_coordinates[coord];
		if (range > longest_range) {
			longest_range = range;
			longest_coord = coord;
		}
	}

	// 错误处理
	if (longest_range <= 0) throw std::runtime_error("The invalid longest range!");
	//if (this->interval == 0) throw std::runtime_error("The invalid interval!");



	//statistics for the longest axis histogram
	int num_interval = int(longest_range / this->interval + 1.0);//向上取整
	std::vector<int> histogram(num_interval, 0);

	Vertex* v2 = this->root_vertex;
	while (v2 != nullptr)
	{
		double pos = v2->positions[longest_coord];
		int pos_his = int((pos - minimum_coordinates[longest_coord]) / this->interval);//向下取整
		histogram[pos_his]++;
		v2 = v2->next_in_region;
	}

	//search for the median of the histogram
	int half_size = int(this->dataset_size * 0.5);
	int current_size = 0;
	int index = 0;
	while (current_size < half_size && index < num_interval)
	{
		current_size += histogram[index];
		index++;
	}
	double median = (index - 0.5) * this->interval + minimum_coordinates[longest_coord];

	//divide the dataset by the median
	Vertex* left_root_v = nullptr;
	Vertex* left_next_v_in_region = nullptr;
	int left_size = 0;
	Vertex* right_root_v = nullptr;
	Vertex* right_next_v_in_region = nullptr;
	int right_size = 0;

	Vertex* v3 = this->root_vertex;
	while (v3 != nullptr)
	{
		if (v3->positions[longest_coord] < median)
		{
			v3->append_sign(false);
			if (left_root_v == nullptr)
			{
				left_root_v = v3;
				left_next_v_in_region = v3;
				//left_next_v_in_region->next_in_region = nullptr;
			}
			else
			{
				left_next_v_in_region->next_in_region = v3;
				left_next_v_in_region = v3;
				//left_next_v_in_region->next_in_region = nullptr;
			}
			left_size++;
			//left_ids.push_back(v->ID);
		}
		else
		{
			v3->append_sign(true);
			if (right_root_v == nullptr)
			{
				right_root_v = v3;
				right_next_v_in_region = v3;
				//right_next_v_in_region->next_in_region = nullptr;
			}
			else
			{
				right_next_v_in_region->next_in_region = v3;
				right_next_v_in_region = v3;
				//right_next_v_in_region->next_in_region = nullptr;
			}
			right_size++;
			//right_ids.push_back(v->ID);
		}
		v3 = v3->next_in_region;
	}
	if (left_next_v_in_region != nullptr) left_next_v_in_region->next_in_region = nullptr;
	if (right_next_v_in_region != nullptr) right_next_v_in_region->next_in_region = nullptr;

	if (left_root_v != nullptr) this->left_child = new Cluster(left_root_v, this->minimum_size, this->position_dimension, left_size);
	if (right_root_v != nullptr) this->right_child = new Cluster(right_root_v, this->minimum_size, this->position_dimension, right_size);

	return std::make_pair(this->left_child, this->right_child);
}