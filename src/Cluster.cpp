// #include "Cluster.h"
// #include <cmath>

// #include <iostream>
// #include <vector>
// #include <thread>
// #include <algorithm>
// #include <limits>
// #include <future>

// // #include "Vertex.h"

// Cluster::Cluster()
// {
// 	this->bit_sign = -1;
// 	this->dataset_size = 0;
// 	this->interval = 1.0;
// 	// this->vertice_id = nullptr;
// 	// this->vertice = nullptr;
// 	this->root_vertex = nullptr;
// 	this->minimum_size = 1;
// 	this->position_dimension = 2;
// 	this->left_child = nullptr;
// 	this->right_child = nullptr;
// }

// Cluster::Cluster(Vertex *root_vertex) : Cluster()
// {
// 	// this->vertice = vertice;
// 	// this->vertice_id = new int[this->dataset_size];//save the ids
// 	// for (int i = 0; i < this->dataset_size; i++)
// 	//{
// 	//	this->vertice_id[i] = ids[i];
// 	// }
// 	this->root_vertex = root_vertex;
// 	// this->interval = interval;
// 	// if (this->interval <= 0) throw std::runtime_error("The invalid interval,whcih cannot be smaller or equal to zero!");
// }

// Cluster::Cluster(Vertex *root_vertex, int minimum_size, int pisition_dimension, int dataset_size) : Cluster(root_vertex)
// {
// 	this->minimum_size = minimum_size;
// 	this->position_dimension = position_dimension;
// 	this->bit_sign = root_vertex->bit_sign;
// 	this->dataset_size = dataset_size;
// }

// Cluster::~Cluster()
// {
// 	// if (this->vertice_id != nullptr) delete[] this->vertice_id;
// 	if (this->left_child != nullptr)
// 		delete this->left_child;
// 	if (this->right_child != nullptr)
// 		delete this->right_child;
// }

// std::pair<Cluster *, Cluster *> Cluster::divide_cluster()
// {
// 	if (this->dataset_size <= 1 || this->dataset_size < this->minimum_size)
// 	{
// 		return {nullptr, nullptr}; // leaf, no children
// 	}
// 	std::vector<double> minimum_coordinates(position_dimension, std::numeric_limits<double>::max());
// 	std::vector<double> maximum_coordinates(position_dimension, std::numeric_limits<double>::min());
// 	std::vector<double> range_coordinates(position_dimension, 0.0);

// 	// search the minimun and maximun coordinates at each axis
// 	Vertex *v1 = this->root_vertex;
// 	int w = 0;
// 	while (v1 != nullptr)
// 	{
// 		for (int coord = 0; coord < position_dimension; coord++)
// 		{
// 			minimum_coordinates[coord] = (v1->positions[coord] < minimum_coordinates[coord]) ? v1->positions[coord] : minimum_coordinates[coord];
// 			maximum_coordinates[coord] = (v1->positions[coord] > maximum_coordinates[coord]) ? v1->positions[coord] : maximum_coordinates[coord];
// 		}
// 		v1 = v1->next_in_region;
// 		w++;
// 	}

// 	if (w != this->dataset_size)
// 	{
// 		throw std::runtime_error("Wrong dataset size!");
// 	}

// 	// calculate the range of axis and search the longest range
// 	double longest_range = std::numeric_limits<double>::min();
// 	int longest_coord = -1;
// 	for (int coord = 0; coord < position_dimension; coord++)
// 	{
// 		range_coordinates[coord] = maximum_coordinates[coord] - minimum_coordinates[coord];
// 		if (range_coordinates[coord] > longest_range)
// 		{
// 			longest_range = range_coordinates[coord];
// 			longest_coord = coord;
// 		}
// 	}

// 	// cannot conduct the operation
// 	if (longest_range <= 0)
// 		throw std::runtime_error("The invalid longest range!");
// 	// if (this->interval == 0) throw std::runtime_error("The invalid interval!");

// 	// statistics for the longest axis histogram
// 	this->interval = longest_range / 100.0;
// 	int num_interval = int(longest_range / this->interval + 1.0); // ����ȡ��
// 	std::vector<int> histogram(num_interval, 0);

// 	Vertex *v2 = this->root_vertex;
// 	while (v2 != nullptr)
// 	{
// 		double pos = v2->positions[longest_coord];
// 		int pos_his = int((pos - minimum_coordinates[longest_coord]) / this->interval); // ����ȡ��
// 		histogram[pos_his]++;
// 		v2 = v2->next_in_region;
// 	}

// 	// search for the median of the histogram
// 	int half_size = int(this->dataset_size * 0.5);
// 	int current_size = 0;
// 	int index = 0;
// 	while (current_size < half_size && index < num_interval)
// 	{
// 		current_size += histogram[index];
// 		index++;
// 	}
// 	double median = (index - 0.5) * this->interval + minimum_coordinates[longest_coord];

// 	// divide the dataset by the median
// 	Vertex *left_root_v = nullptr;
// 	Vertex *left_next_v_in_region = nullptr;
// 	int left_size = 0;
// 	Vertex *right_root_v = nullptr;
// 	Vertex *right_next_v_in_region = nullptr;
// 	int right_size = 0;

// 	Vertex *v3 = this->root_vertex;
// 	while (v3 != nullptr)
// 	{
// 		// Vertex *next = v3->next_in_region; // 1) remember where to go
// 		// v3->next_in_region = nullptr;	   // 2) sever from original list

// 		if (v3->positions[longest_coord] < median)
// 		{
// 			v3->append_sign(false);
// 			if (left_root_v == nullptr)
// 			{
// 				left_root_v = v3;
// 				left_next_v_in_region = v3;
// 				// left_next_v_in_region->next_in_region = nullptr;
// 			}
// 			else
// 			{
// 				left_next_v_in_region->next_in_region = v3;
// 				left_next_v_in_region = v3;
// 				// left_next_v_in_region->next_in_region = nullptr;
// 			}
// 			left_size++;
// 			// left_ids.push_back(v->ID);
// 		}
// 		else
// 		{
// 			v3->append_sign(true);
// 			if (right_root_v == nullptr)
// 			{
// 				right_root_v = v3;
// 				right_next_v_in_region = v3;
// 				// right_next_v_in_region->next_in_region = nullptr;
// 			}
// 			else
// 			{
// 				right_next_v_in_region->next_in_region = v3;
// 				right_next_v_in_region = v3;
// 				// right_next_v_in_region->next_in_region = nullptr;
// 			}
// 			right_size++;
// 			// right_ids.push_back(v->ID);
// 		}
// 		v3 = v3->next_in_region;
// 		// v3 = next;
// 	}

// 	std::cerr
// 		<< "[SPLIT] parent_size=" << this->dataset_size
// 		<< "  left=" << left_size
// 		<< "  right=" << right_size
// 		<< "  axis=" << longest_coord
// 		<< "  median=" << median
// 		<< "\n";
// 	if (left_size + right_size != this->dataset_size)
// 	{
// 		throw std::runtime_error("Split sizes don't add up!");
// 	}

// 	// after your histogram-based left_size/right_size calculation:

// 	if (left_size == 0 || right_size == 0)
// 	{
// 		// 1) Gather all coordinate values along the longest axis
// 		std::vector<double> coords;
// 		coords.reserve(dataset_size);
// 		for (Vertex *v = root_vertex; v; v = v->next_in_region)
// 			coords.push_back(v->positions[longest_coord]);

// 		// 2) Compute true median in O(n)
// 		std::nth_element(
// 			coords.begin(),
// 			coords.begin() + (dataset_size / 2),
// 			coords.end());
// 		double true_median = coords[dataset_size / 2];

// 		// 3) Reset split lists & counters
// 		left_root_v = left_next_v_in_region = nullptr;
// 		right_root_v = right_next_v_in_region = nullptr;
// 		left_size = right_size = 0;

// 		// 4) Re-split every node by true_median, severing from original
// 		Vertex *cursor = this->root_vertex;
// 		while (cursor)
// 		{
// 			Vertex *next = cursor->next_in_region;
// 			cursor->next_in_region = nullptr; // sever link

// 			if (cursor->positions[longest_coord] < true_median)
// 			{
// 				cursor->append_sign(false);
// 				if (!left_root_v)
// 				{
// 					left_root_v = left_next_v_in_region = cursor;
// 				}
// 				else
// 				{
// 					left_next_v_in_region->next_in_region = cursor;
// 					left_next_v_in_region = cursor;
// 				}
// 				++left_size;
// 			}
// 			else
// 			{
// 				cursor->append_sign(true);
// 				if (!right_root_v)
// 				{
// 					right_root_v = right_next_v_in_region = cursor;
// 				}
// 				else
// 				{
// 					right_next_v_in_region->next_in_region = cursor;
// 					right_next_v_in_region = cursor;
// 				}
// 				++right_size;
// 			}

// 			cursor = next;
// 		}

// 		// 5) Log & sanity-check
// 		std::cerr << "[FALLBACK] true_median=" << true_median
// 				  << "  left=" << left_size
// 				  << "  right=" << right_size << "\n";
// 		if (left_size + right_size != this->dataset_size)
// 		{
// 			throw std::runtime_error("Fallback split sizes mismatch");
// 		}
// 	}

// 	// if (left_next_v_in_region != nullptr)
// 	// 	left_next_v_in_region->next_in_region = nullptr;
// 	// if (right_next_v_in_region != nullptr)
// 	// 	right_next_v_in_region->next_in_region = nullptr;

// 	// if (left_root_v != nullptr)
// 	// 	this->left_child = new Cluster(left_root_v, this->minimum_size, this->position_dimension, left_size);
// 	// if (right_root_v != nullptr)
// 	// 	this->right_child = new Cluster(right_root_v, this->minimum_size, this->position_dimension, right_size);

// 	// return std::make_pair(this->left_child, this->right_child);

// 	Cluster *L = nullptr;
// 	Cluster *R = nullptr;
// 	if (left_root_v)
// 	{
// 		L = new Cluster(left_root_v,
// 						this->minimum_size,
// 						this->position_dimension,
// 						left_size);
// 	}
// 	if (right_root_v)
// 	{
// 		R = new Cluster(right_root_v,
// 						this->minimum_size,
// 						this->position_dimension,
// 						right_size);
// 	}
// 	return std::make_pair(L, R);
// }

// MinMax parallel_find_min_max(Vertex *start_vertex, int num_vertices, int position_dimension)
// {
// 	MinMax result;
// 	result.minimum_coordinates.resize(position_dimension, std::numeric_limits<double>::max());
// 	result.maximum_coordinates.resize(position_dimension, std::numeric_limits<double>::min());

// 	Vertex *v = start_vertex;
// 	for (int i = 0; i < num_vertices && v != nullptr; ++i)
// 	{
// 		for (int coord = 0; coord < position_dimension; coord++)
// 		{
// 			result.minimum_coordinates[coord] = std::min(v->positions[coord], result.minimum_coordinates[coord]);
// 			result.maximum_coordinates[coord] = std::max(v->positions[coord], result.maximum_coordinates[coord]);
// 		}
// 		v = v->next_in_region;
// 	}

// 	return result;
// }

// void parallel_histogram(Vertex *start_vertex, int num_vertices, int longest_coord, double minimum_coordinate, double interval, std::vector<int> &histogram)
// {
// 	Vertex *v = start_vertex;
// 	for (int i = 0; i < num_vertices && v != nullptr; ++i)
// 	{
// 		double pos = v->positions[longest_coord];
// 		int pos_his = int((pos - minimum_coordinate) / interval);
// 		histogram[pos_his]++;
// 		v = v->next_in_region;
// 	}
// }

// void parallel_divide_vertices(Vertex *start_vertex, int num_vertices, int longest_coord, double median, Vertex **left_root, Vertex **left_tail, int *left_size, Vertex **right_root, Vertex **right_tail, int *right_size)
// {
// 	Vertex *v = start_vertex;
// 	for (int i = 0; i < num_vertices && v != nullptr; ++i)
// 	{
// 		if (v->positions[longest_coord] < median)
// 		{
// 			v->append_sign(false);
// 			if (*left_root == nullptr)
// 			{
// 				*left_root = v;
// 				*left_tail = v;
// 			}
// 			else
// 			{
// 				(*left_tail)->next_in_region = v;
// 				*left_tail = v;
// 			}
// 			(*left_size)++;
// 		}
// 		else
// 		{
// 			v->append_sign(true);
// 			if (*right_root == nullptr)
// 			{
// 				*right_root = v;
// 				*right_tail = v;
// 			}
// 			else
// 			{
// 				(*right_tail)->next_in_region = v;
// 				*right_tail = v;
// 			}
// 			(*right_size)++;
// 		}
// 		v = v->next_in_region;
// 	}
// }

// std::pair<Cluster *, Cluster *> Cluster::divide_cluster_parallel()
// {
// 	int num_threads = std::thread::hardware_concurrency();

// 	// ���в�����С���������
// 	std::vector<std::future<MinMax>> min_max_futures;
// 	int chunk_size = this->dataset_size / num_threads;
// 	Vertex *current_vertex = this->root_vertex;
// 	for (int i = 0; i < num_threads; ++i)
// 	{
// 		int start_index = i * chunk_size;
// 		int end_index = (i == num_threads - 1) ? this->dataset_size : (i + 1) * chunk_size;
// 		int num_vertices = end_index - start_index;
// 		min_max_futures.push_back(std::async(std::launch::async, parallel_find_min_max, current_vertex, num_vertices, this->position_dimension));
// 		for (int j = 0; j < num_vertices && current_vertex != nullptr; ++j)
// 		{
// 			current_vertex = current_vertex->next_in_region;
// 		}
// 	}

// 	std::vector<double> minimum_coordinates(position_dimension, std::numeric_limits<double>::max());
// 	std::vector<double> maximum_coordinates(position_dimension, std::numeric_limits<double>::min());
// 	for (auto &future : min_max_futures)
// 	{
// 		MinMax result = future.get();
// 		for (int coord = 0; coord < position_dimension; coord++)
// 		{
// 			minimum_coordinates[coord] = std::min(result.minimum_coordinates[coord], minimum_coordinates[coord]);
// 			maximum_coordinates[coord] = std::max(result.maximum_coordinates[coord], maximum_coordinates[coord]);
// 		}
// 	}

// 	// ���㷶Χ���ҵ����Χ
// 	double longest_range = std::numeric_limits<double>::min();
// 	int longest_coord = -1;
// 	for (int coord = 0; coord < position_dimension; coord++)
// 	{
// 		double range = maximum_coordinates[coord] - minimum_coordinates[coord];
// 		if (range > longest_range)
// 		{
// 			longest_range = range;
// 			longest_coord = coord;
// 		}
// 	}

// 	// ������
// 	if (longest_range <= 0)
// 		throw std::runtime_error("The invalid longest range!");
// 	// if (this->interval == 0) throw std::runtime_error("The invalid interval!");

// 	// statistics for the longest axis histogram
// 	int num_interval = int(longest_range / this->interval + 1.0); // ����ȡ��
// 	std::vector<int> histogram(num_interval, 0);

// 	Vertex *v2 = this->root_vertex;
// 	while (v2 != nullptr)
// 	{
// 		double pos = v2->positions[longest_coord];
// 		int pos_his = int((pos - minimum_coordinates[longest_coord]) / this->interval); // ����ȡ��
// 		histogram[pos_his]++;
// 		v2 = v2->next_in_region;
// 	}

// 	// search for the median of the histogram
// 	int half_size = int(this->dataset_size * 0.5);
// 	int current_size = 0;
// 	int index = 0;
// 	while (current_size < half_size && index < num_interval)
// 	{
// 		current_size += histogram[index];
// 		index++;
// 	}
// 	double median = (index - 0.5) * this->interval + minimum_coordinates[longest_coord];

// 	// divide the dataset by the median
// 	Vertex *left_root_v = nullptr;
// 	Vertex *left_next_v_in_region = nullptr;
// 	int left_size = 0;
// 	Vertex *right_root_v = nullptr;
// 	Vertex *right_next_v_in_region = nullptr;
// 	int right_size = 0;

// 	Vertex *v3 = this->root_vertex;
// 	while (v3 != nullptr)
// 	{
// 		if (v3->positions[longest_coord] < median)
// 		{
// 			v3->append_sign(false);
// 			if (left_root_v == nullptr)
// 			{
// 				left_root_v = v3;
// 				left_next_v_in_region = v3;
// 				// left_next_v_in_region->next_in_region = nullptr;
// 			}
// 			else
// 			{
// 				left_next_v_in_region->next_in_region = v3;
// 				left_next_v_in_region = v3;
// 				// left_next_v_in_region->next_in_region = nullptr;
// 			}
// 			left_size++;
// 			// left_ids.push_back(v->ID);
// 		}
// 		else
// 		{
// 			v3->append_sign(true);
// 			if (right_root_v == nullptr)
// 			{
// 				right_root_v = v3;
// 				right_next_v_in_region = v3;
// 				// right_next_v_in_region->next_in_region = nullptr;
// 			}
// 			else
// 			{
// 				right_next_v_in_region->next_in_region = v3;
// 				right_next_v_in_region = v3;
// 				// right_next_v_in_region->next_in_region = nullptr;
// 			}
// 			right_size++;
// 			// right_ids.push_back(v->ID);
// 		}
// 		v3 = v3->next_in_region;
// 	}
// 	if (left_next_v_in_region != nullptr)
// 		left_next_v_in_region->next_in_region = nullptr;
// 	if (right_next_v_in_region != nullptr)
// 		right_next_v_in_region->next_in_region = nullptr;

// 	if (left_root_v != nullptr)
// 		this->left_child = new Cluster(left_root_v, this->minimum_size, this->position_dimension, left_size);
// 	if (right_root_v != nullptr)
// 		this->right_child = new Cluster(right_root_v, this->minimum_size, this->position_dimension, right_size);

// 	return std::make_pair(this->left_child, this->right_child);
// }

#include "Cluster.h"
#include <cmath>

#include <iostream>
#include <vector>
#include <thread>
#include <algorithm>
#include <limits>
#include <future>

// #include "Vertex.h"

Cluster::Cluster()
{
	this->bit_sign = -1;
	this->dataset_size = 0;
	this->interval = 1.0;
	// this->vertice_id = nullptr;
	// this->vertice = nullptr;
	this->root_vertex = nullptr;
	this->minimum_size = 1;
	this->position_dimension = -1; // 这里需要与输入数据的维度保持一致,为了防止调用默认构造函数时出现维度未设置的为题，这里默认设置为-1
	this->left_child = nullptr;
	this->right_child = nullptr;
}

Cluster::Cluster(Vertex *root_vertex) : Cluster()
{
	// this->vertice = vertice;
	// this->vertice_id = new int[this->dataset_size];//save the ids
	// for (int i = 0; i < this->dataset_size; i++)
	//{
	//	this->vertice_id[i] = ids[i];
	// }
	this->root_vertex = root_vertex;
	this->position_dimension = root_vertex->position_dimension;
	// this->interval = interval;
	// if (this->interval <= 0) throw std::runtime_error("The invalid interval,whcih cannot be smaller or equal to zero!");
}

Cluster::Cluster(Vertex *root_vertex, int minimum_size, int pisition_dimension, int dataset_size) : Cluster(root_vertex)
{
	this->minimum_size = minimum_size;
	this->position_dimension = position_dimension;
	this->bit_sign = root_vertex->bit_sign;
	this->dataset_size = dataset_size;
}

Cluster::~Cluster()
{
	// if (this->vertice_id != nullptr) delete[] this->vertice_id;
	// if (this->vertice_id != nullptr) delete[] this->vertice_id;
	if (this->left_child != nullptr)
		delete this->left_child;
	if (this->right_child != nullptr)
		delete this->right_child;
}

// std::pair<Cluster*, Cluster*> Cluster::divide_cluster()
//{
//	if (this->dataset_size <= 1 || this->dataset_size < this->minimum_size)
//	{
//		return { nullptr, nullptr }; // leaf, no children
//	}
//	std::vector<double> minimum_coordinates(position_dimension, std::numeric_limits<double>::max());
//	std::vector<double> maximum_coordinates(position_dimension, std::numeric_limits<double>::min());
//	std::vector<double> range_coordinates(position_dimension, 0.0);
//
//	// search the minimun and maximun coordinates at each axis
//	Vertex* v1 = this->root_vertex;
//	int w = 0;
//	while (v1 != nullptr)
//	{
//		for (int coord = 0; coord < position_dimension; coord++)
//		{
//			minimum_coordinates[coord] = (v1->positions[coord] < minimum_coordinates[coord]) ? v1->positions[coord] : minimum_coordinates[coord];
//			maximum_coordinates[coord] = (v1->positions[coord] > maximum_coordinates[coord]) ? v1->positions[coord] : maximum_coordinates[coord];
//		}
//		v1 = v1->next_in_region;
//		w++;
//	}
//
//	std::cout << w << " " << this->dataset_size << std::endl;
//	if (w != this->dataset_size)
//	{
//		throw std::runtime_error("Wrong dataset size!");
//	}
//
//	// calculate the range of axis and search the longest range
//	double longest_range = std::numeric_limits<double>::min();
//	int longest_coord = -1;
//	for (int coord = 0; coord < position_dimension; coord++)
//	{
//		range_coordinates[coord] = maximum_coordinates[coord] - minimum_coordinates[coord];
//		if (range_coordinates[coord] > longest_range)
//		{
//			longest_range = range_coordinates[coord];
//			longest_coord = coord;
//		}
//	}
//
//	// cannot conduct the operation
//	if (longest_range <= 0)
//		throw std::runtime_error("The invalid longest range!");
//	// if (this->interval == 0) throw std::runtime_error("The invalid interval!");
//
//	// statistics for the longest axis histogram
//	this->interval = longest_range / 100.0;
//	int num_interval = int(longest_range / this->interval + 1.0); //
//	std::vector<int> histogram(num_interval, 0);
//
//	Vertex* v2 = this->root_vertex;
//	while (v2 != nullptr)
//	{
//		double pos = v2->positions[longest_coord];
//		int pos_his = int((pos - minimum_coordinates[longest_coord]) / this->interval); //
//		histogram[pos_his]++;
//		v2 = v2->next_in_region;
//	}
//
//	// search for the median of the histogram
//	int half_size = int(this->dataset_size * 0.5);
//	int current_size = 0;
//	int index = 0;
//	while (current_size < half_size && index < num_interval)
//	{
//		current_size += histogram[index];
//		index++;
//	}
//	double median = (index - 0.5) * this->interval + minimum_coordinates[longest_coord];
//
//	// divide the dataset by the median
//	Vertex* left_root_v = nullptr;
//	Vertex* left_next_v_in_region = nullptr;
//	int left_size = 0;
//	Vertex* right_root_v = nullptr;
//	Vertex* right_next_v_in_region = nullptr;
//	int right_size = 0;
//
//	Vertex* v3 = this->root_vertex;
//	while (v3 != nullptr)
//	{
//		// Vertex *next = v3->next_in_region; // 1) remember where to go
//		// v3->next_in_region = nullptr;	   // 2) sever from original list
//
//		if (v3->positions[longest_coord] < median)
//		{
//			v3->append_sign(false);
//			if (left_root_v == nullptr)
//			{
//				left_root_v = v3;
//				left_next_v_in_region = v3;
//				// left_next_v_in_region->next_in_region = nullptr;
//			}
//			else
//			{
//				left_next_v_in_region->next_in_region = v3;
//				left_next_v_in_region = v3;
//				// left_next_v_in_region->next_in_region = nullptr;
//			}
//			left_size++;
//			// left_ids.push_back(v->ID);
//		}
//		else
//		{
//			v3->append_sign(true);
//			if (right_root_v == nullptr)
//			{
//				right_root_v = v3;
//				right_next_v_in_region = v3;
//				// right_next_v_in_region->next_in_region = nullptr;
//			}
//			else
//			{
//				right_next_v_in_region->next_in_region = v3;
//				right_next_v_in_region = v3;
//				// right_next_v_in_region->next_in_region = nullptr;
//			}
//			right_size++;
//			// right_ids.push_back(v->ID);
//		}
//		v3 = v3->next_in_region;
//		// v3 = next;
//	}
//
//	std::cerr
//		<< "[SPLIT] parent_size=" << this->dataset_size
//		<< "  left=" << left_size
//		<< "  right=" << right_size
//		<< "  axis=" << longest_coord
//		<< "  median=" << median
//		<< "\n";
//	if (left_size + right_size != this->dataset_size)
//	{
//		throw std::runtime_error("Split sizes don't add up!");
//	}
//
//	// after your histogram-based left_size/right_size calculation:
//
//	if (left_size == 0 || right_size == 0)
//	{
//		// 1) Gather all coordinate values along the longest axis
//		std::vector<double> coords;
//		coords.reserve(dataset_size);
//		for (Vertex* v = root_vertex; v; v = v->next_in_region)
//			coords.push_back(v->positions[longest_coord]);
//
//		// 2) Compute true median in O(n)
//		std::nth_element(
//			coords.begin(),
//			coords.begin() + (dataset_size / 2),
//			coords.end());
//		double true_median = coords[dataset_size / 2];
//
//		// 3) Reset split lists & counters
//		left_root_v = left_next_v_in_region = nullptr;
//		right_root_v = right_next_v_in_region = nullptr;
//		left_size = right_size = 0;
//
//		// 4) Re-split every node by true_median, severing from original
//		Vertex* cursor = this->root_vertex;
//		while (cursor)
//		{
//			Vertex* next = cursor->next_in_region;
//			cursor->next_in_region = nullptr; // sever link
//
//			if (cursor->positions[longest_coord] < true_median)
//			{
//				cursor->append_sign(false);
//				if (!left_root_v)
//				{
//					left_root_v = left_next_v_in_region = cursor;
//				}
//				else
//				{
//					left_next_v_in_region->next_in_region = cursor;
//					left_next_v_in_region = cursor;
//				}
//				++left_size;
//			}
//			else
//			{
//				cursor->append_sign(true);
//				if (!right_root_v)
//				{
//					right_root_v = right_next_v_in_region = cursor;
//				}
//				else
//				{
//					right_next_v_in_region->next_in_region = cursor;
//					right_next_v_in_region = cursor;
//				}
//				++right_size;
//			}
//
//			cursor = next;
//		}
//
//		// 5) Log & sanity-check
//		std::cerr << "[FALLBACK] true_median=" << true_median
//			<< "  left=" << left_size
//			<< "  right=" << right_size << "\n";
//		if (left_size + right_size != this->dataset_size)
//		{
//			throw std::runtime_error("Fallback split sizes mismatch");
//		}
//	}
//
//	// if (left_next_v_in_region != nullptr)
//	// 	left_next_v_in_region->next_in_region = nullptr;
//	// if (right_next_v_in_region != nullptr)
//	// 	right_next_v_in_region->next_in_region = nullptr;
//
//	// if (left_root_v != nullptr)
//	// 	this->left_child = new Cluster(left_root_v, this->minimum_size, this->position_dimension, left_size);
//	// if (right_root_v != nullptr)
//	// 	this->right_child = new Cluster(right_root_v, this->minimum_size, this->position_dimension, right_size);
//
//	// return std::make_pair(this->left_child, this->right_child);
//
//	Cluster* L = nullptr;
//	Cluster* R = nullptr;
//	if (left_root_v)
//	{
//		L = new Cluster(left_root_v,
//			this->minimum_size,
//			this->position_dimension,
//			left_size);
//		std::cout <<"L size: " << left_size << " DS size: " << L->dataset_size << std::endl;
//	}
//	if (right_root_v)
//	{
//		R = new Cluster(right_root_v,
//			this->minimum_size,
//			this->position_dimension,
//			right_size);
//		std::cout << "R size: " << right_size << " DS size: " << R->dataset_size << std::endl;
//
//	}
//	return std::make_pair(L, R);
// }

std::pair<Cluster *, Cluster *> Cluster::divide_cluster()
{
	// 方法最大维度的数据排序，然后找到中间的中位数
	std::pair<int, int> DIM_count = this->select_split_DIM(this->root_vertex);

	// 对方差最大维度的数据进行排序
	std::vector<double> nums(DIM_count.second, 0.0);
	Vertex *next = this->root_vertex;
	int id = 0;
	while (next != nullptr)
	{
		nums[id] = next->positions[DIM_count.first];
		next = next->next_in_region;
		++id;
	}
	double median = this->find_median(nums);

	// divide the dataset by the median
	Vertex *left_root_v = nullptr;
	Vertex *left_next_v_in_region = nullptr;
	int left_size = 0;
	Vertex *right_root_v = nullptr;
	Vertex *right_next_v_in_region = nullptr;
	int right_size = 0;

	Vertex *v3 = this->root_vertex;
	while (v3 != nullptr)
	{
		// 这里需要<=median的数据，且应该数量达到总数一半，因为"<=median"时，可能存在大量相同的值，导致数量远远大于一半的情况。
		if (v3->positions[DIM_count.first] <= median && left_size <= int(DIM_count.second * 0.5))
		{
			v3->append_sign(false);
			if (left_root_v == nullptr)
			{
				left_root_v = v3;
				left_next_v_in_region = v3;
			}
			else
			{
				left_next_v_in_region->next_in_region = v3;
				left_next_v_in_region = v3;
			}
			left_size++;
		}
		else
		{
			v3->append_sign(true);
			if (right_root_v == nullptr)
			{
				right_root_v = v3;
				right_next_v_in_region = v3;
			}
			else
			{
				right_next_v_in_region->next_in_region = v3;
				right_next_v_in_region = v3;
			}
			right_size++;
		}
		v3 = v3->next_in_region;
	}
	if (left_next_v_in_region != nullptr)
		left_next_v_in_region->next_in_region = nullptr;
	if (right_next_v_in_region != nullptr)
		right_next_v_in_region->next_in_region = nullptr;

	if (left_root_v != nullptr)
		this->left_child = new Cluster(left_root_v, this->minimum_size, this->position_dimension, left_size);
	if (right_root_v != nullptr)
		this->right_child = new Cluster(right_root_v, this->minimum_size, this->position_dimension, right_size);

	return std::make_pair(this->left_child, this->right_child);
}
/// <summary>
/// 计算root即后续数据中方差最大的维度
/// </summary>
/// <param name="root">root</param>
/// <returns>方差最大的维度,元素数量</returns>
std::pair<int, int> Cluster::select_split_DIM(Vertex *root)
{
	// 计算每个维度的平均值
	int DIM = root->position_dimension;
	int count = 0;
	std::vector<double> means(DIM, 0.0);
	Vertex *next = root;
	while (next != nullptr)
	{
		for (int d = 0; d < DIM; ++d)
		{
			means[d] = (next->positions[d] + means[d] * count) / (count + 1);
		}
		++count;
		next = next->next_in_region;
	}

	// 计算每个维度的方差
	std::vector<double> variances(DIM, 0.0); // 初始化大小为 DIM
	next = root;
	while (next != nullptr)
	{
		for (size_t d = 0; d < DIM; ++d)
		{
			float diff = next->positions[d] - means[d];
			variances[d] += std::abs(diff); // 方差太大了，这里直接用的绝对距离
		}
		next = next->next_in_region;
	}

	// 返回 {最大方差的维度, count}
	auto max_it = std::max_element(variances.begin(), variances.end());
	int max_dim = std::distance(variances.begin(), max_it);
	// std::cout << max_dim << " " << count << std::endl;
	return {max_dim, count};
}

double Cluster::find_median(std::vector<double> &nums)
{
	size_t n = nums.size();
	if (n == 0)
		return 0.0;
	auto mid = nums.begin() + n / 2;
	std::nth_element(nums.begin(), mid, nums.end());

	if (n % 2 != 0)
	{
		return *mid;
	}
	else
	{
		// 对于偶数个元素，需要找到中间两个元素的平均值
		auto left_mid = std::max_element(nums.begin(), mid);
		return (*left_mid + *mid) / 2.0;
	}
}