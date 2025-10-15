// #pragma once
// #include <iostream>
// #include "Vertex.h"
// using namespace std;
// struct MinMax {
// 	std::vector<double> minimum_coordinates;
// 	std::vector<double> maximum_coordinates;
// };
// class Cluster
// {
// public:
// 	Cluster();
// 	Cluster(Vertex* root_vertex);
// 	Cluster(Vertex* root_vertex, int minimum_size, int pisition_dimension, int dataset_size);

// 	~Cluster();
// 	std::pair<Cluster*, Cluster*> divide_cluster();
// 	std::pair<Cluster*, Cluster*> divide_cluster_parallel();

// public:
// 	//int* vertice_id;
// 	//Vertex* vertice;
// 	Vertex* root_vertex;
// 	int dataset_size;//equals to vertex num
// 	int bit_sign;

// private:
// 	double interval;
// 	int position_dimension;
// 	int minimum_size;
// 	Cluster* left_child;
// 	Cluster* right_child;
// 	//MinMax parallel_find_min_max(Vertex* start_vertex, int num_vertices, int position_dimension);
// 	//void parallel_histogram(Vertex* start_vertex, int num_vertices, int longest_coord, double minimum_coordinate, double interval, std::vector<int>& histogram);
// 	//void parallel_divide_vertices(Vertex* start_vertex, int num_vertices, int longest_coord, double median, Vertex** left_root, Vertex** left_tail, int* left_size, Vertex** right_root, Vertex** right_tail, int* right_size);

// };

#pragma once
#include <iostream>
#include "Vertex.h"
using namespace std;
struct MinMax
{
	std::vector<double> minimum_coordinates;
	std::vector<double> maximum_coordinates;
};
class Cluster
{
public:
	Cluster();
	Cluster(Vertex *root_vertex);
	Cluster(Vertex *root_vertex, int minimum_size, int pisition_dimension, int dataset_size);

	~Cluster();
	std::pair<Cluster *, Cluster *> divide_cluster();
	// std::pair<Cluster*, Cluster*> divide_cluster_parallel();

public:
	// int* vertice_id;
	// Vertex* vertice;
	Vertex *root_vertex;
	int dataset_size; // equals to vertex num
	int bit_sign;

private:
	double interval;
	int position_dimension;
	int minimum_size;
	Cluster *left_child;
	Cluster *right_child;
	std::pair<int, int> select_split_DIM(Vertex *root);
	double find_median(std::vector<double> &nums);

	// MinMax parallel_find_min_max(Vertex* start_vertex, int num_vertices, int position_dimension);
	// void parallel_histogram(Vertex* start_vertex, int num_vertices, int longest_coord, double minimum_coordinate, double interval, std::vector<int>& histogram);
	// void parallel_divide_vertices(Vertex* start_vertex, int num_vertices, int longest_coord, double median, Vertex** left_root, Vertex** left_tail, int* left_size, Vertex** right_root, Vertex** right_tail, int* right_size);
};
