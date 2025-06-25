#pragma once
#include <iostream>
#include <vector>

#include "Edge.h"
#include "NNG.h"

class RootFeatures;
class Edge;

enum class State
{
	Activate,
	Suspend,
	Freeze
};

// template <typename T>
class Vertex
{
public:
	Vertex();
	Vertex(int id, int connection_degree, double radius);
	Vertex(int ID, RootFeatures *features);
	Vertex(const Vertex &other);
	~Vertex();

public:
	void push(Edge *ptr_e);
	void push_range(int *first, int *last);
	void add_next_vertex(Vertex *ptr_next);
	void remove_edge(int edge_id);
	void change_direct(Edge *new_direct, NNG *nng);
	void set_position_dimention(int dimention);
	void append_sign(bool sign);
	int get_sign();
	void reset();
	bool operator<(const Vertex &other) const;
	Vertex &operator=(const Vertex &other);
	virtual int GetID() = 0;

public:
	int ID;
	int position_dimension;
	double *positions;
	double radius;
	int connection_degree;
	int color;
	int size;
	State state;
	Vertex *parent;
	Vertex *next;
	Vertex *next_in_region;
	Vertex *last;
	bool is_visited;
	bool is_built_tree = false;
	std::vector<int> edges;
	std::vector<int> merged_vertice;
	Edge *direct;
	RootFeatures *ptr_features;

	int bit_sign;
	int bits_num;
};