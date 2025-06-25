#pragma once
enum class State;
class Vertex;

class Edge
{
public:
	Edge();
	Edge(int, Vertex *, Vertex *);
	Edge(int, double, Vertex *, Vertex *);
	~Edge();

public:
	int ID;
	int shortest_arc_num;
	// int type;
	int length;
	double weight;
	State state;
	Vertex *left_vertex;
	Vertex *right_vertex;
	Edge *next_edge_in_region;
	bool is_checked;
	// Edge* direct;

public:
	Vertex *get_other_vertex(Vertex *);
	virtual double GetWeight() = 0;

	bool compareByValue(const Edge &other) const;
	bool compareByAddress(const Edge &other) const;
	bool compareByEnds(const Edge &other) const;
	bool operator<(const Edge &other) const;
	bool operator>(const Edge &other) const;
	bool operator==(const Edge &other) const;
	bool operator!=(const Edge &other) const;
};
