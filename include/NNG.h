#pragma once
#include<iostream>
//#include <set>
#include<map>
#include"Edge.h"
//#include"Vertex.h"
class Vertex;
class ValueCompare {
	bool operator()(const std::pair<Edge*, double>& a, const std::pair<Edge*, double>& b) const {
		return a.second < b.second;
	}
};
class NNG
{
public:
	NNG();
	~NNG();

public:
	Edge* pop();
	void push(Edge*);
	bool find(Edge*);
	void remove(Edge*);
	int size();
	bool empty();
	void show_map();
	void clear();
private:
	std::multimap<double, Edge*> map;
};

