#pragma once
#include <type_traits>
#include <chrono>
#include <queue>
#include <limits.h>
#include <future>
#include <shared_mutex>
#include "Vertex.h"
#include "Edge.h"
#include "NNG.h"
#include "Cluster.h"

#define RED_TEXT "\033[31m"
#define GREEN_TEXT "\033[32m"
#define RESET_COLOR "\033[0m"

enum class RunMode
{
	None,
	Region,
	Pruning,
	Cycle,
	Region_Pruning,
	Region_Cycle,
	Pruning_Cycle,
	Region_Pruning_Cycle,
	Region_Pruning_Cycle_MultiCores,
	Region_Pruning_Cycle_GPU
};

// �Զ���ȽϺ���
struct CompareAPointers
{
	bool operator()(const Vertex *a, const Vertex *b) const
	{
		return a->edges.size() < b->edges.size(); // ע�⣺������С�ںţ���Ϊ������Ҫ����
	}
};

template <typename T1, typename T2, typename T_Feature>
class RegionAdjacencyGraph
{
	static_assert(std::is_base_of_v<Vertex, T1>, "T1 must be derived from Vertex Class");
	static_assert(std::is_base_of_v<Edge, T2>, "T2 must be derived from Edge Class");
	static_assert(std::is_base_of_v<RootFeatures, T_Feature>, "T_Feature must be derived from RootFeatures Class");

public:
	typedef void (*UpdateFeatures)(T_Feature *, T_Feature *, T2 *);
	typedef double (*MergingCriteria)(T_Feature *, T_Feature *, T2 *);

public:
	RegionAdjacencyGraph(int vertice_num, int edge_num, UpdateFeatures update_features, MergingCriteria MC) : max_priqueue(), update_features_func(update_features), merging_criteria(MC)
	{
		this->vertice = nullptr;
		this->edges = nullptr;
		this->nng = new NNG();
		this->vertice_num = vertice_num;
		this->edge_num = edge_num;
	};

	RegionAdjacencyGraph() : max_priqueue() {

							 };

	~RegionAdjacencyGraph()
	{
		delete this->nng;
	};

	T1 *run(double threshold, RunMode mode, int region_size = 2500)
	{
		T1 *root = nullptr;
		switch (mode)
		{
		case RunMode::None:
			root = this->global_best_optimization(threshold);
			break;
		case RunMode::Region:
			std::cout << "[INFO] Calling func\"region_optimization\"" << std::endl;
			root = this->region_optimization(threshold, region_size);
			break;
		case RunMode::Pruning:
			root = this->purning_optimization(threshold);
			break;
		case RunMode::Cycle:
			root = this->cycle_optimization(threshold);
			break;
		case RunMode::Region_Pruning:
			root = this->region_pruning_optimization(threshold, region_size);
			break;
		case RunMode::Region_Cycle:
			root = this->region_cycle_optimization(threshold, region_size);
			break;
		case RunMode::Pruning_Cycle:
			root = this->pruning_cycle_optimization(threshold);
			break;
		case RunMode::Region_Pruning_Cycle:
			root = this->region_pruning_cycle_optimization(threshold, region_size);
			break;
		case RunMode::Region_Pruning_Cycle_MultiCores:
			root = this->region_pruning_cycle_optimization_MultiCores(threshold, region_size);
			break;
		case RunMode::Region_Pruning_Cycle_GPU:
			break;
		default:
			break;
		}
		return root;
	}

private:
	void construct_nng()
	{
		if (this->vertice == nullptr)
			throw std::runtime_error("The function construct_nng() must be called after the function connect_vertex_and_edge()");
		auto start_time = std::chrono::high_resolution_clock::now();
		if (this->nng != nullptr)
			this->nng->clear();
		for (int i = 0; i < this->edge_num; i++)
		{
			T2 *e = &this->edges[i];
			e->shortest_arc_num = 0;
		}
		for (int i = 0; i < this->vertice_num; i++)
		{
			T1 *v = &this->vertice[i];
			if (v->state != State::Freeze)
			{
				this->initial_direct(v, this->edges);
			}
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
		// std::cout << "vertices����������: " << RED_TEXT << duration.count() << RESET_COLOR << "��" << std::endl;
		if (duration.count() < 10)
			std::cout << "The directions of vertices are calculated: " << GREEN_TEXT << duration.count() << RESET_COLOR << "s" << std::endl;
		else
			std::cout << "The directions of vertices are calculated: " << RED_TEXT << duration.count() << RESET_COLOR << "s" << std::endl;

		if (this->nng != nullptr)
		{
			for (int i = 0; i < this->edge_num; i++)
			{
				T2 *e = &this->edges[i];
				if (e->state == State::Activate && e->shortest_arc_num == 2)
				{
					this->nng->push(&this->edges[i]);
				}
			}
		}
		else
		{
			throw std::runtime_error("nng cannot be null!");
		}
		start_time = std::chrono::high_resolution_clock::now();
		end_time = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
		// std::cout << "NNG��ʼ�����ȶ������: " << RED_TEXT << duration.count() << RESET_COLOR << "��" << std::endl;
		// std::cout << "nng size: " << GREEN_TEXT << this->nng->size() << RESET_COLOR <<std::endl;
		if (duration.count() < 10)
			std::cout << "The priority queue in NNG is initialized: " << GREEN_TEXT << duration.count() << RESET_COLOR << "s" << std::endl;
		else
			std::cout << "The priority queue in NNG is initialized: " << RED_TEXT << duration.count() << RESET_COLOR << "s" << std::endl;
	}

	void construct_nng(NNG *region_nng)
	{
		if (this->vertice == nullptr)
			throw std::runtime_error("The function construct_nng() must be called after the function connect_vertex_and_edge()");
		auto start_time = std::chrono::high_resolution_clock::now();

		if (region_nng != nullptr)
			region_nng->clear();

		for (int i = 0; i < this->edge_num; i++)
		{
			T2 *e = &this->edges[i];
			if (e->state != State::Freeze)
			{
				e->shortest_arc_num = 0;
			}
		}

		for (int i = 0; i < this->vertice_num; i++)
		{
			T1 *v = &this->vertice[i];
			v->direct = nullptr;
			if (v->state != State::Freeze)
			{
				this->initial_direct(v, this->edges);
			}
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
		// std::cout << "vertices����������: " << RED_TEXT << duration.count() << RESET_COLOR << "��" << std::endl;
		// if (duration.count() < 10) std::cout << "The directions of vertices are calculated: " << GREEN_TEXT << duration.count() << RESET_COLOR << "s" << std::endl;
		// else 					   std::cout << "The directions of vertices are calculated: " << RED_TEXT << duration.count() << RESET_COLOR << "s" << std::endl;

		if (region_nng != nullptr)
		{
			for (int i = 0; i < this->edge_num; i++)
			{
				T2 *e = &this->edges[i];
				if (e->state == State::Activate && e->shortest_arc_num == 2)
				{
					region_nng->push(&this->edges[i]);
				}
			}
		}
		else
		{
			throw std::runtime_error("nng cannot be null!");
		}

		start_time = std::chrono::high_resolution_clock::now();
		end_time = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	}

	void initial_direct(T1 *v, T2 *edge_array)
	{
		double min_w = std::numeric_limits<double>::max();
		for (int i = 0; i < v->edges.size(); i++)
		{
			int id = v->edges[i];
			T2 *ptr_e = &edge_array[id];
			if (ptr_e->weight < min_w)
			{
				v->direct = ptr_e;
				min_w = ptr_e->weight;
			}
		}

		if (v->direct != nullptr)
		{
			v->direct->shortest_arc_num++;
		}
	}

	std::pair<T2 *, T1 *> is_edge_connected_vertex(T2 *test_e, T1 *aim_v)
	{
		for (size_t i = 0; i < aim_v->edges.size(); ++i)
		{
			int id = aim_v->edges[i];
			T2 *e_has_in_aim_v = &this->edges[id];
			if (e_has_in_aim_v->compareByEnds(*test_e))
			{
				return std::make_pair(e_has_in_aim_v, (T1 *)e_has_in_aim_v->get_other_vertex(aim_v));
			}
		}
		return std::make_pair(nullptr, nullptr);
	}

	void update_direction_cycle_in_2nd_degree_vertices(T1 *v, NNG *nng_ptr)
	{

		// Edge* min_edge = v->get_min_w_edge(this->edges);
		T2 *min_edge = this->get_min_w_edge((T1 *)v, this->edges);
		if (v->direct != min_edge)
		{
			// if (v->direct->shortest_arc_num == 2) this->nng->remove(v->direct);
			// v->change_direct(min_edge, this->nng);
			v->change_direct(min_edge, nng_ptr);
			// if (v->direct->shortest_arc_num == 2) this->nng->push(v->direct);
		}
		// ���û�䷽��ԭʼ��RAG�е�����cycle����nng�У���������do nothing
	}

	T2 *get_min_w_edge(T1 *v, T2 *edge_array)
	{
		double min_w = std::numeric_limits<double>::max();
		T2 *ptr_e = nullptr;
		for (int i = 0; i < v->edges.size(); i++)
		{
			int id = v->edges[i];
			T2 *ptr_e_it = &edge_array[id];
			if (ptr_e_it->weight < min_w)
			{
				ptr_e = ptr_e_it;
				min_w = ptr_e_it->weight;
			}
		}
		return ptr_e;
	}

	void merge_vertices(Edge *minmum_arc)
	{
		Vertex *seed_v = minmum_arc->left_vertex;
		Vertex *grow_v = minmum_arc->right_vertex;

		// step(2) update features
		// seed_v->mean = (seed_v->area * seed_v->mean + grow_v->area * grow_v->mean) / (seed_v->area + grow_v->area);
		// seed_v->area = seed_v->area + grow_v->area;
		this->update_features_func((T_Feature *)seed_v->ptr_features, (T_Feature *)grow_v->ptr_features, (T2 *)minmum_arc);

		// grow_v->node_type = -1;
		// grow_v->state = State::Freeze;

		// step(1) merge vi and vj
		seed_v->add_next_vertex(grow_v);
		// seed_v->is_changed = true;

		// step(3) remove eij both in Li and Lj
		seed_v->remove_edge(minmum_arc->ID);
		seed_v->state = State::Suspend; // a singnal for change of seed vertex
		grow_v->remove_edge(minmum_arc->ID);
		grow_v->state = State::Freeze;
		minmum_arc->state = State::Freeze;

		// step(4) replace grow_v vj to vi in the arcs pointed by Lj
		for (size_t i = 0; i < grow_v->edges.size(); ++i)
		{
			// ʹ��myVector[i]ִ�в���
			int id = grow_v->edges[i];
			Edge *e_in_grow_v = &this->edges[id];
			if (e_in_grow_v->left_vertex == grow_v)
			{
				e_in_grow_v->left_vertex = seed_v;
			}
			else if (e_in_grow_v->right_vertex == grow_v)
			{
				e_in_grow_v->right_vertex = seed_v;
			}
			else
			{
				std::cout << "Error edge id: " << e_in_grow_v->ID << std::endl;
				throw std::runtime_error("Topological relationship error!");
			}
		}

		// step(5) append Lj to Li, and remove the redundant arcs connecting the same adjacent nodes in Li and update length of edges
		for (size_t i = 0; i < grow_v->edges.size(); ++i)
		{
			int id = grow_v->edges[i];
			T2 *e_in_grow_v = &this->edges[id];
			std::pair<T2 *, T1 *> e_v = this->is_edge_connected_vertex(e_in_grow_v, (T1 *)seed_v);
			if (e_v.first != nullptr)
			{
				Edge *has_connected_edge = e_v.first;
				has_connected_edge->length += e_in_grow_v->length; // update the length of the edge
				has_connected_edge->weight = has_connected_edge->weight <= e_in_grow_v->weight ? has_connected_edge->weight : e_in_grow_v->weight;
				Vertex *neighbor_vertex = e_v.second;
				neighbor_vertex->remove_edge(e_in_grow_v->ID);
				// e_in_grow_v->type = -1;
				e_in_grow_v->state = State::Freeze;
				if (neighbor_vertex->direct == e_in_grow_v) // �ڽӽڵ�ָ�򱻺ϲ��Ľڵ�ʱ����Ϊgrow��seed�γ�cycle������e_in_grow_v�ض�����cycle
				{
					neighbor_vertex->change_direct(e_v.first, this->nng); // ������ָ�򱻺ϲ��ڵ��Ϊָ���½ڵ㣬��Ϊ�½ڵ㻹δ���跽�����Դ˴�����һ���������cycle
				}
			}
			else
			{
				seed_v->push(e_in_grow_v);
			}
		}
		// delete grow_v->edges
		grow_v->edges.clear();
		std::vector<int>().swap(grow_v->edges); // free memory in  grow_v->edges

		// step(6) update the weights of arcs connected to the merged node

		for (size_t i = 0; i < seed_v->edges.size(); ++i)
		{
			int id = seed_v->edges[i];
			Edge *e_in_seed_v = &this->edges[id];
			Vertex *left = e_in_seed_v->left_vertex;
			Vertex *right = e_in_seed_v->right_vertex;
			// double weight = left->area * right->area / (left->area + right->area) * pow(left->mean - right->mean, 2.0);
			double weight = this->merging_criteria((T_Feature *)left->ptr_features, (T_Feature *)right->ptr_features, (T2 *)e_in_seed_v);
			e_in_seed_v->weight = weight;
		}

		this->initial_direct((T1 *)seed_v, this->edges); // ���Ľڵ��ʼ���������ܲ����µ�cycle����������Ҫ���ж��Ƿ����nng

		if (seed_v->direct->shortest_arc_num == 2)
			this->nng->push(seed_v->direct);

		// step(7)update the shortest arc and cycles of adjacent nodes����ʱ��RAG��edge�γɵ����е�cycle������nng�У���������ֻ�ø����ڽӽڵ�ķ���
		for (size_t i = 0; i < seed_v->edges.size(); ++i)
		{
			int id = seed_v->edges[i];
			T2 *e_in_seed_v = &this->edges[id];
			T1 *neighbouring_v = (T1 *)e_in_seed_v->get_other_vertex(seed_v); // �ڽӽڵ�
			this->update_direction_cycle_in_2nd_degree_vertices(neighbouring_v, this->nng);
		}
	}

	void merge_vertices(Edge *minmum_arc, NNG *region_nng)
	{
		Vertex *seed_v = minmum_arc->left_vertex;
		Vertex *grow_v = minmum_arc->right_vertex;

		if (seed_v == grow_v)
		{
			std::cout << minmum_arc->ID << std::endl;
			throw std::runtime_error("Error: try to merge the same vertex!");
		}

		// update centroid positions of seed_v
		for (int i = 0; i < seed_v->position_dimension; i++)
		{
			seed_v->positions[i] = (seed_v->positions[i] * seed_v->size + grow_v->positions[i] * grow_v->size) / (seed_v->size + grow_v->size);
		}

		// step(2) update features
		this->update_features_func((T_Feature *)seed_v->ptr_features, (T_Feature *)grow_v->ptr_features, (T2 *)minmum_arc);

		// step(1) merge vi and vj
		seed_v->add_next_vertex(grow_v);
		// seed_v->is_changed = true;

		// step(3) remove eij both in Li and Lj
		seed_v->remove_edge(minmum_arc->ID);
		seed_v->state = State::Suspend; // a singnal for the change of seed vertex
		grow_v->remove_edge(minmum_arc->ID);
		grow_v->state = State::Freeze;
		minmum_arc->state = State::Freeze;

		// step(4) replace grow_v vj to vi in the arcs pointed by Lj
		for (size_t i = 0; i < grow_v->edges.size(); ++i)
		{
			// ʹ��myVector[i]ִ�в���
			// ���磬��ӡԪ��ֵ
			int id = grow_v->edges[i];
			Edge *e_in_grow_v = &this->edges[id];
			if (e_in_grow_v->left_vertex == grow_v)
			{
				e_in_grow_v->left_vertex = seed_v;
			}
			else if (e_in_grow_v->right_vertex == grow_v)
			{
				e_in_grow_v->right_vertex = seed_v;
			}
			else
			{
				std::cout << e_in_grow_v->left_vertex->ID << " " << e_in_grow_v->right_vertex->ID << std::endl;
				std::cout << seed_v->ID << " " << grow_v->ID << std::endl;
				std::cout << "mini edge: " << minmum_arc->ID;
				std::cout << "error edge: " << e_in_grow_v->ID;
				throw std::runtime_error("Topological relationship error!");
			}
		}

		// step(5) append Lj to Li, and remove the redundant arcs connecting the same adjacent nodes in Li and update length of edges
		for (size_t i = 0; i < grow_v->edges.size(); ++i)
		{
			int id = grow_v->edges[i];
			T2 *e_in_grow_v = &this->edges[id];
			std::pair<T2 *, T1 *> e_v = this->is_edge_connected_vertex(e_in_grow_v, (T1 *)seed_v);
			if (e_v.first != nullptr)
			{
				Edge *has_connected_edge = e_v.first;
				has_connected_edge->length += e_in_grow_v->length; // update the length of the edge
				has_connected_edge->weight = has_connected_edge->weight <= e_in_grow_v->weight ? has_connected_edge->weight : e_in_grow_v->weight;
				Vertex *neighbor_vertex = e_v.second;
				neighbor_vertex->remove_edge(e_in_grow_v->ID);
				// e_in_grow_v->type = -1;
				e_in_grow_v->state = State::Freeze;

				if (neighbor_vertex->direct == e_in_grow_v) // �ڽӽڵ�ָ�򱻺ϲ��Ľڵ�ʱ����Ϊgrow��seed�γ�cycle������e_in_grow_v�ض�����cycle
				{
					neighbor_vertex->change_direct(e_v.first, region_nng); // ������ָ�򱻺ϲ��ڵ��Ϊָ���½ڵ㣬��Ϊ�½ڵ㻹δ���跽�����Դ˴�����һ���������cycle
				}
			}
			else
			{
				seed_v->push(e_in_grow_v);
			}
		}
		// delete grow_v->edges
		grow_v->edges.clear();
		std::vector<int>().swap(grow_v->edges); // free memory in  grow_v->edges

		// step(6) update the weights of arcs connected to the merged node

		for (size_t i = 0; i < seed_v->edges.size(); ++i)
		{
			int id = seed_v->edges[i];
			Edge *e_in_seed_v = &this->edges[id];
			Vertex *left = e_in_seed_v->left_vertex;
			Vertex *right = e_in_seed_v->right_vertex;
			// double weight = left->area * right->area / (left->area + right->area) * pow(left->mean - right->mean, 2.0);
			double weight = this->merging_criteria((T_Feature *)left->ptr_features, (T_Feature *)right->ptr_features, (T2 *)e_in_seed_v);
			e_in_seed_v->weight = weight;
		}

		this->initial_direct((T1 *)seed_v, this->edges); // ���Ľڵ��ʼ���������ܲ����µ�cycle����������Ҫ���ж��Ƿ����nng

		if (seed_v->direct->shortest_arc_num == 2)
		{
			region_nng->push(seed_v->direct);
		}

		// step(7)update the shortest arc and cycles of adjacent nodes����ʱ��RAG��edge�γɵ����е�cycle������nng�У���������ֻ�ø����ڽӽڵ�ķ���
		for (size_t i = 0; i < seed_v->edges.size(); ++i)
		{
			int id = seed_v->edges[i];
			T2 *e_in_seed_v = &this->edges[id];
			T1 *neighbouring_v = (T1 *)e_in_seed_v->get_other_vertex(seed_v); // �ڽӽڵ�
			this->update_direction_cycle_in_2nd_degree_vertices(neighbouring_v, region_nng);
		}
	}

	void simply_merge_vertices(Edge *minmum_arc)
	{
		Vertex *seed_v = minmum_arc->left_vertex;
		Vertex *grow_v = minmum_arc->right_vertex;

		if (seed_v == grow_v)
		{
			std::cout << minmum_arc->ID << std::endl;
			throw std::runtime_error("Error: try to merge the same vertex!");
		}

		// only merge two vertice, important, easy to be wrong
		// update centroid positions of seed_v
		for (int i = 0; i < seed_v->position_dimension; i++)
		{
			seed_v->positions[i] = (seed_v->positions[i] * seed_v->size + grow_v->positions[i] * grow_v->size) / (seed_v->size + grow_v->size);
		}

		// step(2) update features
		this->update_features_func((T_Feature *)seed_v->ptr_features, (T_Feature *)grow_v->ptr_features, (T2 *)minmum_arc);

		// step(1) merge vi and vj
		seed_v->add_next_vertex(grow_v);

		// step(3) remove eij both in Li and Lj
		seed_v->remove_edge(minmum_arc->ID); // Test
		seed_v->state = State::Suspend;		 // a singnal for the change of seed vertex
		grow_v->remove_edge(minmum_arc->ID); // Test
		grow_v->state = State::Freeze;
		minmum_arc->state = State::Freeze;

		// step(4) change the edge state of seed and grow vertex
		for (size_t i = 0; i < seed_v->edges.size(); ++i) // edge state of seed vertex
		{
			int id = seed_v->edges[i];
			Edge *e = &this->edges[id];
			e->state = State::Suspend;
		}
		for (size_t i = 0; i < grow_v->edges.size(); ++i) // edge state of grow vertex
		{
			int id = grow_v->edges[i];
			Edge *e = &this->edges[id];
			e->state = State::Suspend;
		}
	}

	T1 *global_best_optimization(double threshold)
	{
		// double las_min = 0.0;
		this->construct_nng();
		if (!this->nng->empty())
		{
			Edge *minmum_arc = this->nng->pop();
			// las_min = minmum_arc->weight;
			while (minmum_arc->weight < threshold)
			{
				this->merge_vertices(minmum_arc);
				if (this->nng->empty())
					break;
				minmum_arc = this->nng->pop();
			}
			if (!this->nng->empty())
			{
				this->nng->push(minmum_arc);
			}
		}
		if (this->nng->empty())
		{
			std::cout << "NNG empty!" << std::endl;
		}
		else
		{
			Edge *minmum_arc = this->nng->pop();
			std::cout << "Last weight" << minmum_arc->weight << std::endl;
		}
		// std::cout << "Last edge weight: " << las_min << std::endl;
		int activate_v = 0;
		for (int i = 0; i < this->vertice_num; i++)
		{
			Vertex *v = (Vertex *)&this->vertice[i];
			if (v->state != State::Freeze)
			{
				activate_v++;
			}
		}
		std::cout << "Ouput vertex num: " << activate_v << std::endl;
		std::pair<T1 *, int> chains = this->chaining_vertice();
		return chains.first;
	}

	void global_best_optimization(NNG *ptr_nng, double threshold)
	{
		if (!ptr_nng->empty())
		{
			Edge *minmum_arc = ptr_nng->pop();
			while (minmum_arc->weight < threshold)
			{
				this->merge_vertices(minmum_arc, ptr_nng);
				if (ptr_nng->empty())
					break;
				minmum_arc = ptr_nng->pop();
			}
			if (!ptr_nng->empty())
			{
				ptr_nng->push(minmum_arc);
			}
		}
	}

	T1 *purning_optimization(double threshold)
	{
		int iter = 4;
		double *scales = new double[iter]{threshold * 0.001, threshold * 0.01, threshold * 0.1, threshold};

		for (int i = 0; i < iter; i++)
		{
			std::cout << "Iteration: " << GREEN_TEXT << i << RESET_COLOR << " scale: " << scales[i] << " is working..................................................."
					  << std::endl;
			this->freeze_edges(threshold);
			this->suspend_edges(scales[i]);
			this->global_best_optimization(scales[i]);
			this->activate_edges(State::Suspend);
		}

		delete[] scales;

		std::pair<T1 *, int> chains = this->chaining_vertice();
		return chains.first;
	}

	Edge *region_construct_nng(Vertex *root, NNG *region_nng)
	{
		T1 *v = (T1 *)root;
		if (root == nullptr)
			throw std::runtime_error("Null Vertex!");
		auto start_time = std::chrono::high_resolution_clock::now();
		T2 *e_root = nullptr;
		T2 *e_nil = nullptr;
		while (v != nullptr)
		{
			for (auto it = v->edges.begin(); it != v->edges.end(); ++it)
			{
				T2 *e = &this->edges[*it];
				if (e_root == nullptr)
				{
					e->is_checked = true;
					e_root = e;
					e_nil = e_root;
				}
				if (e_root != e && !e->is_checked)
				{
					e->is_checked = true;
					e_nil->next_edge_in_region = e;
					e_nil = e;
				}
				e->shortest_arc_num = 0;
			}
			v = (T1 *)v->next_in_region;
		}
		if (e_nil != nullptr)
			e_nil->next_edge_in_region = nullptr;
		v = (T1 *)root;
		while (v != nullptr)
		{
			this->initial_direct((T1 *)v, this->edges);
			v = (T1 *)v->next_in_region;
		}

		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
		e_nil = e_root;
		while (e_nil != nullptr)
		{
			if (e_nil->state == State::Activate && e_nil->shortest_arc_num == 2)
			{
				region_nng->push(e_nil);
			}
			e_nil = (T2 *)e_nil->next_edge_in_region;
		}

		start_time = std::chrono::high_resolution_clock::now();
		end_time = std::chrono::high_resolution_clock::now();
		duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
		return e_root;
	}

	void region_global_best_optimization(std::vector<std::pair<T1 *, int>> results, double threshold)
	{
		for (auto it = results.begin(); it != results.end(); ++it)
		{
			Vertex *root = it->first;
			int value = it->second;
			NNG *region_nng = new NNG();
			// process vertex and value
			if (root != nullptr)
			{
				this->region_construct_nng(root, region_nng);
			}
			// this->check_direct(root);

			if (!region_nng->empty())
			{
				Edge *minmum_arc = region_nng->pop();

				while (minmum_arc->weight < threshold)
				{
					this->merge_vertices(minmum_arc, region_nng);

					if (region_nng->empty())
						break;
					minmum_arc = region_nng->pop();
				}
			}

			delete region_nng;
		}
	}

	int regionalize(std::vector<std::pair<T1 *, int>> results)
	{
		int suspend_edge_num = 0;
		for (auto it = results.begin(); it != results.end(); ++it)
		{
			T1 *v = it->first;
			T1 *other_v = nullptr;
			v->direct = nullptr;
			// process vertex and related edges
			while (v != nullptr)
			{
				for (auto it = v->edges.begin(); it != v->edges.end();)
				{

					Edge *e = &this->edges[*it];
					T1 *left_vertex = (T1 *)e->left_vertex;
					T1 *right_vertex = (T1 *)e->right_vertex;
					e->shortest_arc_num = 0;
					if (left_vertex == v)
						other_v = right_vertex;
					else
						other_v = left_vertex;

					if ((left_vertex->bits_num == right_vertex->bits_num && left_vertex->bit_sign != right_vertex->bit_sign) ||
						left_vertex->bits_num != right_vertex->bits_num)
					{
						other_v->remove_edge(*it);
						it = v->edges.erase(it); // ʹ�÷���ֵ���µ�����
						e->state = State::Suspend;
						suspend_edge_num++;
					}
					else
					{
						++it; // ֻ�е�û��ɾ��Ԫ��ʱ�ŵ���������
					}
				}
				v->edges.shrink_to_fit(); // ��ѭ��������һ���������ڴ�
				v = (T1 *)v->next_in_region;
			}
		}
		return suspend_edge_num;
	}

	T1 *region_optimization(double threshold, int region_size = 2500)
	{
		int iter = 5;
		T1 *root = &this->vertice[0];
		double *scales = new double[iter]{threshold * 0.001, threshold * 0.01, threshold * 0.1, threshold, threshold};
		// double* scales = new double[iter] {  threshold * 0.001, threshold * 0.01, threshold };
		std::pair<T1 *, int> chains;
		for (int i = 0; i < iter; i++)
		{
			std::cout << "Iter " << i + 1 << ":" << std::endl;
			// std::vector<std::pair<Vertex*, int>> results = this->partition_dataset_by_kd_tree(root, 40000, 10);
			std::vector<std::pair<T1 *, int>> results;
			if (i < iter - 1)
			{
				std::cout << "[INFO] Calling func\"partition_dataset_by_kd_tree\" (member of class RAG) ..." << std::endl;
				results = this->partition_dataset_by_kd_tree(root, region_size);
			}
			else
			{
				results.push_back(chains);
			}

			std::cout << "[INFO] Calling func\"regionalize\" ..." << std::endl;
			this->regionalize(results); // regionalize the vertice
			// int cur_size = this->check_data_size_in_region();
			// std::cout << "Data size after func regionalize(): " << cur_size << std::endl;
			// std::cout << "[INFO] Regionalize Done! " << std::endl;

			std::cout << "[INFO] Calling func \"region_global_best_optimization(params: std::vector<std::pair<T1 *, int>>results, scales[i])\"" << std::endl;
			this->region_global_best_optimization(results, scales[i]); // global optimisation in each region

			std::cout << "[INFO] Iter: " << i + 1 << "'s opt has done! " << std::endl;

			// cur_size = this->check_data_size_in_output();
			// std::cout << "Data size after func optimization(): " << cur_size << std::endl;
			std::cout << "[INFO] Activating suspended edges ... " << std::endl;
			this->activate_edges(State::Suspend);
			// ����Ҫ�����ڼ���״̬��vertex���´�����

			std::cout << "[INFO] Chaining vertices ... " << std::endl;
			chains = this->chaining_vertice();
			root = chains.first;
			int num = chains.second;
			// std::cout << "Data size after chaining_vertice: " << num << std::endl;
			int cur_num = 0;
			Vertex *nil = root;
			while (nil != nullptr)
			{
				nil = nil->next_in_region;
				cur_num++;
			}
			std::cout << "Num: " << cur_num << std::endl;
		}

		std::cout << "End merge" << std::endl;
		delete[] scales;
		return root;
	}

	T1 *cycle_optimization(double threshold)
	{
		int iter_num = 0;
		bool stop = false;

		int num_less_10000 = 0;
		long long merge_time = 0;
		long long nng_time = 0;
		long long activate_time = 0;
		long long update_time = 0;
		while (!stop)
		{
			// construct nng
			auto nng_start_time = std::chrono::high_resolution_clock::now();

			NNG *region_nng = new NNG();
			this->construct_nng(region_nng);
			auto nng_end_time = std::chrono::high_resolution_clock::now();
			auto nng_duration = std::chrono::duration_cast<std::chrono::milliseconds>(nng_end_time - nng_start_time);
			long long nng_time_iter = nng_duration.count();
			nng_time += nng_time_iter;

			auto merge_start_time = std::chrono::high_resolution_clock::now();
			std::pair<bool, int> result = this->cycle_best_optimisation(threshold, region_nng);
			auto merge_end_time = std::chrono::high_resolution_clock::now();
			auto merge_duration = std::chrono::duration_cast<std::chrono::milliseconds>(merge_end_time - merge_start_time);
			long long merge_time_iter = merge_duration.count();
			merge_time += merge_time_iter;

			stop = result.first;
			int merge_num = result.second;
			// std::cout << "Merge num: " << merge_num << std::endl;
			int objects = this->check_data_size_in_region();
			std::cout << "Obj num: " << GREEN_TEXT << objects << RESET_COLOR << " Merge num: " << GREEN_TEXT << merge_num << RESET_COLOR << " ";
			// std::cout << "Obj + Merge: " << objects + merge_num << std::endl;
			// activate edge

			auto activate_start_time = std::chrono::high_resolution_clock::now();
			this->activate_edges(State::Suspend, "123");
			auto activate_end_time = std::chrono::high_resolution_clock::now();
			auto activate_duration = std::chrono::duration_cast<std::chrono::milliseconds>(activate_end_time - activate_start_time);
			long long activate_time_iter = activate_duration.count();
			activate_time += activate_time_iter;

			auto update_start_time = std::chrono::high_resolution_clock::now();
			this->update_weights(State::Suspend);
			auto update_end_time = std::chrono::high_resolution_clock::now();
			auto update_duration = std::chrono::duration_cast<std::chrono::milliseconds>(update_end_time - update_start_time);
			long long update_time_iter = update_duration.count();
			update_time += update_time_iter;

			auto iter_time = std::chrono::duration_cast<std::chrono::milliseconds>(update_end_time - nng_start_time);
			double seconds = iter_time.count() / 1000;
			if (seconds == 0)
				seconds = 1;
			int merge_fre = merge_num / seconds;
			std::cout << "NNG time: " << GREEN_TEXT << nng_time_iter << RESET_COLOR << "ms"
					  << " Merge time: " << GREEN_TEXT << merge_time_iter << RESET_COLOR << "ms"
					  << " Edge time: " << GREEN_TEXT << activate_time_iter << RESET_COLOR << "ms"
					  << " Update time: " << GREEN_TEXT << update_time_iter << RESET_COLOR << "ms"
					  << " Iter time: " << GREEN_TEXT << iter_time.count() << RESET_COLOR << "s"
					  << " Merge fre: " << GREEN_TEXT << merge_fre << RESET_COLOR << " thousand pairs/s" << std::endl;

			/*if (merge_num < 25000)
			{
				num_less_10000 += merge_num;
				merge_num = 0;
				NNG* region_nng11 = new NNG();

				this->construct_nng(region_nng11);
				Edge* minmum_arc;
				if (!region_nng11->empty())
				{
					minmum_arc = region_nng11->pop();
					while (minmum_arc->weight < threshold)
					{
						this->merge_vertices(minmum_arc, region_nng11);
						merge_num++;
						if (region_nng11->empty()) break;
						minmum_arc = region_nng11->pop();
					}
					std::cout << "minmum_arc weight: " << minmum_arc->weight << std::endl;
				}
				std::cout << "Merge num: " << merge_num << std::endl;
				stop = true;
				delete region_nng11;
			}*/
			delete region_nng;
			iter_num++;
		}
		// std::cout << "num_less_10000 num: " << num_less_10000 << std::endl;
		// std::cout << "End merge: iter mun = "<< iter_num << std::endl;

		std::cout << "Total NNG time: " << GREEN_TEXT << nng_time / 1000 << RESET_COLOR << "s" << std::endl
				  << "Total merge time: " << GREEN_TEXT << merge_time / 1000 << RESET_COLOR << "s" << std::endl
				  << "Total edge time: " << GREEN_TEXT << activate_time / 1000 << RESET_COLOR << "s" << std::endl
				  << "Total update time: " << GREEN_TEXT << update_time / 1000 << RESET_COLOR << "s" << std::endl;

		std::pair<T1 *, int> chains = this->chaining_vertice();
		return chains.first;
	}

	T1 *region_pruning_optimization(double threshold, int region_size)
	{
		int iter = 5;
		T1 *root = &this->vertice[0];
		// double* scales = new double[iter] { threshold * 0.001, threshold * 0.01, threshold * 0.1, threshold };
		double *scales = new double[iter]{threshold * 0.001, threshold * 0.01, threshold * 0.1, threshold, threshold};
		std::pair<T1 *, int> chains;
		for (int i = 0; i < iter; i++)
		{
			std::cout << "Iter " << i + 1 << ":" << std::endl;
			// std::vector<std::pair<Vertex*, int>> results = this->partition_dataset_by_kd_tree(root, 40000, 10);
			std::vector<std::pair<T1 *, int>> results;
			if (threshold == scales[i])
				this->freeze_edges(threshold);
			if (i < iter - 1)
			{
				// results = this->partition_dataset_by_kd_tree(root, 40000, 10);
				results = this->partition_dataset_by_kd_tree(root, region_size);
			}
			else
			{
				results.push_back(chains);
			}

			this->regionalize(results); // regionalize the vertice
			// int cur_size = this->check_data_size_in_region();
			// std::cout << "Data size after func regionalize(): " << cur_size << std::endl;
			this->suspend_edges(scales[i]);

			this->region_global_best_optimization(results, scales[i]); // global optimisation in each region
			// cur_size = this->check_data_size_in_output();
			// std::cout << "Data size after func optimization(): " << cur_size << std::endl;

			this->activate_edges(State::Suspend);
			// ����Ҫ�����ڼ���״̬��vertex���´�����
			chains = this->chaining_vertice();
			root = chains.first;
			int num = chains.second;
			// std::cout << "Data size after chaining_vertice: " << num << std::endl;
		}

		delete[] scales;
		return root;
	}

	T1 *region_cycle_optimization(double threshold, int region_size)
	{
		int iter = 2;
		T1 *root = &this->vertice[0];
		double *scales = new double[iter]{threshold * 0.1, threshold};
		std::pair<T1 *, int> chains;
		for (int i = 0; i < 2; i++)
		{
			std::cout << "Iter " << i + 1 << ":" << std::endl;
			std::vector<std::pair<T1 *, int>> results;
			results = this->partition_dataset_by_kd_tree(root, region_size);
			std::cout << "Set num: " << results.size() << std::endl;
			int suspend_edge_num = this->regionalize(results); // regionalize the vertice

			this->region_cycle_best_optimization(results, scales[i]); // global optimisation in each region
			// std::cout << "Regionalize edge num: " << suspend_edge_num << std::endl;
			this->activate_edges(State::Suspend);

			// ����Ҫ�����ڼ���״̬��vertex���´�����
			chains = this->chaining_vertice();
			root = chains.first;
			int num = chains.second;
		}
		delete[] scales;
		return root;
	}

	void region_cycle_best_optimization(std::vector<std::pair<T1 *, int>> results, double threshold, bool is_freeze = false)
	{
		int idx = 1;
		for (auto it = results.begin(); it != results.end(); ++it)
		{
			int total_merge_num = 0;
			auto start_time = std::chrono::high_resolution_clock::now();
			Vertex *v_root = it->first;
			bool stop = false;
			while (!stop)
			{
				// construct nng
				if (!is_freeze)
					this->freeze_edges(threshold);
				NNG *region_nng = new NNG();
				Edge *e_root = this->region_construct_nng(v_root, region_nng);
				int obj_num = this->check_data_size_in_region(v_root);
				std::pair<bool, int> result = this->cycle_best_optimisation(threshold, region_nng);
				stop = result.first;
				int merge_num = result.second;
				total_merge_num += merge_num;
				e_root = this->activate_edges_in_region(e_root);
				this->update_weight_in_region(e_root, State::Suspend);
				v_root = this->update_root_in_region(v_root);
				int obj_new_num = this->check_data_size_in_region(v_root);
				std::cout << "Obj num: " << obj_num << " Merge: " << merge_num
						  << " New num: " << obj_new_num
						  << " Plus: " << merge_num + obj_new_num
						  << std::endl;

				if (merge_num < 20000)
				{
					merge_num = 0;
					NNG *region_nng11 = new NNG();

					this->region_construct_nng(v_root, region_nng11);

					Edge *minmum_arc;
					if (!region_nng11->empty())
					{
						minmum_arc = region_nng11->pop();
						while (minmum_arc->weight < threshold)
						{
							this->merge_vertices(minmum_arc, region_nng11);
							merge_num++;
							if (region_nng11->empty())
								break;
							minmum_arc = region_nng11->pop();
						}
						// std::cout << "minmum_arc weight: " << minmum_arc->weight << std::endl;
					}
					// std::cout << "Merge num: " << merge_num << std::endl;
					total_merge_num += merge_num;
					stop = true;
					delete region_nng11;
				}
				// std::cout << "Merge num: " << merge_num << std::endl;
				delete region_nng;
			}
			auto end_time = std::chrono::high_resolution_clock::now();
			auto activate_duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
			double activate_time_iter = activate_duration.count() / 1000.0;
			std::cout << "Set " << GREEN_TEXT << idx << RESET_COLOR << " Total num: " << GREEN_TEXT << total_merge_num << RESET_COLOR << " Total time: " << GREEN_TEXT << activate_time_iter << RESET_COLOR << "s" << std::endl;
			idx++;
		}
	}

	void process_unit_region_cycle_best_parallel_optimization(std::pair<T1 *, int> &result, int idx, double threshold, bool is_freeze)
	{
		int total_merge_num = 0;
		Vertex *v_root = result.first;
		bool stop = false;
		double mini_W = 0.0;
		Edge *e_root_test = nullptr;
		while (!stop)
		{
			// construct nng
			if (!is_freeze)
				this->freeze_edges(threshold);
			NNG *region_nng = new NNG();
			Edge *e_root = this->region_construct_nng(v_root, region_nng);
			int obj_num = this->check_data_size_in_region(v_root);
			std::pair<bool, int> result = this->cycle_best_optimisation(threshold, region_nng);
			stop = result.first;
			int merge_num = result.second;
			total_merge_num += merge_num;
			e_root = this->activate_edges_in_region(e_root); // extrac edge in activate.
			this->update_weight_in_region(e_root, State::Suspend);
			e_root_test = e_root;
			v_root = this->update_root_in_region(v_root);
			int obj_new_num = this->check_data_size_in_region(v_root);
			// std::cout << "Obj num: " << obj_num << " Merge: " << merge_num
			//	<< " New num: " << obj_new_num
			//	<< " Plus: " << merge_num + obj_new_num
			//	<< std::endl;
			if (merge_num < 2000)
			{
				merge_num = 0;
				NNG *region_nng11 = new NNG();
				this->region_construct_nng(v_root, region_nng11);

				Edge *minmum_arc;
				if (!region_nng11->empty())
				{
					minmum_arc = region_nng11->pop();
					while (minmum_arc->weight < threshold)
					{
						this->merge_vertices(minmum_arc, region_nng11);
						merge_num++;
						mini_W = minmum_arc->weight;
						if (region_nng11->empty())
							break;
						minmum_arc = region_nng11->pop();
					}
					// std::cout << "minmum_arc weight: " << minmum_arc->weight << std::endl;
				}
				// std::cout << "Merge num: " << merge_num << std::endl;
				e_root = this->activate_edges_in_region(e_root); // extrac edge in activate.
				this->update_weight_in_region(e_root, State::Suspend);
				e_root_test = e_root;
				total_merge_num += merge_num;
				stop = true;
				delete region_nng11;
			}
			// std::cout << "Merge num: " << merge_num << std::endl;
			delete region_nng;
		}
		std::cout << "Current set " << GREEN_TEXT << idx << RESET_COLOR << "; minmum weight: " << GREEN_TEXT << mini_W << RESET_COLOR << std::endl;

		// int e_root_test_suspend_num = 0;
		// while (e_root_test != nullptr)
		//{
		//	if (e_root_test->state == State::Suspend)
		//	{
		//		e_root_test_suspend_num++;
		//	}
		//	e_root_test = e_root_test->next_edge_in_region;
		// }
		// std::cout << "Suspend_num " << GREEN_TEXT << e_root_test_suspend_num << RESET_COLOR << std::endl;
	}

	void region_cycle_best_parallel_optimization(std::vector<std::pair<T1 *, int>> results, double threshold, bool is_freeze = false)
	{
		// ��ȡӲ��֧�ֵ�����߳���
		unsigned int max_threads = std::thread::hardware_concurrency();
		// �����̳߳�
		std::vector<std::future<void>> futures;
		vector<T2 *> roots;
		// �����߳�
		int idx = 1;
		for (auto it = results.begin(); it != results.end(); ++it, ++idx)
		{

			// �����߳�����������Ӳ��֧�ֵ�����߳���
			if (futures.size() >= max_threads)
			{
				// �ȴ�һ���߳����
				futures.front().get();
				futures.erase(futures.begin());
			}

			futures.push_back(std::async(std::launch::async,
										 [this, it, idx, threshold, is_freeze]()
										 {
											 this->process_unit_region_cycle_best_parallel_optimization(*it, idx, threshold, is_freeze);
										 }));
		}

		// �ȴ������߳����
		for (auto &future : futures)
		{
			future.get();
		}
	}

	T1 *pruning_cycle_optimization(double threshold)
	{
		int iter_num = 0;
		bool stop = false;
		int num_less_10000 = 0;
		long long merge_time = 0;
		long long nng_time = 0;
		long long activate_time = 0;
		long long update_time = 0;

		while (!stop)
		{
			// construct nng
			auto nng_start_time = std::chrono::high_resolution_clock::now();
			this->freeze_edges(threshold);
			NNG *region_nng = new NNG();
			this->construct_nng(region_nng);
			auto nng_end_time = std::chrono::high_resolution_clock::now();
			auto nng_duration = std::chrono::duration_cast<std::chrono::milliseconds>(nng_end_time - nng_start_time);
			long long nng_time_iter = nng_duration.count();
			nng_time += nng_time_iter;

			auto merge_start_time = std::chrono::high_resolution_clock::now();
			std::pair<bool, int> result = this->cycle_best_optimisation(threshold, region_nng);
			auto merge_end_time = std::chrono::high_resolution_clock::now();
			auto merge_duration = std::chrono::duration_cast<std::chrono::milliseconds>(merge_end_time - merge_start_time);
			long long merge_time_iter = merge_duration.count();
			merge_time += merge_time_iter;

			stop = result.first;
			int merge_num = result.second;
			// std::cout << "Merge num: " << merge_num << std::endl;
			int objects = this->check_data_size_in_region();
			std::cout << "Obj num: " << GREEN_TEXT << objects << RESET_COLOR << " Merge num: " << GREEN_TEXT << merge_num << RESET_COLOR << " ";
			// std::cout << "Obj + Merge: " << objects + merge_num << std::endl;
			// activate edge

			auto activate_start_time = std::chrono::high_resolution_clock::now();
			this->activate_edges(State::Suspend, "123");
			auto activate_end_time = std::chrono::high_resolution_clock::now();
			auto activate_duration = std::chrono::duration_cast<std::chrono::milliseconds>(activate_end_time - activate_start_time);
			long long activate_time_iter = activate_duration.count();
			activate_time += activate_time_iter;

			auto update_start_time = std::chrono::high_resolution_clock::now();
			this->update_weights(State::Suspend);
			auto update_end_time = std::chrono::high_resolution_clock::now();
			auto update_duration = std::chrono::duration_cast<std::chrono::milliseconds>(update_end_time - update_start_time);
			long long update_time_iter = update_duration.count();
			update_time += update_time_iter;

			auto iter_time = std::chrono::duration_cast<std::chrono::milliseconds>(update_end_time - nng_start_time);
			double seconds = iter_time.count() / 1000;
			if (seconds == 0)
				seconds = 1;
			int merge_fre = merge_num / seconds;
			std::cout << "NNG time: " << GREEN_TEXT << nng_time_iter << RESET_COLOR << "ms"
					  << " Merge time: " << GREEN_TEXT << merge_time_iter << RESET_COLOR << "ms"
					  << " Edge time: " << GREEN_TEXT << activate_time_iter << RESET_COLOR << "ms"
					  << " Update time: " << GREEN_TEXT << update_time_iter << RESET_COLOR << "ms"
					  << " Iter time: " << GREEN_TEXT << iter_time.count() << RESET_COLOR << "s"
					  << " Merge fre: " << GREEN_TEXT << merge_fre << RESET_COLOR << " thousand pairs/s" << std::endl;
			delete region_nng;
			iter_num++;

			if (merge_num < 25000)
			{
				num_less_10000 += merge_num;
				merge_num = 0;
				NNG *region_nng11 = new NNG();

				this->construct_nng(region_nng11);
				Edge *minmum_arc;
				if (!region_nng11->empty())
				{
					minmum_arc = region_nng11->pop();
					while (minmum_arc->weight < threshold)
					{
						this->merge_vertices(minmum_arc, region_nng11);
						merge_num++;
						if (region_nng11->empty())
							break;
						minmum_arc = region_nng11->pop();
					}
					std::cout << "minmum_arc weight: " << minmum_arc->weight << std::endl;
				}
				std::cout << "Merge num: " << merge_num << std::endl;
				stop = true;
				delete region_nng11;
			}
		}

		std::cout << "Total NNG time: " << GREEN_TEXT << nng_time / 1000 << RESET_COLOR << "s" << std::endl
				  << "Total merge time: " << GREEN_TEXT << merge_time / 1000 << RESET_COLOR << "s" << std::endl
				  << "Total edge time: " << GREEN_TEXT << activate_time / 1000 << RESET_COLOR << "s" << std::endl
				  << "Total update time: " << GREEN_TEXT << update_time / 1000 << RESET_COLOR << "s" << std::endl;

		std::pair<T1 *, int> chains = this->chaining_vertice();
		return chains.first;
	}

	T1 *region_pruning_cycle_optimization(double threshold, int region_size)
	{
		int iter = 2;
		T1 *root = &this->vertice[0];
		double *scales = new double[iter]{threshold * 0.1, threshold};
		std::pair<T1 *, int> chains;
		for (int i = 0; i < 2; i++)
		{
			std::cout << "Iter " << i + 1 << ":" << std::endl;
			std::vector<std::pair<T1 *, int>> results;
			results = this->partition_dataset_by_kd_tree(root, region_size);
			// results = this->partition_dataset_by_kd_tree(root, 400, 10);

			std::cout << "Set num: " << results.size() << std::endl;
			int suspend_edge_num = this->regionalize(results);				// regionalize the vertice
			this->region_cycle_best_optimization(results, scales[i], true); // global optimisation in each region
			// std::cout << "Regionalize edge num: " << suspend_edge_num << std::endl;
			this->activate_edges(State::Suspend);

			// ����Ҫ�����ڼ���״̬��vertex���´�����
			chains = this->chaining_vertice();
			root = chains.first;
			int num = chains.second;
		}

		delete[] scales;
		return root;
	}

	T1 *region_pruning_cycle_optimization_MultiCores(double threshold, int region_size)
	{
		int iter = 3;
		T1 *root = &this->vertice[0];
		double *scales = new double[iter]{threshold * 0.1, threshold * 0.2, threshold};
		std::pair<T1 *, int> chains;
		long long total_divide_time = 0;
		for (int i = 0; i < iter; i++)
		{
			std::cout << "Iter " << GREEN_TEXT << i + 1 << RESET_COLOR << ":" << std::endl;
			std::vector<std::pair<T1 *, int>> results;
			// results = this->partition_dataset_by_kd_tree(root, 250000);
			// auto start_time = std::chrono::high_resolution_clock::now();
			results = this->partition_dataset_by_kd_tree(root, region_size);
			// auto end_time = std::chrono::high_resolution_clock::now();
			// auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
			// long long t = duration.count();
			// total_divide_time += t;
			std::cout << "Set num: " << GREEN_TEXT << results.size() << RESET_COLOR << std::endl;
			int count = 0;
			std::cout << "Set tables: " << std::endl
					  << "\t" << GREEN_TEXT;
			// GREEN_TEXT << idx << RESET_COLOR
			for (const auto &pair : results)
			{
				T1 *vertex = pair.first;
				int value = pair.second;
				std::cout << value << " ";
				if ((count + 1) % 5 == 0)
				{
					std::cout << std::endl
							  << "\t";
				}
				count++;
			}
			std::cout << std::endl
					  << RESET_COLOR;
			int suspend_edge_num = this->regionalize(results);						 // regionalize the vertice
			this->region_cycle_best_parallel_optimization(results, scales[i], true); // global optimisation in each region
			this->activate_edges(State::Suspend);

			// ����Ҫ�����ڼ���״̬��vertex���´�����
			chains = this->chaining_vertice();
			root = chains.first;
			int num = chains.second;
		}

		// std::cout << "Return root" << std::endl;
		std::cout << "out num: " << this->check_data_size_in_region((Vertex *)root) << std::endl;
		// std::cout << "Divide time: " << total_divide_time << std::endl;

		return root;
	}

	Vertex *update_root_in_region(Vertex *org_root)
	{
		Vertex *new_root = nullptr;
		Vertex *nil = nullptr;
		while (org_root != nullptr)
		{
			if (org_root->state != State::Freeze && new_root == nullptr)
			{
				new_root = org_root;
				nil = new_root;
				org_root = org_root->next_in_region;
				continue;
			}
			if (org_root->state != State::Freeze)
			{
				nil->next_in_region = org_root;
				nil = org_root;
			}
			org_root = org_root->next_in_region;
		}
		if (nil != nullptr)
			nil->next_in_region = nullptr;
		return new_root;
	}

	std::pair<bool, int> cycle_best_optimisation(double threshold, NNG *region_nng)
	{
		int num = 0;
		// std::cout << "cycle_best_optimisation is running" << std::endl;
		if (!region_nng->empty())
		{
			Edge *minimum_edge = region_nng->pop();
			if (minimum_edge->weight >= threshold)
			{
				return std::make_pair(true, num);
			}
			else
			{
				while (minimum_edge->weight < threshold)
				{
					this->simply_merge_vertices(minimum_edge);
					num++;
					if (region_nng->empty())
						break;
					minimum_edge = region_nng->pop();
				}
				return std::make_pair(false, num);
			}
		}
		return std::make_pair(true, num);
	}

	// kd-tree for partitioning the dataset
	std::vector<std::pair<T1 *, int>> partition_dataset_by_kd_tree(T1 *root, int partition_minmum_num)
	{
		std::vector<std::pair<T1 *, int>> results;
		if (root == nullptr)
			return results;
		// int position_dimension = root->position_dimension;
		int pos_dimension = root->position_dimension;
		int current_partition_size = 0; // = this->vertice_num;
		T1 *nil = root;
		while (nil != nullptr)
		{
			current_partition_size++;
			nil = (T1 *)nil->next_in_region;
		}
		Cluster *root_cluster = new Cluster((Vertex *)root, partition_minmum_num, pos_dimension, current_partition_size);

		std::cout << "[INFO] Calling(Into) func\"iter_divide_dataset\" (Recursive)" << std::endl;
		int size = this->iter_divide_dataset(root_cluster, partition_minmum_num, &results);
		// std::cout << "Total size: " << size << std::endl;
		delete root_cluster;
		return results;
	}

	// int iter_divide_dataset(Cluster *root_cluster, int partition_minmum_num, std::vector<std::pair<T1 *, int>> *divide_result)
	// {
	// 	if (root_cluster->dataset_size < partition_minmum_num)
	// 	{
	// 		divide_result->push_back(std::make_pair((T1 *)root_cluster->root_vertex, root_cluster->dataset_size));
	// 		return root_cluster->dataset_size;
	// 	}
	// 	// std::pair<Cluster*, Cluster*> child_clusters = root_cluster->divide_cluster_parallel(); //root_cluster->divide_cluster();
	// 	std::pair<Cluster *, Cluster *> child_clusters = root_cluster->divide_cluster();

	// 	int left_size = iter_divide_dataset(child_clusters.first, partition_minmum_num, divide_result);
	// 	int right_size = iter_divide_dataset(child_clusters.second, partition_minmum_num, divide_result);
	// 	return left_size + right_size;
	// }

	int iter_divide_dataset(Cluster *root_cluster, int partition_minmum_num, std::vector<std::pair<T1 *, int>> *divide_result)
	{
		std::cout << "[INFO] Partition num " << partition_minmum_num << std::endl;
		std::stack<Cluster *> cluster_stack;
		cluster_stack.push(root_cluster);
		int total_size = 0;
		size_t iteration = 0;

		while (!cluster_stack.empty())
		{
			++iteration;
			if (iteration % 10000 == 0)
			{
				std::cout << "[DEBUG] iteration=" << iteration
						  << "  stack_size=" << cluster_stack.size() << "\n";

				std::cout << "[DEBUG] iteration=" << iteration
						  << "  result region vector size=" << divide_result->size() << "\n";

				std::cout << "[DEBUG] iteration=" << iteration
						  << "  total size=" << total_size << "\n";
			}

			Cluster *current = cluster_stack.top();
			cluster_stack.pop();

			if (iteration % 10000 == 0)
			{
				std::cout << "[DEBUG] popping cluster @" << current
						  << "  size=" << current->dataset_size << "\n";
			}

			if (current->dataset_size < partition_minmum_num)
			{
				divide_result->push_back(std::make_pair((T1 *)current->root_vertex, current->dataset_size));
				total_size += current->dataset_size;
				// We release the memory at Other place!
			}
			else
			{
				std::pair<Cluster *, Cluster *> child_clusters = current->divide_cluster();
				// sanity-check: are they really new?
				if (child_clusters.first == current ||
					child_clusters.second == current)
				{
					std::cerr << "[WARNING] divide_cluster returned the same pointer!\n";
					break;
				}
				cluster_stack.push(child_clusters.second);
				cluster_stack.push(child_clusters.first);
			}
		}
		return total_size;
	}

	void check_direct(Vertex *root)
	{
		Vertex *v = root;
		while (v != nullptr)
		{
			if (v->direct == nullptr && v->edges.size() != 0)
				throw std::runtime_error("Error: none direct vertex existed!");
			;
			v = v->next_in_region;
		}
	}

	int check_data_size_in_region()
	{

		int num = 0;
		for (int i = 0; i < this->vertice_num; i++)
		{
			Vertex *v = &this->vertice[i];
			Vertex *nil = v;
			if (v->state != State::Freeze && v->is_visited == false)
			{
				// while (nil != nullptr)
				//{
				//	nil->is_visited = true;
				//	num += nil->size;
				//	nil = nil->next_in_region;
				// }
				num++;
			}
		}

		for (int i = 0; i < this->vertice_num; i++)
		{
			Vertex *v = &this->vertice[i];
			v->is_visited = false;
		}

		return num;
	}

	int check_data_size_in_region(Vertex *root)
	{

		int num = 0;
		while (root != nullptr)
		{
			Vertex *v = root;
			Vertex *nil = v;
			if (v->state != State::Freeze && v->is_visited == false)
			{
				num++;
			}
			root = root->next_in_region;
		}
		while (root != nullptr)
		{
			root->is_visited = false;
			root = root->next_in_region;
		}

		return num;
	}

	int check_data_size_in_output()
	{
		int num = 0;
		for (int i = 0; i < this->vertice_num; i++)
		{
			Vertex *v = &this->vertice[i];
			if (v->state != State::Freeze)
			{
				v->is_visited = true;
				num += v->size;
				Vertex *nt = v->next;
				int ttt = 0;
				Vertex *pre = nullptr;
				while (nt != nullptr)
				{

					if (nt->is_visited == true)
					{
						ttt++;
					}
					if (nt->ID == 193591)
					{
						pre;
						ttt++;
					}
					ttt++;
					nt->is_visited = true;
					pre = nt;
					nt = nt->next;
				}
			}
		}

		for (int i = 0; i < this->vertice_num; i++)
		{
			Vertex *v = &this->vertice[i];
			v->is_visited = false;
		}
		return num;
	}

	std::pair<T1 *, int> chaining_vertice()
	{
		T1 *root = nullptr;
		T1 *nil = root;
		int num = 0;
		for (int i = 0; i < this->vertice_num; i++)
		{
			T1 *v = &this->vertice[i];
			if (v->state != State::Freeze)
			{
				if (v->parent != nullptr)
				{
					std::cout << "wrong" << std::endl;
				}
				if (root == nullptr)
				{
					root = v;
					nil = root;
				}
				else
				{
					nil->next_in_region = v;
					nil = v;
				}
				nil->next_in_region = nullptr;
				v->bits_num = 0;
				v->bit_sign = 0;
				v->direct = nullptr;
				v->is_visited = false;
				v->state = State::Activate;
				num += v->size;
			}
		}
		return std::make_pair(root, num);
	}

	std::vector<Vertex *> partition_dataset_by_BFS(int partition_minmum_num, int sign)
	{
		std::vector<Vertex *> roots;
		if (this->vertice == nullptr)
			return roots;
		int total_num = 0;
		int bit_sign = sign;
		for (int i = 0; i < this->vertice_num; i++)
		{
			Vertex *v = &this->vertice[i];
			while (v->parent != nullptr)
			{
				v = v->parent;
			}
			if (bit_sign == v->bit_sign)
			{
				int num = this->broad_fist_search(v, partition_minmum_num);
				total_num += num;
				roots.push_back(v);
			}
		}
		std::cout << "Vertice num: " << total_num << std::endl;
		return roots;
	}

	int broad_fist_search(Vertex *root, int limit_num)
	{
		std::queue<Vertex *> queue;
		queue.push(root);
		root->is_visited = true;
		int search_num = 0;
		int bit_sign = root->bit_sign + 1;
		Vertex *ptr_next = root;
		root->next_in_region = nullptr;
		while (queue.size() > 0)
		{
			Vertex *v = queue.front();
			queue.pop();

			if (v != root)
			{
				ptr_next->next_in_region = v;
				ptr_next = v;
				v->next_in_region = nullptr;
			}
			v->append_sign(true);
			search_num++;
			if (search_num >= limit_num)
				break;
			for (int i = 0; i < v->edges.size(); i++)
			{
				int e_id = v->edges[i];
				Edge *e = &this->edges[e_id];
				Vertex *v_adj = e->get_other_vertex(v);
				if (v_adj->bit_sign != bit_sign && v_adj->is_visited == false)
				{
					// v_adj->append_sign(true);
					queue.push(v_adj);
					v_adj->is_visited = true;
				}
			}
		}

		return search_num;
	}

	void check_ouput_num()
	{
		int ouput_num = 0;
		if (this->vertice != nullptr)
		{
			// std::cout << "Ouput vertice: ";
			for (int i = 0; i < this->vertice_num; i++)
			{
				Vertex *v = &this->vertice[i];
				if (v->parent == nullptr)
				{
					ouput_num++;
					// std::cout << v->ID << " ";
				}
			}
		}
		std::cout << std::endl
				  << "Check: Num of Ouputs: " << GREEN_TEXT << ouput_num << RESET_COLOR << std::endl;
	}

	void suspend_edges(double threshold)
	{
		int suspend = 0;
		if (this->edges != nullptr)
		{
			for (int i = 0; i < this->edge_num; i++)
			{
				T2 *e_ptr = &this->edges[i];
				if (e_ptr->weight > threshold && e_ptr->state == State::Activate)
				{
					suspend++;
					this->suspend_edge(e_ptr);
				}
			}
		}
		std::cout << "Number of suspend edges: " << GREEN_TEXT << suspend << RESET_COLOR << std::endl;
	}

	void suspend_edge(Edge *e_ptr)
	{
		Vertex *left_vertex = e_ptr->left_vertex;
		left_vertex->remove_edge(e_ptr->ID);
		Vertex *right_vertex = e_ptr->right_vertex;
		right_vertex->remove_edge(e_ptr->ID);
		e_ptr->state = State::Suspend;
	}

	void freeze_edges(double threshold)
	{
		int freeze = 0;

		if (this->edges != nullptr)
		{
			for (int i = 0; i < this->edge_num; i++)
			{
				T2 *e_ptr = &this->edges[i];
				if (e_ptr->state == State::Activate && e_ptr->weight >= threshold)
				{
					freeze++;
					this->freeze_edge(e_ptr);
				}
			}
		}
		std::cout << "Number of freeze edges: " << GREEN_TEXT << freeze << RESET_COLOR << std::endl;
	}

	void freeze_edge(Edge *e_ptr)
	{
		Vertex *left_vertex = e_ptr->left_vertex;
		left_vertex->remove_edge(e_ptr->ID);
		Vertex *right_vertex = e_ptr->right_vertex;
		right_vertex->remove_edge(e_ptr->ID);
		e_ptr->state = State::Freeze;
	}

	void update_weight_in_region(Edge *root, State state)
	{
		while (root != nullptr)
		{
			this->update_edge_weight(root, state);
			root = root->next_edge_in_region;
		}
	}

	void update_weights(State state)
	{
		if (this->edges == nullptr)
			return;
		// multiThreadProcessing();
		for (int i = 0; i < this->edge_num; i++)
		{
			T2 *e_ptr = &this->edges[i];
			this->update_edge_weight(e_ptr, state);
		}
	}

	void update_edge_weight(Edge *e_ptr, State state)
	{
		e_ptr->is_checked = false;
		if (e_ptr->state == state)
		{
			e_ptr->weight = this->merging_criteria(
				(T_Feature *)e_ptr->left_vertex->ptr_features,
				(T_Feature *)e_ptr->right_vertex->ptr_features,
				(T2 *)e_ptr);
			e_ptr->state = State::Activate;
		}
	}

	void multi_thread_update_weights(State state)
	{
		// �߳����������Ը���ʵ���������
		int num_threads = std::thread::hardware_concurrency();
		// �����񻮷ֳɶ����
		int chunk_size = this->edge_num / num_threads;
		if (this->edge_num % num_threads != 0)
		{
			chunk_size++;
		}

		// �����̳߳�
		std::vector<std::future<void>> futures;

		// �����߳�
		for (int i = 0; i < num_threads; i++)
		{
			int start = i * chunk_size;
			int end = std::min(start + chunk_size, edge_num);
			futures.push_back(std::async(std::launch::async,
										 [this, start, end, state]()
										 {
											 for (int j = start; j < end; j++)
											 {
												 this->update_edge_weight(&this->edges[j], state);
											 }
										 }));
		}

		// �ȴ������߳����
		for (auto &future : futures)
		{
			future.get();
		}
	}

	void activate_edges(State state, string info = "")
	{
		// int activate = 0;
		// int suspend = 0;
		// int freeze = 0;
		if (this->edges == nullptr)
			return;
		for (int i = 0; i < this->edge_num; i++)
		{
			Edge *e_ptr = &this->edges[i];
			e_ptr->is_checked = false;
			if (e_ptr->state == state)
			{
				State final_state;
				if (info == "")
					final_state = this->activate_edge(e_ptr);
				else
					final_state = this->activate_edge_in_cycle(e_ptr);
			}
		}
		/*switch (final_state)
		{
		case State::Activate:
			activate++;
			break;
		case State::Suspend:
			suspend++;
			break;
		case State::Freeze:
			freeze++;
			break;
		default:
			break;
		}*/
		// std::cout << "Number of activate edges: " << GREEN_TEXT << activate << RESET_COLOR << std::endl;
		// std::cout << "Number of still suspend edges: " << GREEN_TEXT << suspend << RESET_COLOR << std::endl;
		// std::cout << "Number of new freeze edges: " << GREEN_TEXT << freeze << RESET_COLOR << std::endl;
	}

	State activate_edge(Edge *e_ptr)
	{
		// ԭ��edge���ߵ�vertexû�б仯�����
		if (e_ptr->left_vertex->state == State::Activate && e_ptr->right_vertex->state == State::Activate)
		{
			e_ptr->left_vertex->push(e_ptr);
			e_ptr->right_vertex->push(e_ptr);
			e_ptr->state = State::Activate;

			return e_ptr->state;
		}
		Vertex *left_vertex = e_ptr->left_vertex;
		Vertex *right_vertex = e_ptr->right_vertex;
		while (left_vertex->parent != nullptr)
		{
			left_vertex = left_vertex->parent;
		}
		while (right_vertex->parent != nullptr)
		{
			right_vertex = right_vertex->parent;
		}

		e_ptr->left_vertex = left_vertex;
		e_ptr->right_vertex = right_vertex;

		// ԭ��edge���ߵĽڵ�ͨ���������Ӻϲ�Ϊһ��������
		if (e_ptr->left_vertex == e_ptr->right_vertex)
		{
			e_ptr->state = State::Freeze;
			return e_ptr->state;
		}
		// �ж��Ƿ���������
		for (int i = 0; i < e_ptr->left_vertex->edges.size(); i++)
		{
			T2 *test_e = &this->edges[e_ptr->left_vertex->edges[i]];
			if (e_ptr->compareByEnds(*test_e)) // ������������
			{
				test_e->length += e_ptr->length;
				e_ptr->weight = e_ptr->weight <= test_e->weight ? e_ptr->weight : test_e->weight;
				e_ptr->state = State::Freeze;
				return e_ptr->state;
			}
		}

		// �������ȷ��e_ptr���������ڵ�Ψһ������
		e_ptr->left_vertex->push(e_ptr);
		e_ptr->right_vertex->push(e_ptr);
		e_ptr->weight = this->merging_criteria((T_Feature *)e_ptr->left_vertex->ptr_features, (T_Feature *)e_ptr->right_vertex->ptr_features, (T2 *)e_ptr);
		e_ptr->state = State::Activate;
		return e_ptr->state;
	}

	std::shared_mutex vertexMutex; // ���ڱ��� Vertex ����
	std::shared_mutex edgeMutex;   // ���ڱ��� Edge ����
	std::shared_mutex edgesMutex;  // ���ڱ��� this->edges ����
	State activate_edge_parallel(Edge *e_ptr)
	{

		if (e_ptr->left_vertex->state == State::Activate && e_ptr->right_vertex->state == State::Activate)
		{
			e_ptr->left_vertex->push(e_ptr);
			e_ptr->left_vertex->push(e_ptr);
			e_ptr->right_vertex->push(e_ptr);
			e_ptr->state = State::Activate;
			return e_ptr->state;
		}

		Vertex *left_vertex = e_ptr->left_vertex;
		Vertex *right_vertex = e_ptr->right_vertex;

		while (left_vertex->parent != nullptr)
		{
			left_vertex = left_vertex->parent;
		}
		while (right_vertex->parent != nullptr)
		{
			right_vertex = right_vertex->parent;
		}
		e_ptr->left_vertex = left_vertex;
		e_ptr->right_vertex = right_vertex;

		// ԭ��edge���ߵĽڵ�ͨ���������Ӻϲ�Ϊһ�������� - ��ȡ����
		if (e_ptr->left_vertex == e_ptr->right_vertex)
		{
			e_ptr->state = State::Freeze;
			return e_ptr->state;
		}

		// �ж��Ƿ��������� - ��ȡ this->edges �� e_ptr->left_vertex->edges

		for (int i = 0; i < e_ptr->left_vertex->edges.size(); i++)
		{
			T2 *test_e = &this->edges[e_ptr->left_vertex->edges[i]];
			if (e_ptr->compareByEnds(*test_e))
			{
				test_e->length += e_ptr->length;
				e_ptr->state = State::Freeze;
				return e_ptr->state;
			}
		}

		// �������ȷ��e_ptr���������ڵ�Ψһ������

		e_ptr->left_vertex->push(e_ptr);
		e_ptr->right_vertex->push(e_ptr);

		e_ptr->weight = this->merging_criteria((T_Feature *)e_ptr->left_vertex->ptr_features, (T_Feature *)e_ptr->right_vertex->ptr_features, (T2 *)e_ptr);
		e_ptr->state = State::Activate;

		return e_ptr->state;
	}

	Edge *activate_edges_in_region(Edge *e_root)
	{
		Edge *new_e_root = nullptr;
		Edge *nil = nullptr;
		while (e_root != nullptr)
		{
			e_root->is_checked = false;

			if (e_root->state != State::Freeze)
			{

				State final_state = this->activate_edge_in_cycle(e_root);
				if (new_e_root == nullptr && final_state != State::Freeze)
				{
					new_e_root = e_root;
					nil = new_e_root;
					e_root = e_root->next_edge_in_region;

					continue;
				}
				if (final_state != State::Freeze)
				{
					nil->next_edge_in_region = e_root;
					nil = e_root;
				}
			}
			e_root = e_root->next_edge_in_region;
		}
		if (nil != nullptr)
			nil->next_edge_in_region = nullptr;
		/*	std::cout << "Number of activate edges: " << GREEN_TEXT << activate << RESET_COLOR << std::endl;
			std::cout << "Number of still suspend edges: " << GREEN_TEXT << suspend << RESET_COLOR << std::endl;
			std::cout << "Number of new freeze edges: " << GREEN_TEXT << freeze << RESET_COLOR << std::endl;*/
		return new_e_root;
	}

	// work
	State activate_edge_in_cycle(Edge *e_ptr)
	{
		// ԭ��edge���ߵ�vertexû�б仯�����
		if (e_ptr->left_vertex->state == State::Activate && e_ptr->right_vertex->state == State::Activate)
		{
			// Vertex didnot remove any edges in cycle optimisation method.Thus, No need to push edge to vertex here.
			e_ptr->state = State::Activate;
			return e_ptr->state;
		}

		// The vertex ends of the edge changed
		Vertex *left_vertex = e_ptr->left_vertex;
		Vertex *right_vertex = e_ptr->right_vertex;

		int test = 0;
		while (left_vertex->parent != nullptr)
		{
			if (left_vertex->state != State::Freeze)
			{

				test++;
				throw std::runtime_error("State error!");
			}
			left_vertex = left_vertex->parent;
		}
		while (right_vertex->parent != nullptr)
		{
			if (right_vertex->state != State::Freeze)
			{
				test++;
				throw std::runtime_error("State error!");
			}
			right_vertex = right_vertex->parent;
		}

		e_ptr->left_vertex = left_vertex;
		e_ptr->right_vertex = right_vertex;

		// ԭ��edge���ߵĽڵ�ͨ���������Ӻϲ�Ϊһ��������
		if (e_ptr->left_vertex == e_ptr->right_vertex)
		{
			// Though vertex doesnot contation e_ptr, the code is safe. It means removing nothing.
			e_ptr->left_vertex->remove_edge(e_ptr->ID);
			e_ptr->right_vertex->remove_edge(e_ptr->ID);
			e_ptr->state = State::Freeze;
			return e_ptr->state;
		}

		// �ж��Ƿ��������� ���ܻ�������
		V_E_Topo_correct(e_ptr, e_ptr->left_vertex);
		V_E_Topo_correct(e_ptr, e_ptr->right_vertex);

		if (e_ptr->state == State::Freeze)
			return e_ptr->state;

		// �������ȷ��e_ptr���������ڵ�Ψһ������
		// e_ptr->weight = this->merging_criteria((T_Feature*)e_ptr->left_vertex->ptr_features, (T_Feature*)e_ptr->right_vertex->ptr_features, (T2*)e_ptr);
		// e_ptr->state = State::Activate;
		e_ptr->state = State::Suspend;
		return e_ptr->state;
	}

	void V_E_Topo_correct(Edge *e_ptr, Vertex *v_of_E_ends)
	{
		if (e_ptr->state == State::Freeze)
		{
			v_of_E_ends->remove_edge(e_ptr->ID);
			return;
		}
		int equal_num = 0;
		Edge *existed_other_edge = nullptr;
		bool is_existed_same_edge = false;
		for (int i = 0; i < v_of_E_ends->edges.size(); i++)
		{
			Edge *test_e = &this->edges[v_of_E_ends->edges[i]];
			if (e_ptr->compareByEnds(*test_e)) // ������������
			{
				if (e_ptr->compareByAddress(*test_e))
					is_existed_same_edge = true;
				else
					existed_other_edge = test_e;
				equal_num++;
			}
		}

		switch (equal_num)
		{
		case 0:
			v_of_E_ends->push(e_ptr);
			break;
		case 1:
			if (existed_other_edge != nullptr)
			{
				existed_other_edge->length += e_ptr->length;
				existed_other_edge->weight = existed_other_edge->weight <= e_ptr->weight ? existed_other_edge->weight : e_ptr->weight;
				e_ptr->state = State::Freeze;
				v_of_E_ends->remove_edge(e_ptr->ID); // In pratical, v_of_E_ends has no edge* e_ptr
			}
			break;
		case 2:
			if (is_existed_same_edge == true && existed_other_edge != nullptr)
			{
				existed_other_edge->length += e_ptr->length;
				existed_other_edge->weight = existed_other_edge->weight <= e_ptr->weight ? existed_other_edge->weight : e_ptr->weight;
				e_ptr->state = State::Freeze;
				v_of_E_ends->remove_edge(e_ptr->ID);
			}
			else
			{
				throw std::runtime_error("Wrong linked edges!");
			}
			break;
		default:
			throw std::runtime_error("Wrong linked edges!");
			break;
		}
	}

	void propagating_color(T1 *v)
	{
		int color = v->color;
		T1 *SIL = v;
		while (SIL != nullptr)
		{
			SIL->color = color;
			SIL = (T1 *)SIL->next;
			// SIL = (T1*)SIL->next_in_region;
		}
	}

	// int test = 0;
public:
	void initial_vertex_color() // WP coloring method
	{
		// this->check_ouput_num();
		// this->check_single_vertex();
		this->activate_edges(State::Freeze);
		this->activate_edges(State::Suspend);

		this->check_ouput_num();

		for (int i = 0; i < this->vertice_num; i++)
		{
			T1 *v = &this->vertice[i];
			if (v->state != State::Freeze)
			{
				// test_size += v->size;
				this->max_priqueue.push(v);
			}
		}
		// std::cout << "output size: " << test_size << std::endl;
		const int color_num = max_priqueue.top()->edges.size();
		int *color_arr = new int[color_num];

		while (!this->max_priqueue.empty())
		{
			// std::cout << this->max_priqueue.size() << std::endl;
			for (int i = 0; i < color_num; i++)
			{
				color_arr[i] = -1;
			}
			const T1 *const_ptr = max_priqueue.top();
			const T1 &const_v = *const_ptr;
			T1 *v = &this->vertice[const_v.ID];
			int color = -1;
			if (v->size == 0)
				color = 0;
			int edge_size = v->edges.size();
			for (int i = 0; i < edge_size; i++)
			{
				Edge *e = &this->edges[v->edges[i]];
				T1 *other_v = (T1 *)e->get_other_vertex(v);
				if (other_v->color > -1)
				{
					color_arr[other_v->color]++;
				}
			}
			for (int i = 0; i < color_num; i++)
			{
				if (color_arr[i] == -1)
				{
					color = i;
					break;
				}
			}

			v->color = color;
			this->propagating_color(v);
			this->max_priqueue.pop();
		}
		delete[] color_arr;
	}

	virtual void connect_vertex_and_edge() = 0;
	virtual void output(const char *) = 0;

public:
	int test = 0;
	UpdateFeatures update_features_func;
	MergingCriteria merging_criteria;
	T1 *vertice;
	T2 *edges;
	NNG *nng = nullptr;
	std::priority_queue<T1 *, std::vector<T1 *>, CompareAPointers> max_priqueue;
	int edge_num;
	int vertice_num;
};
