#pragma once
#include <map>
#include <queue>
#include "Vertex_Cloud.h"
#include "Edge_Cloud.h"
#include "Point_Cloud.h"
#include "BaseFeatures.h"
#include "NNG.h"
#include "RegionAdjacencyGraph.h"
#include "RAGCallbacks.h"

using PC_GetFeat = RAGCallbacks::GetFeatures<PointCloudFeatures, Edge_Cloud, Vertex_Cloud>;
using PC_Update = RAGCallbacks::UpdateFeatures<PointCloudFeatures, Edge_Cloud, Vertex_Cloud>;
using PC_Merge = RAGCallbacks::MergingCriteria<PointCloudFeatures, Edge_Cloud, Vertex_Cloud>;

class RAG_Cloud : public RegionAdjacencyGraph<Vertex_Cloud, Edge_Cloud, PointCloudFeatures>
{
public:
	PC_GetFeat get_features_func;
	RAG_Cloud(Point_Cloud *pc, double alpha, double search_radius,
			  // PointCloudNamespace::GetFeatures get_vertex_features,
			  PC_Update update_features,
			  PC_Merge MC);

	RAG_Cloud(Vertex_Cloud *root, int size, double alpha, double search_radius,
			  PC_Update update_features,
			  PC_Merge MC);

	RAG_Cloud(Vertex_Cloud *root, int size, double alpha,
			  PC_Update update_features,
			  PC_Merge MC);

	~RAG_Cloud();

	void output(const char *) override;
	std::pair<Vertex *, int> get_outputs();

private:
	PointCloudFeatures *features = nullptr;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud = nullptr;
	Point_Cloud *PC = nullptr;
	double alpha;
	double search_radius;

private:
	void connect_vertex_and_edge() override;
	void connect_vertex_and_edge2(Vertex_Cloud *root);
	void connect_vertex_and_edge3(Vertex_Cloud *root);
	// double merge_E_Sij(Vertex_Cloud*, Vertex_Cloud*);
	int kd_tree_radius_consturct_topo(double);
	int kd_tree_knn_consturct_topo(int);
	int full_linkage_topo();
	void initial_edges(int edge_index);
	void initial_vertice();
	void initial_vertice(Vertex_Cloud *);
	pcl::PointCloud<pcl::PointXYZ>::Ptr initial_cloud(Vertex_Cloud *);
	int test = 0;

	//
	// pcl::PointCloud<pcl::Normal>::Ptr normal_cloud = nullptr;
	//
};
