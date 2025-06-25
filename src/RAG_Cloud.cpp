#include "RAG_Cloud.h"
#include <pcl/point_cloud.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <vector>
#include <cmath>
#include <pcl/search/kdtree.h>
// #include <pcl/segmentation/mean_shift.h>

RAG_Cloud::RAG_Cloud(Point_Cloud *pc, double alpha, double search_radius,
					 // PointCloudNamespace::GetFeatures get_vertex_features,
					 PC_Update update_features,
					 PC_Merge MC) : RegionAdjacencyGraph(pc->count,
														 -1,
														 update_features,
														 MC) //,// get_features_func(get_vertex_features)
{
	auto start_time = std::chrono::high_resolution_clock::now();
	this->alpha = alpha;
	this->search_radius = search_radius;
	this->PC = pc;

	this->cloud = pc->pcl_cloud;
	this->vertice_num = this->cloud->size();
	this->features = new PointCloudFeatures[this->vertice_num];
	this->vertice = new Vertex_Cloud[this->vertice_num];
	this->connect_vertex_and_edge();

	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	if (duration.count() < 10)
		std::cout << "RAG is construted: " << GREEN_TEXT << duration.count() << RESET_COLOR << "s" << std::endl;
	else
		std::cout << "RAG is construted: " << RED_TEXT << duration.count() << RESET_COLOR << "s" << std::endl;
}

RAG_Cloud::RAG_Cloud(Vertex_Cloud *root, int size, double alpha, double search_radius,
					 PC_Update update_features,
					 PC_Merge MC) : RegionAdjacencyGraph(size,
														 -1,
														 update_features,
														 MC)
{
	auto start_time = std::chrono::high_resolution_clock::now();
	this->alpha = alpha;
	this->search_radius = search_radius;
	this->cloud = this->initial_cloud(root);
	this->vertice_num = size;
	this->features = new PointCloudFeatures[this->vertice_num];
	this->vertice = new Vertex_Cloud[this->vertice_num];
	this->connect_vertex_and_edge2(root);

	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
	// if (duration.count() < 10) std::cout << "RAG is construted: " << GREEN_TEXT << duration.count() << RESET_COLOR << "s" << std::endl;
	// else 					   std::cout << "RAG is construted: " << RED_TEXT << duration.count() << RESET_COLOR << "s" << std::endl;
}

RAG_Cloud::RAG_Cloud(Vertex_Cloud *root, int size, double alpha,
					 PC_Update update_features,
					 PC_Merge MC) : RegionAdjacencyGraph(size,
														 -1,
														 update_features,
														 MC)
{
	auto start_time = std::chrono::high_resolution_clock::now();
	this->alpha = alpha;
	// this->search_radius = search_radius;
	this->cloud = this->initial_cloud(root);
	this->vertice_num = size;
	this->features = new PointCloudFeatures[this->vertice_num];
	this->vertice = new Vertex_Cloud[this->vertice_num];
	this->connect_vertex_and_edge3(root);

	auto end_time = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
}

RAG_Cloud::~RAG_Cloud()
{
	if (this->vertice != nullptr)
		delete[] this->vertice;
	if (this->edges != nullptr)
		delete[] this->edges;
	if (this->features != nullptr)
		delete[] this->features;
}

void RAG_Cloud::output(const char *output_path)
{
	// save data as las file
	int output_num = 0; // for counting
	int total_num = 0;
	for (int i = 0; i < this->vertice_num; i++)
	{
		Vertex_Cloud *v = &this->vertice[i];
		if (v->state != State::Freeze)
		{
			Vertex_Cloud *next = v;
			while (next != nullptr)
			{
				next->color = output_num;
				next = (Vertex_Cloud *)next->next;
				total_num++;
			}
			output_num++;
		}
	}
	std::cout << "Output num: " << output_num << std::endl;
	std::cout << "Total num: " << total_num << std::endl;

	LASreader *reader = this->PC->get_las_reader();
	LASheader header = this->PC->get_header();
	// header.point_data_format = 5;
	// header.point_data_record_length += 4; // ���� 4 ���ֽڣ����ڴ洢 int ���͵Ķ�������
	const char *LAS_ATTRIBUTE_NEW_ATTRIBUTE_NAME = "SegLabel";
	LASattribute new_attribute(LAS_ATTRIBUTE_I32, LAS_ATTRIBUTE_NEW_ATTRIBUTE_NAME);
	std::cout << "data_format: " << header.point_data_format << std::endl;
	std::cout << "data length" << header.point_data_record_length << std::endl;
	// header.point_data_format = 5; // ����Ϊ 5 �����ж��������
	header.point_data_record_length += 21; // ���Ӷ������Եĳ��� (U16 ռ�� 2 ���ֽ�)
	header.add_attribute(new_attribute);

	// header.point_data_format

	LASwriteOpener laswriteopener;
	laswriteopener.set_file_name(output_path); // ��������ļ���
	int attribute_index = reader->header.get_attribute_index("treeID");
	std::cout << "attribute_index" << attribute_index << std::endl;
	LASwriter *writer = laswriteopener.open(&header);
	if (writer == 0)
	{
		fprintf(stderr, "ERROR: could not open LAS writer.\n");
	}
	int vertex_ID = 0;
	while (reader->read_point())
	{
		Vertex_Cloud *v = &this->vertice[vertex_ID];
		// I32 new_attribute_value = v->color;
		// if (v->color == 2)
		//{
		//	LASpoint* point = &(reader->point);
		//	point->set_attribute(attribute_index, new_attribute_value);
		//	//*((I32*)point->extra_bytes) = new_attribute_value;
		//	writer->write_point(point);
		// }
		// vertex_ID++;

		I32 new_attribute_value = v->color;
		LASpoint *point = &(reader->point);
		point->set_attribute(attribute_index, new_attribute_value);
		//*((I32*)point->extra_bytes) = new_attribute_value;
		writer->write_point(point);
		vertex_ID++;
	}
	reader->close();
	writer->close();
}

std::pair<Vertex *, int> RAG_Cloud::get_outputs()
{
	int none_color = 0;
	for (int i = 0; i < this->vertice_num; i++)
	{
		Vertex_Cloud *v = &this->vertice[i];
		if (v->state != State::Freeze)
		{
			none_color++;
		}
	}
	return std::pair<Vertex *, int>();
}

void RAG_Cloud::connect_vertex_and_edge()
{
	this->initial_vertice();
	this->edge_num = this->kd_tree_radius_consturct_topo(this->search_radius);
	this->initial_edges(this->edge_num);
}

void RAG_Cloud::connect_vertex_and_edge2(Vertex_Cloud *root)
{
	this->initial_vertice(root);
	this->edge_num = this->kd_tree_radius_consturct_topo(this->search_radius);
	this->initial_edges(this->edge_num);
}

void RAG_Cloud::connect_vertex_and_edge3(Vertex_Cloud *root)
{
	this->initial_vertice(root);
	this->edge_num = this->full_linkage_topo();
	this->initial_edges(this->edge_num);
}

int RAG_Cloud::kd_tree_radius_consturct_topo(double bandwidth)
{
	if (this->vertice == nullptr)
		return -1;
	double radius = bandwidth;
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	kdtree.setInputCloud(this->cloud);

	int edge_index = 0;
	for (size_t i = 0; i < this->vertice_num; i++) // this->vertice_num; i++)
	{
		std::vector<int> pointIdxRadiusSearch;
		std::vector<float> pointRadiusSquaredDistance;
		pcl::PointXYZ pt = this->cloud->points[i];
		Vertex_Cloud *v = &this->vertice[i];
		if (kdtree.radiusSearch(pt, radius, pointIdxRadiusSearch, pointRadiusSquaredDistance) > 0)
		{
			for (size_t j = 1; j < pointIdxRadiusSearch.size(); j++)
			{
				int neighbor_index = pointIdxRadiusSearch[j];
				Vertex_Cloud *right_v = &this->vertice[neighbor_index];
				if (neighbor_index > i)
				{
					v->edges.push_back(edge_index);
					right_v->edges.push_back(edge_index);
					edge_index++;
				}
			}
			// indivisual point
			if (pointIdxRadiusSearch.size() == 1)
			{
				std::vector<int> pointIdxNKNSearch(2);
				std::vector<float> pointNKNSquaredDistance(2);
				if (kdtree.nearestKSearch(pt, 2, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
				{
					int neighbor_index = pointIdxNKNSearch[1];
					Vertex_Cloud *right_v = &this->vertice[neighbor_index];
					v->edges.push_back(edge_index);
					right_v->edges.push_back(edge_index);
					edge_index++;
				}
			}
		}
	}
	// this->edge_num = edge_index;
	// std::cout << "edge num: " << edge_index <<std::endl;
	return edge_index;
}

int RAG_Cloud::kd_tree_knn_consturct_topo(int K)
{
	if (this->vertice == nullptr)
		return -1;
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
	kdtree.setInputCloud(this->cloud);
	int edge_index = 0;
	std::cout << "knn..." << std::endl;
	for (size_t i = 0; i < this->vertice_num; i++) // this->vertice_num; i++)
	{
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);
		pcl::PointXYZ pt = this->cloud->points[i];

		Vertex_Cloud *v = &this->vertice[i];
		if (kdtree.nearestKSearch(pt, K, pointIdxNKNSearch, pointNKNSquaredDistance) > 0)
		{
			// skip the 1st point which is the pt
			for (size_t j = 1; j < pointIdxNKNSearch.size(); ++j)
			{
				int neighbor_index = pointIdxNKNSearch[j];
				Vertex_Cloud *right_v = &this->vertice[neighbor_index];
				this->edges[edge_index].ID = edge_index;
				this->edges[edge_index].left_vertex = v;
				this->edges[edge_index].right_vertex = right_v;
				// this->edges[edge_index].E_Sij = 1.0 / pow(2.0, ((PointCloudFeatures*)v->ptr_features)->alpha) *
				//	(pow(v->positions[0] - right_v->positions[0], 2.0) + pow(v->positions[1] - right_v->positions[1], 2.0));
				// double weight = this->merging_criteria(
				//	(PointCloudFeatures*)this->edges[edge_index].left_vertex->ptr_features,
				//	(PointCloudFeatures*)this->edges[edge_index].left_vertex->ptr_features,
				//	&this->edges[edge_index]);
				// this->edges[edge_index].weight = weight;
				v->edges.push_back(edge_index);
				// right_v->edges.push_back(edge_index);
				edge_index++;
			}
		}
	}
	std::cout << "check topo..." << std::endl;

	for (size_t i = 0; i < this->edge_num; i++)
	{
		Edge_Cloud *e = &this->edges[i];
		Vertex_Cloud *left_v = (Vertex_Cloud *)e->left_vertex;
		Vertex_Cloud *right_v = (Vertex_Cloud *)e->right_vertex;
		if (left_v == nullptr || right_v == nullptr)
		{
			std::cout << "Error ID: " << e->ID << " i: " << i;
			test++;
		}
		bool is_left_contain_edge = std::find(left_v->edges.begin(), left_v->edges.end(), e->ID) != left_v->edges.end();
		bool is_right_contain_edge = std::find(right_v->edges.begin(), right_v->edges.end(), e->ID) != right_v->edges.end();
		if (!is_left_contain_edge)
		{
			if (right_v->edges.size() > 1)
			{
				right_v->edges.erase(std::remove(right_v->edges.begin(), right_v->edges.end(), e->ID), right_v->edges.end());
				e->state = State::Freeze;
			}
			else
			{
				left_v->edges.push_back(e->ID);
			}
		}
		if (!is_right_contain_edge)
		{
			if (left_v->edges.size() > 1)
			{
				left_v->edges.erase(std::remove(left_v->edges.begin(), left_v->edges.end(), e->ID), left_v->edges.end());
				e->state = State::Freeze;
			}
			else
			{
				right_v->edges.push_back(e->ID);
			}
		}
	}
	return 0;
}

int RAG_Cloud::full_linkage_topo()
{
	if (this->vertice == nullptr)
		return -1;
	int edge_index = 0;
	for (size_t i = 0; i < this->vertice_num - 1; i++)
	{
		Vertex_Cloud *v_i = &this->vertice[i];
		for (size_t j = i + 1; j < this->vertice_num; j++)
		{
			Vertex_Cloud *v_j = &this->vertice[j];
			v_j->edges.push_back(edge_index);
			v_i->edges.push_back(edge_index);
			edge_index++;
		}
	}
	return edge_index;
}

void RAG_Cloud::initial_edges(int edge_num)
{
	// if (tag != 0)
	//{
	//	test++;
	// }
	this->edges = new Edge_Cloud[edge_num];
	for (size_t i = 0; i < this->vertice_num; i++)
	{
		Vertex_Cloud *v = &this->vertice[i];
		for (size_t j = 0; j < v->edges.size(); ++j)
		{
			int edge_id = v->edges[j];
			Edge_Cloud *e = &this->edges[edge_id];
			e->ID = edge_id;
			if (e->left_vertex == nullptr)
			{
				e->left_vertex = v;
			}
			else if (e->right_vertex == nullptr)
			{
				e->right_vertex = v;
			}
			else
			{
				throw std::runtime_error("Repeated edge error!");
			}
		}
	}

	for (size_t i = 0; i < this->edge_num; i++)
	{

		Edge_Cloud *e = &this->edges[i];
		if (e->left_vertex != nullptr || e->right_vertex != nullptr)
		{
			Vertex_Cloud *left_vertex = (Vertex_Cloud *)e->left_vertex;
			Vertex_Cloud *right_vertex = (Vertex_Cloud *)e->right_vertex;

			PointCloudFeatures *l_feature = (PointCloudFeatures *)(e->left_vertex->ptr_features);
			PointCloudFeatures *r_feature = (PointCloudFeatures *)(e->right_vertex->ptr_features);

			double weight = this->merging_criteria(l_feature, r_feature, e);
			e->weight = weight;
		}
		else
		{
			throw std::runtime_error("Edge error!");
		}
	}
}

void RAG_Cloud::initial_vertice()
{
	for (size_t i = 0; i < this->vertice_num; i++)
	{
		pcl::PointXYZ pt = this->cloud->points[i];
		Vertex_Cloud *v = &this->vertice[i];
		v->ID = i;
		this->features[v->ID].set_feature(v, 1, 0.0, this->alpha, pt.x, pt.y, pt.z);
		v->ptr_features = (RootFeatures *)&this->features[v->ID];
		v->set_position_dimention(3);
		v->positions[0] = pt.x;
		v->positions[1] = pt.y;
		v->positions[2] = pt.z;

		if (i + 1 < this->vertice_num)
		{
			this->vertice[i].next_in_region = &this->vertice[i + 1];
		} // chain
	}
}

void RAG_Cloud::initial_vertice(Vertex_Cloud *root)
{
	Vertex_Cloud *nil = root;
	int i = 0;
	while (nil != nullptr)
	{
		pcl::PointXYZ pt = this->cloud->points[i];
		Vertex_Cloud *v = &this->vertice[i];
		v->ID = i;
		v->Org_ID = nil->Org_ID;
		PointCloudFeatures *pt_features = (PointCloudFeatures *)nil->ptr_features;
		v->size = pt_features->N;
		this->features[v->ID].set_feature(v, pt_features->N, pt_features->E_Si, pt_features->alpha,
										  pt_features->x_bar, pt_features->y_bar, pt_features->z_bar);
		v->ptr_features = (RootFeatures *)&this->features[v->ID];
		v->set_position_dimention(3);
		v->positions[0] = pt.x;
		v->positions[1] = pt.y;
		v->positions[2] = pt.z;

		// ! setting the normal
		// !
		// Eigen::Vector3d sum_n(0, 0, 0);
		// Vertex_Cloud *t = nil;
		// while (t)
		// {
		// 	const auto &pn = normal_cloud->points[t->ID];
		// 	sum_n += Eigen::Vector3d(pn.normal_x, pn.normal_y, pn.normal_z);
		// 	t = (Vertex_Cloud *)t->next_in_region;
		// }
		// sum_n.normalize();
		// features[v->ID].normal = sum_n;
		// v->normal = sum_n;
		// !
		// ! end of setting the normal

		if (i + 1 < this->vertice_num)
		{
			this->vertice[i].next_in_region = &this->vertice[i + 1];
		}
		nil = (Vertex_Cloud *)nil->next_in_region;
		i++;
	}
}

pcl::PointCloud<pcl::PointXYZ>::Ptr RAG_Cloud::initial_cloud(Vertex_Cloud *v)
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr C(new pcl::PointCloud<pcl::PointXYZ>);
	Vertex_Cloud *nil = v;
	while (nil != nullptr)
	{
		pcl::PointXYZ point;
		point.x = nil->positions[0];
		point.y = nil->positions[1];
		point.z = nil->positions[2];
		C->points.push_back(point);
		nil = (Vertex_Cloud *)nil->next_in_region;
	}
	return C;
}
