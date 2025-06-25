#include "HRM3D.h"
#include "BaseFeatures.h"

// PointCloudFeatures* get_cloud_features(int*, double*, int, int, int);
void update_cloud_features(PointCloudFeatures *, PointCloudFeatures *, Edge_Cloud *);
double calculate_total_E_Sij(PointCloudFeatures *, PointCloudFeatures *);
double hrm3d_cloud_mc(PointCloudFeatures *, PointCloudFeatures *, Edge_Cloud *);

HRM3D::HRM3D(const char *cloud_file_char)
{
	// this->get_features_func = get_cloud_features;
	this->update_features_func = update_cloud_features;
	this->merging_criteria = hrm3d_cloud_mc;
	this->cloud = new Point_Cloud(cloud_file_char);
}

HRM3D::~HRM3D()
{
	delete this->cloud;
}

void HRM3D::run(double scale_parameter, double alpha, double search_radius, int region_size, RunMode mode)
{
	this->rag = new RAG_Cloud(this->cloud, alpha, search_radius, this->update_features_func, this->merging_criteria);
	std::cout << "Running..." << std::endl;
	int iter = 4;

	RAG_Cloud *temp_rag = this->rag;
	Vertex_Cloud *pre_root = nullptr;
	Vertex_Cloud *root = nullptr;
	double *radius_scale_factor = new double[3]{20.0, 40, 80.0};

	long long total_run_time = 0;
	for (int i = 0; i < iter; i++)
	{
		std::cout << "Round " << GREEN_TEXT << i + 1 << RESET_COLOR << " begin: " << std::endl;
		auto start_time = std::chrono::high_resolution_clock::now();
		root = temp_rag->run(scale_parameter, mode, region_size);
		auto end_time = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
		long long t = duration.count();
		total_run_time += t;
		std::cout << "Round " << GREEN_TEXT << i + 1 << RESET_COLOR << " time: " << GREEN_TEXT << t << RESET_COLOR << "s" << std::endl;

		if (pre_root == nullptr)
		{
			this->assignment_org_ID(root);
		}
		else
		{
			this->assignment_value_to_org_vertex(root);
		}
		pre_root = root;
		int size = this->check_num(root);
		std::cout << "HRM3D Round" << i + 1 << ": " << size << std::endl;
		if (i != iter - 1)
		{
			temp_rag = new RAG_Cloud(root, size, alpha, search_radius * radius_scale_factor[i], this->update_features_func, this->merging_criteria);
		}
		else
		{
			// break;
			temp_rag = new RAG_Cloud(root, size, alpha, this->update_features_func, this->merging_criteria);

			auto start_time = std::chrono::high_resolution_clock::now();
			root = temp_rag->run(7000.0, mode, region_size);
			auto end_time = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
			long long t = duration.count();
			total_run_time += t;

			int new_size4 = this->check_num(root);
			std::cout << "HRM3D Round5: " << new_size4 << std::endl;
			this->assignment_value_to_org_vertex(root);
		}
	}

	std::cout << "Run time: " << GREEN_TEXT << total_run_time << RESET_COLOR << "ms" << std::endl;
}

void HRM3D::assignment_org_ID(Vertex_Cloud *root)
{
	Vertex_Cloud *nil = root;
	while (nil != nullptr)
	{
		nil->Org_ID = nil->ID;
		nil = (Vertex_Cloud *)nil->next_in_region;
	}
}

void HRM3D::assignment_value_to_org_vertex(Vertex_Cloud *root)
{
	Vertex_Cloud *nil = root;
	while (nil != nullptr)
	{
		Vertex_Cloud *search_next = (Vertex_Cloud *)nil->next;
		Vertex_Cloud *org_vertex = &this->rag->vertice[nil->Org_ID];
		while (search_next != nullptr)
		{
			Vertex_Cloud *org_next = &this->rag->vertice[search_next->Org_ID];
			org_vertex->add_next_vertex(org_next);
			org_next->state = State::Freeze;
			search_next = (Vertex_Cloud *)search_next->next;
		}
		nil = (Vertex_Cloud *)nil->next_in_region;
	}
}

int HRM3D::check_num(Vertex_Cloud *root)
{
	Vertex_Cloud *nil = root;
	int new_size = 0;
	while (nil != nullptr)
	{
		new_size++;
		nil = (Vertex_Cloud *)nil->next_in_region;
	}
	return new_size;
}

void HRM3D::output(const char *las_file_char)
{
	std::cout << "output..." << std::endl;
	rag->output(las_file_char);
}

void update_cloud_features(PointCloudFeatures *b1, PointCloudFeatures *b2, Edge_Cloud *e)
{
	PointCloudFeatures *left = static_cast<PointCloudFeatures *>(b1);
	PointCloudFeatures *right = static_cast<PointCloudFeatures *>(b2);
	Vertex_Cloud *highest_vertex = (left->highest_v->positions[2] > right->highest_v->positions[2]) ? left->highest_v : right->highest_v;
	int total_N = left->N + right->N;
	left->x_bar = (left->x_bar * left->N + right->x_bar * right->N) / total_N;
	left->y_bar = (left->y_bar * left->N + right->y_bar * right->N) / total_N;
	left->z_bar = (left->z_bar * left->N + right->z_bar * right->N) / total_N;
	left->N = total_N;
	left->highest_v = highest_vertex;
	left->E_Si = e->E_Sij;
}

double hrm3d_cloud_mc(PointCloudFeatures *left_feature, PointCloudFeatures *right_feature, Edge_Cloud *e)
{
	double alpha = left_feature->alpha;
	double E_Sij = calculate_total_E_Sij(left_feature, right_feature);
	e->E_Sij = E_Sij; // will be used in the merger
	double G_Si = 1.0 / pow(1.0 * left_feature->N, alpha) * left_feature->E_Si;
	double G_Sj = 1.0 / pow(1.0 * right_feature->N, alpha) * right_feature->E_Si;
	double G_Sij = 1.0 / pow((1.0 * left_feature->N + 1.0 * right_feature->N), alpha) * e->E_Sij;

	double mc = max(0.0, G_Sij - G_Si - G_Sj);
	return mc;
}

double calculate_total_E_Sij(PointCloudFeatures *b1, PointCloudFeatures *b2)
{
	PointCloudFeatures *left = static_cast<PointCloudFeatures *>(b1);
	PointCloudFeatures *right = static_cast<PointCloudFeatures *>(b2);
	Vertex_Cloud *highest_vertex = (left->highest_v->positions[2] > right->highest_v->positions[2]) ? left->highest_v : right->highest_v;

	Vertex_Cloud *left_high = left->highest_v; // RAG merge function make sure left_high doesnot contaion right_high
	Vertex_Cloud *right_high = right->highest_v;

	double total_E_Si = 0.0;
	if (highest_vertex == left_high)
	{
		Vertex_Cloud *a = right_high;
		Vertex_Cloud *b = highest_vertex;
		PointCloudFeatures *a_feature = right;
		double delta_x = b->positions[0] - a->positions[0];
		double delta_y = b->positions[1] - a->positions[1];
		// double delta_z = b->positions[2] - a->positions[2];
		total_E_Si = left->E_Si + right->E_Si -
					 2.0 * a_feature->N * (delta_x * (a_feature->x_bar - a->positions[0]) + delta_y * (a_feature->y_bar - a->positions[1])) +
					 a_feature->N * (pow(delta_x, 2.0) + pow(delta_y, 2.0));
	}
	if (highest_vertex == right_high)
	{
		Vertex_Cloud *a = left_high;
		Vertex_Cloud *b = highest_vertex;
		PointCloudFeatures *a_feature = left;
		double delta_x = b->positions[0] - a->positions[0];
		double delta_y = b->positions[1] - a->positions[1];
		// double delta_z = b->positions[2] - a->positions[2];
		total_E_Si = left->E_Si + right->E_Si -
					 2.0 * a_feature->N * (delta_x * (a_feature->x_bar - a->positions[0]) + delta_y * (a_feature->y_bar - a->positions[1])) +
					 a_feature->N * (pow(delta_x, 2.0) + pow(delta_y, 2.0));
	}
	return total_E_Si;
}

// double normal_angle_difference(
// 	PointCloudFeatures *A,
// 	PointCloudFeatures *B,
// 	Edge_Cloud *e)
// {
// 	// 1) dot → clamp → acos → radians
// 	double d = A->normal.dot(B->normal);
// 	d = std::clamp(d, -1.0, 1.0);
// 	double angle = std::acos(d);

// 	// 2) store if you like
// 	e->weight = angle;

// 	// 3) return it
// 	return angle;
// }

// void merge_region_normals(
// 	PointCloudFeatures *A,
// 	PointCloudFeatures *B,
// 	Edge_Cloud *e)
// {
// 	// — existing fusion (counts, centroids, E_Si) —
// 	int totalN = A->N + B->N;
// 	A->x_bar = (A->x_bar * A->N + B->x_bar * B->N) / totalN;
// 	A->y_bar = (A->y_bar * A->N + B->y_bar * B->N) / totalN;
// 	A->z_bar = (A->z_bar * A->N + B->z_bar * B->N) / totalN;
// 	A->N = totalN;
// 	A->E_Si = e->E_Sij;

// 	// — new normal fusion —
// 	Eigen::Vector3d n1 = A->normal;
// 	Eigen::Vector3d n2 = B->normal;
// 	Eigen::Vector3d nm = (n1 * double(A->N) + n2 * double(B->N)) / double(totalN);
// 	if (nm.norm() > 1e-6)
// 		nm.normalize();
// 	else
// 		nm = n1; // fallback

// 	A->normal = nm;

// 	// — also write into the surviving Vertex_Cloud —
// 	// RegionAdjacencyGraph typically makes A’s vertex the “merged” one.
// 	// You can fetch it from the edge:
// 	Vertex_Cloud *surv = static_cast<Vertex_Cloud *>(e->left_vertex);
// 	surv->normal = nm;
// }