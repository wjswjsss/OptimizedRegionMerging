#pragma once
#include <iostream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <Eigen/Core>
// #include "RabbaniVertex.h"
using namespace std;
using namespace boost::multiprecision;

class Vertex_Cloud;
class RabbaniVertex;
class RootFeatures
{
public:
	RootFeatures();
	virtual void FeatureInfo() = 0;
	// ������麯��
	virtual ~RootFeatures() {} // ����������
public:
	std::string name;
};

class ImageFeatures : RootFeatures
{
public:
	double area;
	double mean;
	double std;
	double mean_arry[3];
	double std_arry[3];
	int band_num;

public:
	ImageFeatures();
	ImageFeatures(int band_num);
	ImageFeatures(int *, int);
	virtual ~ImageFeatures() {};
	double Sum(int *, int);
	double Mean(int *, int);
};

class HSWOFeatures : public ImageFeatures
{
public:
	HSWOFeatures(int *, int);
	HSWOFeatures(double, double, int);
	~HSWOFeatures();
	void FeatureInfo() override;
};

class MRSFeatures : public ImageFeatures
{
public:
	MRSFeatures();
	MRSFeatures(int *, double *, int, int, int);
	MRSFeatures(
		double area,
		double mean,
		double std,
		double length,
		double color,
		int band_num,
		double *band_weights,
		int top, int bottom, int left, int right, int box,
		double smooth, double compt);
	~MRSFeatures();
	void FeatureInfo() override;
	void set_feature(int *DNs, double *args, int args_size, int x, int y);

public:
	double shape_weight, compt_weight;
	double band_weights[3];
	double color;
	double length;
	int top, bottom, left, right;
	int box;
	double smooth, compt;
};

class PointCloudFeatures : RootFeatures
{
public:
	PointCloudFeatures();
	~PointCloudFeatures();
	void set_feature(Vertex_Cloud *v, int N, double E_Si, double alpha, double x_bar, double y_bar, double z_bar);
	Vertex_Cloud *highest_v;
	int N;
	double E_Si;
	double alpha;
	double x_bar;
	double y_bar;
	double z_bar;
	void FeatureInfo() override;
};

class BloomingTreeFeatures : RootFeatures
{
public:
	BloomingTreeFeatures();
	~BloomingTreeFeatures();
	void set_feature(int mass,
					 cpp_dec_float_50 G,
					 cpp_dec_float_50 p,
					 cpp_dec_float_50 parallax,
					 cpp_dec_float_50 ra,
					 cpp_dec_float_50 dec,
					 cpp_dec_float_50 pmra,
					 cpp_dec_float_50 pmdec,
					 cpp_dec_float_50 phot_g_mean_mag);

public:
	int mass;
	cpp_dec_float_50 G; // gravitational constant
	cpp_dec_float_50 p; // portion
	cpp_dec_float_50 parallax;
	cpp_dec_float_50 ra;			  // right ascension
	cpp_dec_float_50 dec;			  // declination
	cpp_dec_float_50 pmra;			  // proper motion in right ascension
	cpp_dec_float_50 pmdec;			  // proper motion in declination
	cpp_dec_float_50 phot_g_mean_mag; // light degree
	cpp_dec_float_50 r;

	void FeatureInfo() override;
};

// + ----------------------------- +
// | Creating the Rabbani features |
// + ----------------------------- +
class RabbaniFeatures : public RootFeatures
{
public:
	RabbaniFeatures() : highest_v(nullptr),
						x_bar(0),
						y_bar(0),
						z_bar(0),
						size(0),
						normal(Eigen::Vector3d::Zero())
	{
	}

	~RabbaniFeatures() {}

	// call this when you first create a single-point region
	void set_feature(RabbaniVertex *v,
					 double x_bar_,
					 double y_bar_,
					 double z_bar_,
					 const Eigen::Vector3d &normal_)
	{
		highest_v = v;
		x_bar = x_bar_;
		y_bar = y_bar_;
		z_bar = z_bar_;
		normal = normal_;
		size = 1; // start with ONE point per Cluster.
	}

	Eigen::Vector3d normal;

	RabbaniVertex *highest_v;

	int size;
	// int N;						// number of points in region
	// double E_Si;				// internal energy
	// double alpha;				// user‐supplied
	double x_bar, y_bar, z_bar; // centroid

	void FeatureInfo() override
	{
		std::cout << "Rabbani Feature:\n"
				  << "  centroid=(" << x_bar << "," << y_bar << "," << z_bar << ")\n"
				  << "  normal=(" << normal.transpose() << ")\n";
	}
};
