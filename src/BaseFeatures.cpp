#include<iostream>
#include <cstdlib>
#include <cmath>
#include "BaseFeatures.h"

#include <boost/multiprecision/cpp_dec_float.hpp>
using namespace std;
using namespace boost::multiprecision;


RootFeatures::RootFeatures()
{
	this->name = "features";
}

ImageFeatures::ImageFeatures():RootFeatures()
{
	//this->band_num = 3;
	//this->mean_arry = new double[this->band_num];
	//this->std_arry = new double[this->band_num];
}

ImageFeatures::ImageFeatures(int band_num):RootFeatures()
{
	this->band_num = band_num;
	//this->mean_arry = new double[this->band_num];
	//this->std_arry = new double[this->band_num];
}
ImageFeatures::ImageFeatures(int* dn_values,int band_num) :RootFeatures()
{
	this->area = 1.0;
	this->mean = this->Mean(dn_values, band_num);   
	this->std = 0.0;
	//this->mean_arry = new double[band_num];
	//this->std_arry = new double[band_num];
	for (int i = 0; i < band_num; i++)
	{
		this->mean_arry[i] = dn_values[i];
		this->std_arry[i] = 0.0;
	}
	this->band_num = band_num;
}

//ImageFeatures::~ImageFeatures()
//{
//	if(this->mean_arry != nullptr) delete[] this->mean_arry;
//	if(this->std_arry != nullptr) delete[] this->std_arry;
//}
double ImageFeatures::Sum(int* dn_values, int size)
{
	double total = 0.0;
	for (int i = 0; i < size; i++)
	{
		total += dn_values[i];
	}
	return total;
}
double ImageFeatures::Mean(int* pixels, int size)
{
	double mean = this->Sum(pixels, size) / size;
	return mean;
}

HSWOFeatures::HSWOFeatures(int* pixels, int size):ImageFeatures(pixels,size)
{
	this->area = 1.0;
	this->mean = this->Mean(pixels, size);
}

HSWOFeatures::HSWOFeatures(double area, double mean, int band_num):ImageFeatures(band_num)
{
	this->area = area;
	this->mean = mean;
}

void HSWOFeatures::FeatureInfo()
{
	std::cout << "do nothing" << endl;
}

HSWOFeatures::~HSWOFeatures()
{
	
}

//double Merge_stds(double a1, double a2, double std1, double std2, double u1, double u2);
MRSFeatures::MRSFeatures():ImageFeatures()
{
	double total = 0.0;
	//this->band_weights = new double[3];
	this->color = total;
	this->length = 4.0;
	this->top = 0.0;
	this->bottom = 0.0;
	this->left = 0.0;
	this->right = 0.0;
	this->box = 0.0;
	this->smooth = 0.0;
	this->compt = 0.0;
}

MRSFeatures::MRSFeatures(int* DNs,  double* args,int args_size,int x, int y) :ImageFeatures(DNs, args_size - 2)
{
	double total = 0.0;
	//this->band_weights = new double[this->band_num];
	for (int i = 0; i < band_num; i++)
	{
		this->band_weights[i] = args[i];
		total += this->band_weights[i] * this->area * this->std_arry[i];
	}
	this->shape_weight = args[band_num];
	this->compt_weight = args[band_num+1];
	this->color = total;
	this->length = 4.0;
	this->top = y;
	this->bottom = y;
	this->left = x;
	this->right = x;
	this->box = (abs(this->top - this->bottom) + 1) * (abs(this->left - this->right) + 1);
	this->smooth = this->area * this->length / this->box;
	this->compt  = this->area * this->length / sqrt(this->area);
}

MRSFeatures::MRSFeatures(
	double area, 
	double mean,
	double std, 
	double length,
	double color,
	int band_num,
	double* band_weights,
	int top, int bottom, int left, int right, int box,  
	double smooth, double compt):ImageFeatures(band_num)
{
	this->area = area;
	this->mean = mean;
	this->std = std;
	this->length = length;
	//this->band_weights = new double[band_num];
	for (int i = 0; i < band_num;i++)
	{
		this->band_weights[i] = band_weights[i];
	}
	this->color = color;
	this->top = top;
	this->bottom = bottom;
	this->left = left;
	this->right = right;
	this->box = box;
	this->smooth = smooth;
	this->compt = compt;
}

void MRSFeatures::set_feature(int* DNs, double* args, int args_size, int x, int y)
{
	this->band_num = args_size - 2;
	this->area = 1.0;
	this->mean = this->Mean(DNs, this->band_num);
	this->std = 0.0;
	//this->mean_arry = new double[this->band_num];
	//this->std_arry = new double[this->band_num];
	//this->band_weights = new double[this->band_num];

	for (int i = 0; i < this->band_num; i++)
	{
		this->mean_arry[i] = DNs[i];
		this->std_arry[i] = 0.0;
	}

	double total = 0.0;
	for (int i = 0; i < band_num; i++)
	{
		this->band_weights[i] = args[i];
		total += this->band_weights[i] * this->area * this->std_arry[i];
	}
	this->shape_weight = args[band_num];
	this->compt_weight = args[band_num + 1];
	this->color = total;
	this->length = 4.0;
	this->top = y;
	this->bottom = y;
	this->left = x;
	this->right = x;
	this->box = (abs(this->top - this->bottom) + 1) * (abs(this->left - this->right) + 1);
	this->smooth = this->area * this->length / this->box;
	this->compt = this->area * this->length / sqrt(this->area);
}

MRSFeatures::~MRSFeatures()
{
	//delete[] this->band_weights;
	//delete[] this->mean_arry;
	//delete[] this->std_arry;
}

void MRSFeatures::FeatureInfo()
{
	std::cout << "do nothing" << endl;
}


PointCloudFeatures::PointCloudFeatures() :RootFeatures()
{
	this->highest_v = nullptr;
	//this->dis_p = 0.0;
	this->N = 1;
	this->E_Si = 0.0;
	this->alpha = 1.0;
	this->x_bar = 0.0;
	this->y_bar = 0.0;
	this->z_bar = 0.0;
}

PointCloudFeatures::~PointCloudFeatures()
{
}

void PointCloudFeatures::set_feature(Vertex_Cloud* v, int N, double E_Si, double alpha, double x_bar, double y_bar, double z_bar)
{
	this->highest_v = v;
	this->N = N;
	this->E_Si = E_Si;
	this->alpha = alpha;
	this->x_bar = x_bar;
	this->y_bar = y_bar;
	this->z_bar = z_bar;
}

void PointCloudFeatures::FeatureInfo()
{
	std::cout << "PointCloudFeatures" << std::endl;
}

BloomingTreeFeatures::BloomingTreeFeatures() :RootFeatures()
{
	this->mass = 0;
	this->G = cpp_dec_float_50("0");
	this->p = cpp_dec_float_50("0");
	this->parallax = cpp_dec_float_50("0");
	this->ra = cpp_dec_float_50("0");
	this->dec = cpp_dec_float_50("0");
	this->pmra = cpp_dec_float_50("0");
	this->pmdec = cpp_dec_float_50("0");
	this->phot_g_mean_mag = cpp_dec_float_50("0");
	this->r = cpp_dec_float_50("0");
}

void BloomingTreeFeatures::set_feature(int mass,
	cpp_dec_float_50 G, cpp_dec_float_50 p,
	cpp_dec_float_50 parallax, cpp_dec_float_50 ra, 
	cpp_dec_float_50 dec, cpp_dec_float_50 pmra, cpp_dec_float_50 pmdec, cpp_dec_float_50 phot_g_mean_mag)
{
	this->mass = mass;
	this->G = cpp_dec_float_50(G);
	this->p = cpp_dec_float_50(p);
	this->parallax = cpp_dec_float_50(parallax);
	this->ra = cpp_dec_float_50(ra);
	this->dec = cpp_dec_float_50(dec);
	this->pmra = cpp_dec_float_50(pmra);
	this->pmdec = cpp_dec_float_50(pmdec);
	this->phot_g_mean_mag = cpp_dec_float_50(phot_g_mean_mag);
	this->r = 1.0 / this->parallax;
}


void BloomingTreeFeatures::FeatureInfo()
{
	std::cout << "BloomingTreeFeatures" << std::endl;
}

BloomingTreeFeatures::~BloomingTreeFeatures()
{

}



