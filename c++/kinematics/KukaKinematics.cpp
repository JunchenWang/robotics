#include "KukaKinematics.h"
#include "Eigen/Dense"
#include <complex>
//#define _USE_MATH_DEFINES
//#include <math.h>
using namespace Eigen;
using namespace std;

const double pi = acos(-1);
KukaKinematics::KukaKinematics()
{
	d1 = 340;
	d3 = d5 = 400;
	d7 = 126;
	eps1 = 1e-12;
	eps0 = 1e-7;
	updateDHTable();
}

KukaKinematics::~KukaKinematics()
{
}

void KukaKinematics::setDHParameters(double _d1, double _d3, double _d5, double _d7)
{
	d1 = _d1;
	d3 = _d3;
	d5 = _d5;
	d7 = _d7;
	updateDHTable();
}

void KukaKinematics::forwardKinematics(const double* angles, double* T)
{
	Map<Matrix4d> ret(T);
	ret = Matrix4d::Identity();
	for (int i = 0; i < 7; i++)
	{
		double ct = cos(angles[i]);
		double st = sin(angles[i]);
		double cpha = cos(dh_table[i][0]);
		double spha = sin(dh_table[i][0]);
		Matrix4d m;
		m << ct, -st, 0, dh_table[i][1],
			st * cpha, ct * cpha, -spha, -dh_table[i][2] * spha,
			st * spha, ct * spha, cpha, dh_table[i][2] * cpha,
			0, 0, 0, 1;
		ret = ret * m;
	}
}

void KukaKinematics::inverseKinematics(const double* T, double* angles)
{
}

void KukaKinematics::updateDHTable()
{
	double table[][4] = { 0, 0, d1, 0,
						 -pi / 2, 0, 0, 0,
						 pi / 2, 0, d3, 0,
						 pi / 2, 0, 0, 0,
						 -pi / 2, 0, d5, 0,
						 -pi / 2, 0, 0, 0,
						 pi / 2, 0, d7, 0 };
	memcpy(dh_table, table, 28 * sizeof(double));
}
