#include "KukaKinematics.h"
#include "Eigen/Dense"
#include <complex>
#include <iostream>
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
	eps1 = 1e-12; // singularity thershold for theta_1
	eps0 = 1e-7;  // singularity threshold for theta_2 and 6
	cfg[0] = 1;
	cfg[1] = -1;
	cfg[2] = 1;
	lowers[0] = -170;
	lowers[1] = -120;
	lowers[2] = -170;
	lowers[3] = -120;
	lowers[4] = -170;
	lowers[5] = -120;
	lowers[6] = -175;

	uppers[0] = 170;
	uppers[1] = 120;
	uppers[2] = 170;
	uppers[3] = 120;
	uppers[4] = 170;
	uppers[5] = 120;
	uppers[6] = 175;

	for (int i = 0; i < 7; i++)
	{
		lowers[i] *= pi / 180;
		uppers[i] *= pi / 180;
	}
	default_tol = 0.1;
	kesai = 0;
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

void KukaKinematics::forwardKinematics(const double *angles, double *T, int n)
{
	Map<Matrix4d> ret(T);
	ret = Matrix4d::Identity();
	for (int i = 0; i < n; i++)
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

void KukaKinematics::jointLimitMapping(std::vector<double> &kesai_range)
{
	Map<Matrix3d> As(mAs), Bs(mBs), Cs(mCs), Aw(mAw), Bw(mBw), Cw(mCw);
	vector<double> kesai_range_1, kesai_range_2, kesai_range_3, kesai_range_123;
	double kesai_s1 = kesai_range_hj(As(2, 2), Bs(2, 2), Cs(2, 2), cfg[0], lowers[1], uppers[1], kesai_range_2);
	kesai_range_pj(As(1, 2), Bs(1, 2), Cs(1, 2), As(0, 2), Bs(0, 2), Cs(0, 2), cfg[0], lowers[0], uppers[0], kesai_s1, kesai_range_1);
	kesai_range_pj(As(2, 1), Bs(2, 1), Cs(2, 1), -As(2, 0), -Bs(2, 0), -Cs(2, 0), cfg[0], lowers[2], uppers[2], kesai_s1, kesai_range_3);
	range_intersection(kesai_range_1, kesai_range_2, kesai_range_3, kesai_range_123);
	if (kesai_range_123.empty())
	{
		if (kesai_s1 != INT_MAX)
			kesai_range_123.push_back(kesai_s1), kesai_range_123.push_back(kesai_s1);
		else
			return;
	}

	vector<double> kesai_range_5, kesai_range_6, kesai_range_7, kesai_range_567;
	kesai_s1 = kesai_range_hj(Aw(2, 2), Bw(2, 2), Cw(2, 2), cfg[2], lowers[5], uppers[5], kesai_range_6);
	kesai_range_pj(Aw(1, 2), Bw(1, 2), Cw(1, 2), Aw(0, 2), Bw(0, 2), Cw(0, 2), cfg[2], lowers[4], uppers[4], kesai_s1, kesai_range_5);
	kesai_range_pj(Aw(2, 1), Bw(2, 1), Cw(2, 1), -Aw(2, 0), -Bw(2, 0), -Cw(2, 0), cfg[2], lowers[6], uppers[6], kesai_s1, kesai_range_7);
	range_intersection(kesai_range_5, kesai_range_6, kesai_range_7, kesai_range_567);
	if (kesai_range_567.empty())
	{
		if (kesai_s1 != INT_MAX)
			kesai_range_567.push_back(kesai_s1), kesai_range_567.push_back(kesai_s1);
		else
			return;
	}
	range_intersection(kesai_range_123, kesai_range_567, kesai_range);
}

int KukaKinematics::inverseKinematicsWithPsi(double kesai, double *angles)
{
	Map<Matrix3d> As(mAs), Bs(mBs), Cs(mCs), Aw(mAw), Bw(mBw), Cw(mCw);
	double t2 = As(2, 2) * sin(kesai) + Bs(2, 2) * cos(kesai) + Cs(2, 2);
	complex<double> z{t2, 0};
	angles[1] = cfg[0] * acos(z).real();
	if (angles[1] < lowers[1] || angles[1] > uppers[1])
		return 0;
	if (fabs(fabs(t2) - 1) < eps0)
	{
		if (t2 > 0)
		{
			double theta1and3 = atan2(As(1, 0) * sin(kesai) + Bs(1, 0) * cos(kesai) + Cs(1, 0),
									  As(0, 0) * sin(kesai) + Bs(0, 0) * cos(kesai) + Cs(0, 0));
			vector<double> bd;
			range_intersection(vector<double>{lowers[2], uppers[2]}, vector<double>{theta1and3 - uppers[0], theta1and3 - lowers[0]}, bd);
			if (bd.empty())
				return 0;
			angles[2] = (bd[0] + bd[1]) / 2;
			angles[0] = theta1and3 - angles[2];
		}
		else
		{
			double theta3_1 = atan2(As(0, 1) * sin(kesai) + Bs(0, 1) * cos(kesai) + Cs(0, 1),
									As(1, 1) * sin(kesai) + Bs(1, 1) * cos(kesai) + Cs(1, 1));
			vector<double> bd;
			range_intersection(vector<double>{lowers[2], uppers[2]}, vector<double>{theta3_1 + lowers[0], theta3_1 + uppers[0]}, bd);
			if (bd.empty())
				return 0;
			angles[2] = (bd[0] + bd[1]) / 2;
			angles[0] = angles[2] - theta3_1;
		}
	}
	else
	{
		angles[0] = atan2(cfg[0] * (As(1, 2) * sin(kesai) + Bs(1, 2) * cos(kesai) + Cs(1, 2)),
						  cfg[0] * (As(0, 2) * sin(kesai) + Bs(0, 2) * cos(kesai) + Cs(0, 2)));
		angles[2] = atan2(cfg[0] * (As(2, 1) * sin(kesai) + Bs(2, 1) * cos(kesai) + Cs(2, 1)),
						  -cfg[0] * (As(2, 0) * sin(kesai) + Bs(2, 0) * cos(kesai) + Cs(2, 0)));
		if (angles[2] < lowers[2] || angles[2] > uppers[2])
			return 0;
		if (angles[0] < lowers[0] || angles[0] > uppers[0])
			return 0;
	}
	double t6 = Aw(2, 2) * sin(kesai) + Bw(2, 2) * cos(kesai) + Cw(2, 2);
	z.real(t6);
	angles[5] = cfg[2] * acos(z).real();
	if (angles[5] < lowers[5] || angles[5] > uppers[5])
		return 0;
	if (fabs(fabs(t6) - 1) < eps0)
	{
		if (t6 > 0)
		{
			double theta5and7 = atan2(Aw(1, 0) * sin(kesai) + Bw(1, 0) * cos(kesai) + Cw(1, 0),
									  Aw(0, 0) * sin(kesai) + Bw(0, 0) * cos(kesai) + Cw(0, 0));
			vector<double> bd;
			range_intersection(vector<double>{lowers[6], uppers[6]}, vector<double>{theta5and7 - uppers[4], theta5and7 - lowers[4]}, bd);
			if (bd.empty())
				return 0;
			angles[6] = (bd[0] + bd[1]) / 2;
			angles[4] = theta5and7 - angles[6];
		}
		else
		{
			double theta7_5 = atan2(Aw(0, 1) * sin(kesai) + Bw(0, 1) * cos(kesai) + Cw(0, 1),
									Aw(1, 1) * sin(kesai) + Bw(1, 1) * cos(kesai) + Cw(1, 1));
			vector<double> bd;
			range_intersection(vector<double>{lowers[6], uppers[6]}, vector<double>{theta7_5 + lowers[4], theta7_5 + uppers[4]}, bd);
			if (bd.empty())
				return 0;
			angles[6] = (bd[0] + bd[1]) / 2;
			angles[4] = angles[6] - theta7_5;
		}
	}
	else
	{
		angles[4] = atan2(cfg[2] * (Aw(1, 2) * sin(kesai) + Bw(1, 2) * cos(kesai) + Cw(1, 2)),
						  cfg[2] * (Aw(0, 2) * sin(kesai) + Bw(0, 2) * cos(kesai) + Cw(0, 2)));
		angles[6] = atan2(cfg[2] * (Aw(2, 1) * sin(kesai) + Bw(2, 1) * cos(kesai) + Cw(2, 1)),
						  -cfg[2] * (Aw(2, 0) * sin(kesai) + Bw(2, 0) * cos(kesai) + Cw(2, 0)));
		if (angles[6] < lowers[6] || angles[6] > uppers[6])
			return 0;
		if (angles[4] < lowers[4] || angles[4] > uppers[4])
			return 0;
	}
	return 1;
}

int KukaKinematics::inverseKinematics(const double *T, double *angles)
{
	angles[3] = calABCMatrix(T);
	if (angles[3] < lowers[3] || angles[3] > uppers[3])
		return 0;
	vector<double> kesai_range;
	jointLimitMapping(kesai_range);
	if (kesai_range.empty())
		return 0;
	return inverseKinematicsWithPsi(choosePsi(kesai_range), angles);
}

void KukaKinematics::setCfg(int gc2, int gc4, int gc6)
{
	cfg[0] = gc2;
	cfg[1] = gc4;
	cfg[2] = gc6;
}

void KukaKinematics::genRandomPose(double *T, double *angles)
{
	double tol = 0.05;
	double s;
	for (int i = 0; i < 7; i++)
	{
		if (i == 3)
			continue;
		s = (double)rand() / RAND_MAX;
		angles[i] = lowers[i] + tol + s * (uppers[i] - lowers[i] - 2 * tol);
	}
	s = (double)rand() / RAND_MAX;
	if (s > 0.5)
		angles[3] = tol + (s - 0.5) / 0.5 * (uppers[3] - 2*tol);
	else
		angles[3] = -tol + (s - 0.5) / 0.5 * (-lowers[3] - 2*tol);
	cfg[0] = angles[1] > 0 ? 1 : -1;
	cfg[1] = angles[3] > 0 ? 1 : -1;
	cfg[2] = angles[5] > 0 ? 1 : -1;
	forwardKinematics(angles, T);
}

double KukaKinematics::choosePsi(const vector<double> &kesai_range)
{
	double j = 0, len = kesai_range[1] - kesai_range[0];
	for (int i = 1; i < kesai_range.size() / 2; i++)
	{
		if (kesai_range[2 * i + 1] - kesai_range[2 * i] > len)
		{
			len = kesai_range[2 * i + 1] - kesai_range[2 * i];
			j = i;
		}
	}
	return (kesai_range[2 * j] + kesai_range[2 * j + 1]) / 2;
}

void KukaKinematics::updateDHTable()
{
	double table[][4] = {0, 0, d1, 0,
						 -pi / 2, 0, 0, 0,
						 pi / 2, 0, d3, 0,
						 pi / 2, 0, 0, 0,
						 -pi / 2, 0, d5, 0,
						 -pi / 2, 0, 0, 0,
						 pi / 2, 0, d7, 0};
	memcpy(dh_table, table, 28 * sizeof(double));
}

double KukaKinematics::calABCMatrix(const double *T)
{
	Map<const Matrix4d> Th(T);
	Matrix3d R = Th.topLeftCorner(3, 3);
	Vector3d t = Th.col(3).topRows(3);
	Vector3d z(0, 0, 1);
	Vector3d p02(0, 0, d1);
	Vector3d p67(0, 0, d7);
	Vector3d p26 = t - p02 - R * p67;
	double l_p26 = p26.norm();
	double l2_p26 = p26.squaredNorm();
	Vector3d p26_hat = p26.normalized();
	double theta3 = 0;
	complex<double> zz{(l2_p26 - d3 * d3 - d5 * d5) / (2 * d3 * d5), 0};
	double theta4 = cfg[1] * acos(zz).real();
	Matrix3d R34;
	double theta1;
	R34 << cos(theta4), 0, -sin(theta4), 0, 1, 0, sin(theta4), 0, cos(theta4);
	if (fabs(fabs(p26_hat.dot(z)) - 1) < eps1)
		theta1 = 0;
	else
		theta1 = atan2(p26(1), p26(0));
	zz.real((d3 * d3 + l2_p26 - d5 * d5) / (2 * d3 * l_p26));
	double phi = acos(zz).real();
	double theta2 = atan2(sqrt(p26(0) * p26(0) + p26(1) * p26(1)), p26(2)) + cfg[1] * phi;
	Matrix4d T03;
	double angles[] = {theta1, theta2, theta3};
	forwardKinematics(angles, T03.data(), 3);
	Matrix3d R03 = T03.topLeftCorner(3, 3);
	Map<Matrix3d> As(mAs), Bs(mBs), Cs(mCs), Aw(mAw), Bw(mBw), Cw(mCw);
	Matrix3d p26_hat_skew;
	p26_hat_skew << 0, -p26_hat(2), p26_hat(1),
		p26_hat(2), 0, -p26_hat(0),
		-p26_hat(1), p26_hat(0), 0;
	As = p26_hat_skew * R03;
	Bs = -p26_hat_skew * As;
	Cs = p26_hat * p26_hat.transpose() * R03;
	Matrix3d R34_tr = R34.transpose();
	Aw = R34_tr * As.transpose() * R;
	Bw = R34_tr * Bs.transpose() * R;
	Cw = R34_tr * Cs.transpose() * R;
	return theta4;
}

double KukaKinematics::kesai_range_hj(double a, double b, double c, int cfg, double lower, double upper, std::vector<double> &kesai_range)
{
	double kesai_s = INT_MAX;
	if (a == 0 && b == 0)
	{
		double theta = cfg * acos(complex<double>{c, 0}).real();
		if (theta >= lower && theta <= upper)
			kesai_range.push_back(-pi), kesai_range.push_back(pi);
	}
	else
	{
		double kesai1 = atan2(a, b);
		double t1 = a * sin(kesai1) + b * cos(kesai1) + c;
		double theta1 = cfg * acos(complex<double>{t1, 0}).real();
		double kesai2 = atan2(-a, -b);
		double t2 = a * sin(kesai2) + b * cos(kesai2) + c;
		double theta2 = cfg * acos(complex<double>{t2, 0}).real();
		double thr1 = abs(abs(t1) - 1);
		double thr2 = abs(abs(t2) - 1);
		if (thr1 < thr2 && thr1 < eps0)
			kesai_s = kesai1;
		else if (thr2 < thr1 && thr2 < eps0)
			kesai_s = kesai2;
		double l = min(theta1, theta2);
		double u = max(theta1, theta2);

		if (upper <= l || lower >= u)
			return kesai_s;

		vector<double> bdu{-pi, pi}, bdl{-pi, pi};
		double ks1, ks2;
		if (upper < u)
		{
			bdu.clear();
			theta2kesai_hj(upper, a, b, c, ks1, ks2);
			if ((b * sin(ks1) - a * cos(ks1)) * upper > 0)
			{
				bdu.push_back(-pi);
				bdu.push_back(ks1);
				bdu.push_back(ks2);
				bdu.push_back(pi);
			}
			else
			{
				bdu.push_back(ks1);
				bdu.push_back(ks2);
			}
		}
		if (lower > l)
		{
			bdl.clear();
			theta2kesai_hj(lower, a, b, c, ks1, ks2);
			if ((b * sin(ks1) - a * cos(ks1)) * lower > 0)
			{
				bdl.push_back(ks1);
				bdl.push_back(ks2);
			}
			else
			{
				bdl.push_back(-pi);
				bdl.push_back(ks1);
				bdl.push_back(ks2);
				bdl.push_back(pi);
			}
		}
		range_intersection(bdl, bdu, kesai_range);
	}
	return kesai_s;
}

void KukaKinematics::theta2kesai_hj(double theta, double a, double b, double c, double &kesai1, double &kesai2)
{
	double delta = a * a + b * b - pow(c - cos(theta), 2);
	kesai1 = 2 * atan((a - sqrt(delta)) / (cos(theta) + b - c));
	kesai2 = 2 * atan((a + sqrt(delta)) / (cos(theta) + b - c));
	if (kesai1 > kesai2)
	{
		double tem = kesai1;
		kesai1 = kesai2;
		kesai2 = tem;
	}
}

double KukaKinematics::kesai2theta_pj(double kesai, double an, double bn, double cn, double ad, double bd, double cd, int cfg)
{
	return atan2(cfg * (an * sin(kesai) + bn * cos(kesai) + cn),
				 cfg * (ad * sin(kesai) + bd * cos(kesai) + cd));
}

void KukaKinematics::theta2kesai_pj_2(double theta, double an, double bn, double cn, double ad, double bd, double cd, double &ks1, double &ks2)
{
	double ap = (cd - bd) * tan(theta) + (bn - cn);
	double bp = 2 * (ad * tan(theta) - an);
	double cp = (bd + cd) * tan(theta) - (bn + cn);
	double delta2 = bp * bp - 4 * ap * cp;
	ks1 = 2 * atan((-bp - sqrt(delta2)) / (2 * ap));
	ks2 = 2 * atan((-bp + sqrt(delta2)) / (2 * ap));
	if (ks1 > ks2)
	{
		double tem = ks1;
		ks1 = ks2;
		ks2 = tem;
	}
}

double KukaKinematics::theta2kesai_pj_1(double theta, double an, double bn, double cn, double ad, double bd, double cd, int cfg)
{
	double ap = (cd - bd) * tan(theta) + (bn - cn);
	double bp = 2 * (ad * tan(theta) - an);
	double cp = (bd + cd) * tan(theta) - (bn + cn);
	if (abs(ap) < abs(bp) * 1e-12)
		return 2 * atan(-cp / bp);

	double delta2 = bp * bp - 4 * ap * cp;
	double kesai1 = 2 * atan((-bp - sqrt(delta2)) / (2 * ap));
	double kesai2 = 2 * atan((-bp + sqrt(delta2)) / (2 * ap));
	double e1 = abs(kesai2theta_pj(kesai1, an, bn, cn, ad, bd, cd, cfg) - theta);
	double e2 = abs(kesai2theta_pj(kesai2, an, bn, cn, ad, bd, cd, cfg) - theta);
	if (e1 < e2)
		return kesai1;
	else
		return kesai2;
}

double KukaKinematics::theta2kesai_pj_1_s(double theta, double an, double bn, double cn, double ad, double bd, double cd, double kesai_s)
{
	double ap = (cd - bd) * tan(theta) + (bn - cn);
	double bp = 2 * (ad * tan(theta) - an);
	double cp = (bd + cd) * tan(theta) - (bn + cn);
	if (abs(ap) < abs(bp) * 1e-12)
		return 2 * atan(-cp / bp);

	double delta2 = bp * bp - 4 * ap * cp;
	double kesai1 = 2 * atan((-bp - sqrt(delta2)) / (2 * ap));
	double kesai2 = 2 * atan((-bp + sqrt(delta2)) / (2 * ap));
	double d1 = abs(kesai1 - kesai_s);
	double d2 = abs(kesai2 - kesai_s);
	double th1 = min(d1, 2 * pi - d1);
	double th2 = min(d2, 2 * pi - d2);
	if (th1 > th2)
		return kesai1;
	else
		return kesai2;
}

void range_intersection(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> &c)
{
	for (int i = 0; i < a.size() / 2; i++)
		for (int j = 0; j < b.size() / 2; j++)
		{
			if (a[2 * i + 1] < b[2 * j] || b[2 * j + 1] < a[2 * i])
				continue;
			c.push_back(max(a[2 * i], b[2 * j]));
			c.push_back(min(a[2 * i + 1], b[2 * j + 1]));
		}
}

void KukaKinematics::kesai_range_pj(double an, double bn, double cn,
									double ad, double bd, double cd,
									int cfg, double lower, double upper,
									double kesai_s, std::vector<double> &kesai_range)
{
	auto y = [=](double kesai) {
		return kesai2theta_pj(kesai, an, bn, cn, ad, bd, cd, cfg);
	};
	auto k2 = [=](double theta, double &k1, double &k2) {
		theta2kesai_pj_2(theta, an, bn, cn, ad, bd, cd, k1, k2);
	};
	auto k1 = [=](double theta) {
		return theta2kesai_pj_1(theta, an, bn, cn, ad, bd, cd, cfg);
	};
	auto ks1 = [=](double theta) {
		return theta2kesai_pj_1_s(theta, an, bn, cn, ad, bd, cd, kesai_s);
	};

	double at = cn * bd - bn * cd;
	double bt = an * cd - cn * ad;
	double ct = an * bd - bn * ad;
	if (at == 0 && bt == 0 && ct == 0)
	{
		double theta = atan2(cfg * (bn + cn), cfg * (bd + cd));
		if (upper >= theta && lower <= theta)
			kesai_range.push_back(-pi), kesai_range.push_back(pi);
		return;
	}
	auto d = [=](double kesai) {
		return at * sin(kesai) + bt * cos(kesai) + ct;
	};
	double tol = default_tol;
	double delta = (at * at + bt * bt - ct * ct);
	double delta_n = delta / (at * at + bt * bt + ct * ct);
	if (kesai_s != INT_MAX)
	{
		if (abs(delta_n) > 1e-3 && bt != ct)
		{
			complex<double> z{delta, 0};
			complex<double> kesai1 = 2.0 * atan((at + sqrt(z)) / (bt - ct));
			complex<double> kesai2 = 2.0 * atan((at - sqrt(z)) / (bt - ct));
			tol = max(default_tol, min(abs(kesai2 - kesai1), 2 * pi - abs(kesai2 - kesai1)));
		}
		double kesai1 = kesai_s - tol;
		double kesai2 = kesai_s + tol;
		double theta1 = y(kesai1);
		double theta2 = y(kesai2);
		double l = min(theta1, theta2);
		double u = max(theta1, theta2);
		vector<double> bd_s, bdu, bdl, bd;
		if (kesai1 < -pi)
			bd_s.push_back(kesai2), bd_s.push_back(2 * pi + kesai1);
		else if (kesai2 > pi)
			bd_s.push_back(kesai2 - 2 * pi), bd_s.push_back(kesai1);
		else
		{
			bd_s.push_back(-pi), bd_s.push_back(kesai1);
			bd_s.push_back(kesai2), bd_s.push_back(pi);
		}
		double d1 = d(kesai1);
		double kesai_u, kesai_l;
		if ((theta1 - theta2) * d1 > 0) //N type
		{
			if (upper <= l || lower >= u)
				return;
			if (upper >= u || upper <= l)
				kesai_u = kesai_s;
			else
				kesai_u = ks1(upper);
			if (lower <= l || lower >= u)
				kesai_l = kesai_s;
			else
				kesai_l = ks1(lower);
			if (theta1 > theta2)
			{
				bdu.push_back(-pi), bdu.push_back(kesai_u);
				bdl.push_back(kesai_l), bdl.push_back(pi);
			}
			else
			{
				bdu.push_back(kesai_u), bdu.push_back(pi);
				bdl.push_back(-pi), bdl.push_back(kesai_l);
			}
			if ((kesai_u - kesai_l) * (theta1 - theta2) > 0) // % not >= !!
			{
				range_intersection(bdu, bdl, bd);
			}
			else
			{
				bd.insert(bd.end(), bdl.begin(), bdl.end());
				bd.insert(bd.end(), bdu.begin(), bdu.end());
			}
		}
		else // H type
		{
			if (upper <= u && upper >= l)
				kesai_u = kesai_s;
			else
				kesai_u = ks1(upper);

			if (lower >= l && lower <= u)
				kesai_l = kesai_s;
			else
				kesai_l = ks1(lower);

			if (theta1 > theta2)
			{
				bdu.push_back(kesai_u), bdu.push_back(pi);
				bdl.push_back(-pi), bdl.push_back(kesai_l);
			}
			else
			{
				bdu.push_back(-pi), bdu.push_back(kesai_u);
				bdl.push_back(kesai_l), bdl.push_back(pi);
			}

			if ((kesai_u - kesai_l) * (theta1 - theta2) <= 0)
				range_intersection(bdu, bdl, bd);
			else
			{
				bd.insert(bd.end(), bdl.begin(), bdl.end());
				bd.insert(bd.end(), bdu.begin(), bdu.end());
			}
		}
		range_intersection(bd, bd_s, kesai_range);
	}
	else if (delta > 0)
	{
		double kesai1 = 2 * atan((at + sqrt(delta)) / (bt - ct));
		double kesai2 = 2 * atan((at - sqrt(delta)) / (bt - ct));
		if (at * cos(kesai1) - bt * sin(kesai1) < 0) // % local maxmia
		{
			double tem = kesai2;
			kesai2 = kesai1;
			kesai1 = tem;
		}
		double y1 = y(kesai1);
		double y2 = y(kesai2);
		if (y1 == y2)
		{
			if (upper >= y1 && lower <= y1)
				kesai_range.push_back(-pi), kesai_range.push_back(pi);
			return;
		}
		if (y1 < y2) // no jump
		{
			vector<double> bdl{-pi, pi}, bdu{-pi, pi};
			if (lower >= y2)
				return;
			else if (lower < y2 && lower > y1)
			{
				k2(lower, kesai1, kesai2);
				bdl.clear();
				if (d(kesai1) < 0)
				{
					bdl.push_back(-pi), bdl.push_back(kesai1);
					bdl.push_back(kesai2), bdl.push_back(pi);
				}
				else
					bdl.push_back(kesai1), bdl.push_back(kesai2);
			}

			if (upper <= y1)
				return;
			else if (upper < y2 && upper > y1)
			{
				k2(upper, kesai1, kesai2);
				bdu.clear();
				if (d(kesai1) < 0)
				{
					bdu.push_back(kesai1), bdu.push_back(kesai2);
				}
				else
				{
					bdu.push_back(-pi), bdu.push_back(kesai1);
					bdu.push_back(kesai2), bdu.push_back(pi);
				}
			}

			range_intersection(bdl, bdu, kesai_range);
		}
		else // with jump
		{
			double ks1, ks2;
			vector<double> bdl, bdu;
			if (upper > y1 || upper < y2)
			{
				k2(upper, ks1, ks2);
				if ((ks1 - kesai1) * (ks2 - kesai1) > 0)
				{
					bdu.push_back(-pi), bdu.push_back(ks1);
					bdu.push_back(ks2), bdu.push_back(pi);
				}
				else
					bdu.push_back(ks1), bdu.push_back(ks2);
			}
			if (lower > y1 || lower < y2)
			{
				k2(lower, ks1, ks2);
				if ((ks1 - kesai2) * (ks2 - kesai2) > 0)
				{
					bdl.push_back(-pi), bdl.push_back(ks1);
					bdl.push_back(ks2), bdl.push_back(pi);
				}
				else
					bdl.push_back(ks1), bdl.push_back(ks2);
			}
			if (lower > y1 || upper < y2)
				range_intersection(bdl, bdu, kesai_range);
			else
			{
				kesai_range.insert(kesai_range.end(), bdl.begin(), bdl.end());
				kesai_range.insert(kesai_range.end(), bdu.begin(), bdu.end());
			}
		}
	}
	else
	{
		double kl = k1(lower);
		double ku = k1(upper);
		double kesai1 = min(kl, ku);
		double kesai2 = max(kl, ku);
		double theta = y((kesai1 + kesai2) / 2);
		if (theta > lower && theta < upper) // % no jump
			kesai_range.push_back(kesai1), kesai_range.push_back(kesai2);
		else // % with jump
		{
			kesai_range.push_back(-pi), kesai_range.push_back(kesai1);
			kesai_range.push_back(kesai2), kesai_range.push_back(pi);
		}
	}
}

void range_intersection(const std::vector<double> &a, const std::vector<double> &b, const std::vector<double> &c, std::vector<double> &d)
{
	for (int i = 0; i < a.size() / 2; i++)
		for (int j = 0; j < b.size() / 2; j++)
			for (int k = 0; k < c.size() / 2; k++)
			{
				if (a[2 * i + 1] < b[2 * j] || b[2 * j + 1] < a[2 * i])
					continue;
				double t1 = max(a[2 * i], b[2 * j]);
				double t2 = min(a[2 * i + 1], b[2 * j + 1]);
				if (c[2 * k + 1] < t1 || c[2 * k] > t2)
					continue;
				d.push_back(max(c[2 * k], t1));
				d.push_back(min(c[2 * k + 1], t2));
			}
}
