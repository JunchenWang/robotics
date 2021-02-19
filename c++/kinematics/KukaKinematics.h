#ifndef KUKA_KINEMATICS_H
#define KUKA_KINEMATICS_H
// matrix is in colum-major order
#include <vector>

void range_intersection(const std::vector<double> & a, const std::vector<double> & b, std::vector<double> & c);
class KukaKinematics
{
public:
	KukaKinematics();
	~KukaKinematics();
	void setDHParameters(double _d1, double _d3, double _d5, double _d7);
	void forwardKinematics(const double *angles, double *T, int n);
	void jointLimitMapping(std::vector<double> kesai_range);
	int inverseKinematics(double kesai, double *angles);
	int inverseKinematics(const double *T, double *angles);
	double calABCMatrix(const double *T);
	void setCfg(int gc2, int gc4, int gc6);
protected:
	void updateDHTable();
	double kesai_range_hj(double a, double b, double c, int cfg, double lower, double upper, std::vector<double> & kesai_range);
	void theta2kesai_hj(double theta, double a, double b, double c, double & kesai1, double & kesai2);
	double kesai2theta_pj(double kesai, double an, double bn, double cn, double ad, double bd, double cd, int cfg);
	void theta2kesai_pj_2(double theta, double an, double bn, double cn, double ad, double bd, double cd, double & ks1, double & ks2);
	double theta2kesai_pj_1(double theta, double an, double bn, double cn, double ad, double bd, double cd, int cfg);
	double theta2kesai_pj_1_s(double theta, double an, double bn, double cn, double ad, double bd, double cd, double kesai_s);
	void kesai_range_pj(double an, double bn, double cn,
	                    double ad, double bd, double cd, 
						int cfg, double lower, double upper, 
						double kesai_s, std::vector<double> & kesai_range);
private:
	double d1, d3, d5, d7;
	double dh_table[7][4];
	double mAs[9], mBs[9], mCs[9];
	double mAw[9], mBw[9], mCw[9];
	double eps1, eps0;
	int cfg[3];
	double default_tol;
	double theta[7];
	double lowers[7];
	double uppers[7];
};

#endif
