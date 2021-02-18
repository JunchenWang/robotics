#ifndef KUKA_KINEMATICS_H
#define KUKA_KINEMATICS_H
// matrix is in row-major order
#include <vector>
class KukaKinematics
{
public:
	KukaKinematics();
	~KukaKinematics();
	void setDHParameters(double _d1, double _d3, double _d5, double _d7);
	void forwardKinematics(const double *angles, double *T);
	void jointLimitMapping(const double *T, std::vector<double> kesai_range);
	void inverseKinematics(double kesai, double *angles);
	void inverseKinematics(const double *T, double *angles);
	void setCfg(int gc1, int gc2, int gc3);
protected:
	void updateDHTable();
private:
	double d1, d3, d5, d7;
	double dh_table[7][4];
	double As[3][3], Bs[3][3], Cs[3][3];
	double Aw[3][3], Bw[3][3], Cw[3][3];
	double eps1, eps0;
	int cfg[3];
};

#endif
