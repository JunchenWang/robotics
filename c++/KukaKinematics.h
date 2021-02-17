#ifndef KUKA_KINEMATICS_H
#define KUKA_KINEMATICS_H
class KukaKinematics
{
public:
	KukaKinematics();
	~KukaKinematics();
	void setDHParameters(double _d1, double _d3, double _d5, double _d7);
	void forwardKinematics(const double* angles, double* T);
	void inverseKinematics(const double* T, double* angles);
protected:
	void updateDHTable();
private:
	double d1, d3, d5, d7;
	double dh_table[7][4];
};


#endif
