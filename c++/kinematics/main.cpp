#include<iostream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
#include "KukaKinematics.h"
const double pi = acos(-1);
int main()
{
    Matrix<double, 7, 1> angle;
    angle << 20, 40, 30, 50, 10, 20, 45;
    angle = angle / 180 * pi;
    Matrix4d T;
    KukaKinematics kin;
    kin.forwardKinematics(angle.data(), T.data());
    cout << T << endl;
    return 0;
}