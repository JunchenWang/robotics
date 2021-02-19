#include<iostream>
#include <complex>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
#include "KukaKinematics.h"
const double pi = acos(-1);
int main()
{
    Matrix3d R;
    R <<    -0.6612 ,  -0.4121 ,  -0.6269,
   -0.6742 ,  -0.0400 ,   0.7375,
   -0.3290 ,   0.9103,   -0.2513;
    Vector3d t(200, 200, 600);
    Matrix4d T = Matrix4d::Identity();
    Matrix<double, 7, 1> angles;
    T.topLeftCorner(3,3) = R;
    T.col(3).topRows(3) = t;
    KukaKinematics kin;
    kin.calABCMatrix(T.data());
    kin.inverseKinematics(0.1, angles.data());
    cout << angles << endl;
    return 0;
}