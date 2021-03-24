#include<iostream>
#include <complex>
#include <ctime>
#include <fstream>
#include"Eigen/Dense"
using namespace Eigen;
using namespace std;
#include "KukaKinematics.h"
const double pi = acos(-1);
int main()
{
    ofstream fout("data.txt");
    srand(time(NULL));
    double joints[7];
    double angles[7];
    double T[16], T2[16];
    Matrix3d R, R2;
    Vector3d t, t2;
    Vector3d euler, euler2;
    KukaKinematics kin;
    clock_t time1 = clock();
    for(int i = 0; i < 100000; i++)
    {
        kin.genRandomPose(T, joints);
        int ret = kin.inverseKinematics(T, angles);
        Map<Matrix4d> mT(T);
        kin.forwardKinematics(angles, T2);
        Map<Matrix4d> mT2(T2);
        R = mT.topLeftCorner(3,3);
        t = mT.col(3).topRows(3);
        euler = R.eulerAngles(2,1,0);
        R2 = mT2.topLeftCorner(3,3);
        t2 = mT2.col(3).topRows(3);
        euler2 = R2.eulerAngles(2,1,0);
        if ((euler - euler2).norm() > 1e-3 || (t - t2).norm() > 1e-3)
        {
            cout << (euler - euler2).norm() << " " << (t - t2).norm() << endl;
            fout << joints[0] << " " << joints[1] << " " << joints[2] << " " << joints[3] << " "
                 << joints[4] << " "<< joints[5] << " " << joints[6] << endl;
        }
    }
    clock_t time2 = clock();
    cout << (time2 - time1)/(double)CLOCKS_PER_SEC << " s" << endl;
    return 0;
}