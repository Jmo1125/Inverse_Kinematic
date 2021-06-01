#include "test.h"
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Core>

//template<typename _Matrix_Type_>
//_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon =
//    std::numeric_limits<double>::epsilon())
//{
//    Eigen::JacobiSVD< _Matrix_Type_ > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
//    double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
//    return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
//}

using namespace mr;
using namespace std;
int main()
{
    Eigen::MatrixXd Jacob(6,6);
    Eigen::Vector3d Eular;
    Eigen::VectorXd thetalist(6, 1); //当前关节角度值
    thetalist << 0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
    cout << "initial joint = " << endl;
    cout << thetalist << endl;
    cout << "Initial = " << endl;
    cout << FKinSpace(thetalist, Jacob) << endl;

//    Eigen::MatrixXd Initial(4, 4); //初始的位姿矩阵
//    Initial <<  1, 0, 0, 0,
//                0, 1, 0, 0,
//                0, 0, 1, 0,
//                0, 0, 0, 1;
//    cout << "Initial = " << endl;
//    cout << Initial << endl;

    Eigen::VectorXd TargetTheta(6, 1); //目标关节角度值
    TargetTheta << 0.2,0.2,0.2,0.2,0.2,0.2;
    Eigen::MatrixXd Tar(4, 4); //目标位姿矩阵
    Tar = FKinSpace(TargetTheta, Jacob);
    cout << "Tar = " << endl;
    cout << Tar << endl;

    thetalist = IKinSpace(Tar, thetalist);
    Eigen::MatrixXd inv(4, 4);
    inv = FKinSpace(thetalist, Jacob);
    cout << "FKinSpace = " << endl;
    cout << inv << endl;
    cout << "Error = " << endl;
    cout << Tar - inv << endl;
    cout << "joint result = " << endl;
    cout << thetalist << endl;

}

