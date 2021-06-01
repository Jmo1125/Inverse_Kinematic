#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace mr {
Eigen::Vector3d Rotate2Eular(const Eigen::MatrixXd T_final);

Eigen::MatrixXd FKinSpace(const Eigen::VectorXd& thetaList,  Eigen::MatrixXd& Jacob);

Eigen::MatrixXd IKinSpace(const Eigen::MatrixXd& T_final, const Eigen::VectorXd& Initial_theta);

Eigen::MatrixXd pinv(const Eigen::MatrixXd&  A);

Eigen::MatrixXd JacobSpace(const Eigen::MatrixXd& T01, const Eigen::MatrixXd& T12, const Eigen::MatrixXd& T23, const Eigen::MatrixXd& T34, const Eigen::MatrixXd& T45, const Eigen::MatrixXd& T56);

bool Is_OK(const Eigen::VectorXd& diffe);

void Normalize(Eigen::VectorXd& angle);
}
