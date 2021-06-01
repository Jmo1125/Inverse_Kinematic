#include "test.h"
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <vector>

namespace mr {
    void Normalize(Eigen::VectorXd& angle){
        for(int i=0;i<angle.size();i++){
            while(angle[i]>M_PI || angle[i]<-M_PI){
                if(angle[i]>M_PI) angle[i]-=M_PI;
                if(angle[i]<-M_PI) angle[i]+=M_PI;
            }
        }
    }
    Eigen::MatrixXd FKinSpace(const Eigen::VectorXd& thetalist, Eigen::MatrixXd& Jacob) {
            Eigen::VectorXd thetal = thetalist;
            Eigen::MatrixXd T_final;

            double theta2 = thetal(1) + M_PI/2.0;
            double d3 = 5.021;
            double theta5 = thetal(4) - M_PI/2.0;
            double theta6 = thetal(5);

            Eigen::MatrixXd DH(6,4);
            DH << 0, 0, thetal(0), 0,
                0.2, M_PI/2.0, theta2, 0,
                0.26, M_PI/2.0, thetal(2), d3,
                0, 0, thetal(3), 0,
                0, -M_PI/2.0, theta5, -0.632,
                0.145, -M_PI/2.0, theta6, 0.2575;

            Eigen::VectorXd a = DH.col(0);
            Eigen::VectorXd alpha = DH.col(1);
            Eigen::VectorXd theta = DH.col(2);
            Eigen::VectorXd d = DH.col(3);

            Eigen::Matrix4d T01, T12, T23, T34, T45, T56;
            int i = 0;
            T01 << cos(theta(i)), -sin(theta(i)), 0, a(i),
                    sin(theta(i))*cos(alpha(i)), cos(theta(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i),
                    sin(theta(i))*sin(alpha(i)), cos(theta(i))*sin(alpha(i)), cos(alpha(i)), cos(alpha(i))*d(i),
                    0,0,0,1;
            i = i + 1;
            T12 << cos(theta(i)), -sin(theta(i)), 0, a(i),
                    sin(theta(i))*cos(alpha(i)), cos(theta(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i),
                    sin(theta(i))*sin(alpha(i)), cos(theta(i))*sin(alpha(i)), cos(alpha(i)), cos(alpha(i))*d(i),
                    0,0,0,1;
            i = i + 1;
            T23 << cos(theta(i)), -sin(theta(i)), 0, a(i),
                    sin(theta(i))*cos(alpha(i)), cos(theta(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i),
                    sin(theta(i))*sin(alpha(i)), cos(theta(i))*sin(alpha(i)), cos(alpha(i)), cos(alpha(i))*d(i),
                    0,0,0,1;
            i = i + 1;
            T34 << cos(theta(i)), -sin(theta(i)), 0, a(i),
                    sin(theta(i))*cos(alpha(i)), cos(theta(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i),
                    sin(theta(i))*sin(alpha(i)), cos(theta(i))*sin(alpha(i)), cos(alpha(i)), cos(alpha(i))*d(i),
                    0,0,0,1;
            i = i + 1;
            T45 << cos(theta(i)), -sin(theta(i)), 0, a(i),
                    sin(theta(i))*cos(alpha(i)), cos(theta(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i),
                    sin(theta(i))*sin(alpha(i)), cos(theta(i))*sin(alpha(i)), cos(alpha(i)), cos(alpha(i))*d(i),
                    0,0,0,1;
            i = i + 1;
            T56 << cos(theta(i)), -sin(theta(i)), 0, a(i),
                    sin(theta(i))*cos(alpha(i)), cos(theta(i))*cos(alpha(i)), -sin(alpha(i)), -sin(alpha(i))*d(i),
                    sin(theta(i))*sin(alpha(i)), cos(theta(i))*sin(alpha(i)), cos(alpha(i)), cos(alpha(i))*d(i),
                    0,0,0,1;

            T_final = T01 * T12 * T23 * T34 * T45 * T56;
            Jacob = JacobSpace(T01, T12, T23, T34, T45, T56);

            return T_final;
    }

    Eigen::Vector3d Rotate2Eular(const Eigen::MatrixXd T_final){
        Eigen::Vector3d rotate(3);
        double roll = 0, pitch = 0, yaw = 0;
        pitch = atan2(-T_final(2,0), sqrt(T_final(2,1) * T_final(2,1) + T_final(2,2) * T_final(2,2)));
        if(pitch >= -M_PI/2.0||pitch <= M_PI/2.0){
            roll = atan2(T_final(2,1), T_final(2,2));
            yaw = atan2(T_final(1,0), T_final(0,0));
        }
        else{
            roll = atan2(-T_final(2,1), -T_final(2,2));
            yaw = atan2(-T_final(1,0), -T_final(0,0));
        }
        rotate << roll, pitch, yaw;
        return rotate;
    }

    Eigen::MatrixXd IKinSpace(const Eigen::MatrixXd& T_final, const Eigen::VectorXd& Initial_theta){
        Eigen::MatrixXd T = T_final;
        Eigen::MatrixXd Final;
        Eigen::VectorXd diffe(6);
        Eigen::MatrixXd Jacob(6,6), Jacob_Ni(6,6);
        Eigen::Vector3d Eular;
        Eigen::VectorXd thetalist(6), dq;

        Eigen::Vector3d EularTarget = Rotate2Eular(T);
        //Eigen::Vector3d EularTarget = T.block<3,3>(0,0).eulerAngles(0,1,2);
        Final = FKinSpace(Initial_theta, Jacob);// the position of now
        Jacob_Ni = pinv(Jacob);
        //Jacob_Ni = Jacob.inverse();
        //Eular = Final.block<3,3>(0,0).eulerAngles(2,1,0);

        Eular = Rotate2Eular(Final);
        diffe.block<3, 1>(0, 0) = (T - Final).block<3, 1>(0, 3);
        diffe.block<3, 1>(3, 0) = EularTarget - Eular;
//        F = pinv(Final) * T - Eigen::Matrix4d::Identity();
//        std::cout << F << std::endl;
//        diffe.block<3,1>(3,0);
        dq = Jacob_Ni * diffe;
        std::cout << "diffe" << std::endl;
        std::cout << diffe << std::endl;
        std::cout << "dq" << std::endl;
        std::cout << dq << std::endl;
        thetalist = dq + Initial_theta;

        int num = 0;
        bool flag = false;

        while(!flag && num < 5000){
            num++;
            Final = FKinSpace(thetalist, Jacob);
            //Eular = Final.block<3,3>(0,0).eulerAngles(0,1,2);
            Eular = Rotate2Eular(Final);
            Jacob_Ni = pinv(Jacob);
           // Jacob_Ni = Jacob.inverse();
            diffe.block<3, 1>(0, 0) = (T - Final).block<3, 1>(0, 3);
            diffe.block<3, 1>(3, 0) = EularTarget - Eular;
            flag = Is_OK(diffe);

            dq = Jacob_Ni * diffe;
            thetalist = (dq + thetalist);
            Normalize(thetalist);
        }

        std::cout<<flag<<std::endl;

        return thetalist;
    }

/* Function: Get the generalized inverse matrix
 * Inputs: the matrix A
 * Returns: the generalized inverse matrix X
 */
    Eigen::MatrixXd pinv(const Eigen::MatrixXd&  A){
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
        double  pinvtoler = 1.e-8; //tolerance
        int row = A.rows();
        int col = A.cols();
        int k = std::min(row,col);
        Eigen::MatrixXd X = Eigen::MatrixXd::Zero(col,row);
        Eigen::MatrixXd singularValues_inv = svd.singularValues();//奇异值
        Eigen::MatrixXd singularValues_inv_mat = Eigen::MatrixXd::Zero(col, row);

        for (long i = 0; i < k; ++i) {
            if (singularValues_inv(i) > pinvtoler)
                singularValues_inv(i) = 1.0 / singularValues_inv(i);
            else singularValues_inv(i) = 0;
        }
        for (long i = 0; i < k; ++i)
        {
            singularValues_inv_mat(i, i) = singularValues_inv(i);
        }
        X=(svd.matrixV())*(singularValues_inv_mat)*(svd.matrixU().transpose());

        return X;
    }

/* Function: Get the Jacob matrix
 * Inputs: Each homogeneous transformation matrix
 * Returns: The Jacob matrix
 */
    Eigen::MatrixXd JacobSpace(const Eigen::MatrixXd& T01, const Eigen::MatrixXd& T12, const Eigen::MatrixXd& T23, const Eigen::MatrixXd& T34, const Eigen::MatrixXd& T45, const Eigen::MatrixXd& T56){
        Eigen::MatrixXd Jacob(6,6);
        Eigen::MatrixXd T(6,6);
        Eigen::MatrixXd Final_T = T01 * T12 * T23 * T34 * T45 * T56;
        Eigen::Vector3d pN = Final_T.block<3,1>(0,3);

        T = T01;
        Jacob.block<3,1>(0,0) = T.block<3,1>(0,2).cross(pN - T.block<3,1>(0,3));
        Jacob.block<3,1>(3,0) = T.block<3,1>(0,2);
        T = T*T12;
        Jacob.block<3,1>(0,1) = T.block<3,1>(0,2).cross(pN - T.block<3,1>(0,3));
        Jacob.block<3,1>(3,1) = T.block<3,1>(0,2);
        T = T*T23;
        Jacob.block<3,1>(0,2) = T.block<3,1>(0,2).cross(pN - T.block<3,1>(0,3));
        Jacob.block<3,1>(3,2) = T.block<3,1>(0,2);
        T = T*T34;
        Jacob.block<3,1>(0,3) = T.block<3,1>(0,2).cross(pN - T.block<3,1>(0,3));
        Jacob.block<3,1>(3,3) =T.block<3,1>(0,2);
        T = T*T45;
        Jacob.block<3,1>(0,4) = T.block<3,1>(0,2).cross(pN - T.block<3,1>(0,3));
        Jacob.block<3,1>(3,4) = T.block<3,1>(0,2);
        T = T*T56;
        Jacob.block<3,1>(0,5) = T.block<3,1>(0,2).cross(pN - T.block<3,1>(0,3));
        Jacob.block<3,1>(3,5) = T.block<3,1>(0,2);
        return Jacob;
    }

/* Function: Determine whether the difference is within the allowable range
 * Inputs: the difference
 * Returns: the result(true means the joint angels are right and can be output)
 */
    bool Is_OK(const Eigen::VectorXd& diffe){
        for(int j = 0; j < 6; j++){
            if(fabs(diffe(j)) > 0.01){
                return false;
            }
        }
        return true;
    }
}

