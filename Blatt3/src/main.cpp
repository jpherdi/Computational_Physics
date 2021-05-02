#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <cmath>
#include <fstream>
#include <numeric>

// Eigen::ArrayXcd springchain(Eigen::VectorXd& m, Eigen::VectorXd& k, Eigen::VectorXd& l, int N){
//     Eigen::VectorXd x(N);
//     x(0) = 0.;
//     std::partial_sum(l.begin(), l.end(), &x(1), std::plus<double>());

//     Eigen::MatrixXd M_inv = m.asDiagonal().inverse();
//     Eigen::MatrixXd K(N, N);
//     for (int i= 0; i < N; i++){
//         for(int j= 0; j <= i; j++){
//             K(i,j) = 0;
//             if(i==j){
//                 if(i == 0){
//                     K(i,i) = k(i);
//                     K(i, i+1) = -k(i);
//                 }
//                 else if(i== N-1){
//                     K(i,i) = k(i-1);
//                     K(i, i-1) = -k(i-1);
//                 }
//                 else{
//                     K(i,i) = k(i-1)+k(i);
//                     K(i, i-1) = -k(i-1);
//                     K(i, i+1) = -k(i);
//                 }
//             }
//         }
//     }
//     Eigen::MatrixXd A = M_inv*K;
//     // std::cout << K << std::endl;    
//     return A.eigenvalues();
// }

void lanczos(Eigen::MatrixXd& A){

    Eigen::MatrixXd q(A.rows(), A.cols()+2);
    q(0,0) = 0;
    q(0,1) = 1;
    Eigen::VectorXd gamma(A.cols()+2);
    gamma(1) = 1;
    Eigen::VectorXd delta(A.cols()+1);
    
    for (int i=1; i<=A.rows(); i++){
        gamma(i+1) = std::sqrt(q.col(i).transpose()* q.col(i));
        delta(i) = q.col(i).transpose()*A*q.col(i);
        q.col(i+1) = 1/gamma(i+1) * (A- delta(1)*Eigen::MatrixXd::Identity(A.rows(), A.cols()))*q.col(i) - gamma(i)*q.col(i-1);
    }

    Eigen::MatrixXd T(A.rows(), A.cols());
    for (int i= 0; i < A.rows()-1; i++){
        T(i,i) = delta(i+1);
        T(i, i+1) = gamma(i+3);
        T(i+1, i) = gamma(i+3);
    }
    T(A.rows()-1, A.cols()-1) = delta(A.cols());
    Eigen::MatrixXd T_old = T;
    //std::cout << T << std::endl; 
    Eigen::MatrixXd Z = Eigen::MatrixXd::Identity(T.rows(), T.cols());

    for (int i=0; i < T.rows()-3; i++){
        Eigen::JacobiRotation<double> J;
        Eigen::MatrixXd J_1 = Eigen::MatrixXd::Identity(T.rows(), T.cols());
        J.makeJacobi(T(Eigen::seq(i,i+1), Eigen::seq(i,i+1)), 0, 1);
        J_1(Eigen::seq(i,i+1), Eigen::seq(i,i+1)) << J.c(), J.s(), -J.s(), J.c();
        T(Eigen::seq(i,i+1), Eigen::seq(i,i+1)).applyOnTheLeft(0, 1, J.adjoint());
        T(Eigen::seq(i,i+1), Eigen::seq(i,i+1)).applyOnTheRight(0, 1, J);
        std::cout << J_1 << std::endl;

        // J.makeJacobi(T(Eigen::seq(i,i+1), Eigen::seq(i,i+1)), 1, 0);
        // Eigen::MatrixXd J_2 = Eigen::MatrixXd::Identity(T.rows(), T.cols());
        // J_2(Eigen::seq(i,i+1), Eigen::seq(i,i+1)) << J.c(), J.s(), -J.s(), J.c();
        // T(Eigen::seq(i,i+1), Eigen::seq(i,i+1)).applyOnTheLeft(1, 0, J.adjoint());
        // T(Eigen::seq(i,i+1), Eigen::seq(i,i+1)).applyOnTheRight(1, 0, J);
        Z = Z * J_1;// *J_2;
    }

    std::cout << T << std::endl;
    // std::cout << T_old << std::endl;
    std::cout << Z.transpose() * T_old * Z << std::endl;

}

int main()
{
	// Ex. 1
    // 
    // int N = 10;
    // Eigen::VectorXd m(N), k(N-1), l(N-1);

    // for (int i=0; i<N; i++){
    //     m(i) = i+1;
    // }
    // for (int j=0; j<N-1; j++){
    //     k(j) = N-j-1;
    //     l(j) = std::abs(5.-j-1) + 1;
    // }

    // Eigen::ArrayXcd omega_squared = springchain(m, k, l, N);
    // Eigen::ArrayXd omega = omega_squared.real().abs().sqrt();
    // // Position x
    // std::cout << omega << std::endl;

    // Ex. 2
    int N = 4;
    Eigen::MatrixXd A= Eigen::MatrixXd::Random(N,N);
    lanczos(A);

    // std::ofstream file;
    // file.open("output/basis.txt", std::ofstream::trunc);
    // file << "\nDas Ergebnis" << result << std::endl;
    // file.close();

    return 0;
}