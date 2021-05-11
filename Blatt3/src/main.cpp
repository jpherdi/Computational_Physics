#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <cmath>
#include <fstream>
#include <numeric>

// Input: m Vector, k Vector, l Vector (even though it is not needed as we later realized), dimension N
Eigen::ArrayXcd springchain(Eigen::VectorXd& m, Eigen::VectorXd& k, Eigen::VectorXd& l, int N){
    // Initialisation of relative coordinates x (in protocol eta)
    Eigen::VectorXd x(N);
    x(0) = 0.;
    std::partial_sum(l.begin(), l.end(), &x(1), std::plus<double>());

    // Definition of M^-1 and K matrix
    Eigen::MatrixXd M_inv = m.asDiagonal().inverse();
    Eigen::MatrixXd K(N, N);
    for (int i= 0; i < N; i++){
        for(int j= 0; j <= i; j++){
            K(i,j) = 0;
            if(i==j){
                if(i == 0){
                    K(i,i) = k(i);
                    K(i, i+1) = -k(i);
                }
                else if(i== N-1){
                    K(i,i) = k(i-1);
                    K(i, i-1) = -k(i-1);
                }
                else{
                    K(i,i) = k(i-1)+k(i);
                    K(i, i-1) = -k(i-1);
                    K(i, i+1) = -k(i);
                }
            }
        }
    }
    Eigen::MatrixXd A = M_inv*K;
    // Calculate Eigenvalues   
    return A.eigenvalues();
}

// Lanczos algorithm, sadly not finished. Returns only the tridiag matrix
Eigen::MatrixXd lanczos(Eigen::MatrixXd& A){

    // Implementation as defined in the lecture

    // q matrix with the new basis q
    Eigen::MatrixXd q(A.rows(), A.cols()+2);
    q(0,0) = 0;
    q(0,1) = 1;
    // gamma vector with all the normalisations
    Eigen::VectorXd gamma(A.cols()+2);
    gamma(1) = 1;
    // delta vector with diagonal elements of the new tridiag matrix
    Eigen::VectorXd delta(A.cols()+1);
    
    // the new basis build with a for loop
    for (int i=1; i<=A.rows(); i++){
        gamma(i+1) = std::sqrt(q.col(i).transpose()* q.col(i));
        delta(i) = q.col(i).transpose()*A*q.col(i);
        q.col(i+1) = 1/gamma(i+1) * (A- delta(1)*Eigen::MatrixXd::Identity(A.rows(), A.cols()))*q.col(i) - gamma(i)*q.col(i-1);
    }

    // make new tridiag matrix based on gamma-, and delta-vector, with some error
    Eigen::MatrixXd T(A.rows(), A.cols());
    for (int i= 0; i < A.rows()-1; i++){
        T(i,i) = delta(i+1);
        T(i, i+1) = gamma(i+3);
        T(i+1, i) = gamma(i+3);
    }
    T(A.rows()-1, A.cols()-1) = delta(A.cols());

    return T;
    
    // Begin of Jacobi-Rotation Implementation, sadly not finished
    Eigen::MatrixXd Z = Eigen::MatrixXd::Identity(T.rows(), T.cols());

    for (int i=0; i < T.rows()-1; i++){
        Eigen::JacobiRotation<double> J;
        Eigen::MatrixXd J_1 = Eigen::MatrixXd::Identity(T.rows(), T.cols());
        J.makeJacobi(T(Eigen::seq(i,i+1), Eigen::seq(i,i+1)), 0, 1);
        J_1(Eigen::seq(i,i+1), Eigen::seq(i,i+1)) << J.c(), J.s(), -J.s(), J.c();
        T(Eigen::seq(i,i+1), Eigen::seq(i,i+1)).applyOnTheLeft(0, 1, J.adjoint());
        T(Eigen::seq(i,i+1), Eigen::seq(i,i+1)).applyOnTheRight(0, 1, J);
        Z = Z * J_1;
    }

}

int main()
{
	// Ex. 1
    // N is the dimension, here 10 for part b)
    int N = 10;
    Eigen::VectorXd m(N), k(N-1), l(N-1);
    // Initialisation of m, k and l
    for (int i=0; i<N; i++){
        m(i) = i+1;
    }
    for (int j=0; j<N-1; j++){
        k(j) = N-j-1;
        l(j) = std::abs(5.-j-1) + 1;
    }
    // Calculate omega_squared with springchain function
    Eigen::ArrayXcd omega_squared = springchain(m, k, l, N);
    // Calculate omega with taking the square root of all omega_squared values
    Eigen::ArrayXd omega = omega_squared.real().abs().sqrt();
    
    // Output in file ./output/omega.txt
    std::ofstream file;
    file.open("output/omega.txt", std::ofstream::trunc);
    file << "\nDas Ergebnis:\n" << omega << std::endl;
    file.close();

    // Ex. 2
    int dim = 5;
    Eigen::MatrixXd A= Eigen::MatrixXd::Random(dim,dim);
    Eigen::MatrixXd T = lanczos(A);

    std::ofstream file2;
    file2.open("output/tridiag.txt", std::ofstream::trunc);
    file2 << "Die Zufallsmatrix A:\n" << A << std::endl;
    file2 << "\nwird zur Tridiagonalmatrix:\n" << T << std::endl;
    file2.close();

    return 0;
}


