#include <iostream>
#include <Eigen/Dense>
#include <string>
#include <cmath>
#include <fstream>
#include <numeric>

Eigen::ArrayXcd springchain(Eigen::VectorXd& m, Eigen::VectorXd& k, Eigen::VectorXd& l, int N){
    Eigen::VectorXd x(N);
    x(0) = 0.;
    std::partial_sum(l.begin(), l.end(), &x(1), std::plus<double>());

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
    // std::cout << K << std::endl;    
    return A.eigenvalues();
}

int main()
{
	// Ex. 1
    // 
    int N = 10;
    Eigen::VectorXd m(N), k(N-1), l(N-1);

    for (int i=0; i<N; i++){
        m(i) = i+1;
    }
    for (int j=0; j<N-1; j++){
        k(j) = N-j-1;
        l(j) = std::abs(5.-j-1) + 1;
    }

    Eigen::ArrayXcd omega_squared = springchain(m, k, l, N);
    Eigen::ArrayXd omega = omega_squared.real().abs().sqrt();
    // Position x
    std::cout << omega << std::endl;

    // std::ofstream file;
    // file.open("output/basis.txt", std::ofstream::trunc);
    // file << "\nDas Ergebnis" << result << std::endl;
    // file.close();

    return 0;
}