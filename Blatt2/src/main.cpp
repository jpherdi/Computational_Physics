#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <string>
#include <fstream>
#include <service.cpp>
#include <profiler.h>

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n");

int main()
{
	// Ex. 1
    // Load and store data
    Eigen::MatrixXd matA(512,512);
    std::string filename = "data/Bild";
    size_t row = 512;
    size_t col = 512;
    int load = loadData(matA, filename, row, col);
    std::ofstream file;
    file.open("output/data.txt", std::ofstream::out | std::ofstream::trunc);
	file << "# Picture with k:\n" << matA.format(CSVFormat) << std::endl; 
	file.close();

    // Compute U, W, V matrix 
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matA, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd W = svd.singularValues().asDiagonal();
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV().transpose();

    // Perform k-approx with values in range array
    int range[3]={10,20,50};
	for(int i : range){
        // Use Eigen::seq to decrease dimension of W, and of U and V respectively
        Eigen::MatrixXd W_new = W(Eigen::seq(0,i), Eigen::seq(0,i));
        Eigen::MatrixXd U_new = U(Eigen::all, Eigen::seq(0,i));
        Eigen::MatrixXd V_new = V(Eigen::seq(0,i), Eigen::all);
        Eigen::MatrixXd mat_new = U_new * W_new* V_new;
        // Save new picture in txt file, to make a simple read out in python
        std::ofstream file;
        file.open("output/data"+std::to_string(i)+".txt", std::ofstream::out | std::ofstream::trunc);
        file << "# Picture with k:"+std::to_string(i)+"\n" << mat_new.format(CSVFormat) << std::endl; 
        file.close();
    }

    // Ex. 2

    // Number of dimension N linear from 1 to 1000
    Profiler::init(3);
    int anzahl = 1000;
    int start = 1;

    // Vectors for the building of the random NxN matrix, for the LU decomposition and the solve function
    Eigen::VectorXd random_times(anzahl - start);
    Eigen::VectorXd LU_times(anzahl - start);
    Eigen::VectorXd solve_times(anzahl - start);

    for (int N = start; N < anzahl; N++)
    {
        //All timers are reseted
        Profiler::resetAll();
        //Timer for random initialisation and saving
        Profiler::start(0);
        Eigen::MatrixXd M = Eigen::MatrixXd::Random(N, N);
          Profiler::stop(0);
        random_times[N - start] = Profiler::getTimeInS(0);
        Eigen::VectorXd b = Eigen::VectorXd::Random(N);
        //Timer for LU decomposition without saving
        Profiler::start(1);
        M.partialPivLu();
        Profiler::stop(1);
        LU_times[N - start] = Profiler::getTimeInS(1);
        Eigen::PartialPivLU<Eigen::MatrixXd> mat(M);
        //Timer for solve function without saving
        Profiler::start(2);
        Eigen::VectorXd x = mat.solve(b);
        Profiler::stop(2);
        solve_times[N - start] = Profiler::getTimeInS(2);
    }
    // Store vectors in file
    Eigen::MatrixXd store(anzahl-start, 3);
    store << random_times, LU_times, solve_times;
    std::ofstream file2;
    file2.open("output/times.txt", std::ofstream::out | std::ofstream::trunc);
	file2 << "# Times: Random, LU, Solve\n" << store.format(CSVFormat) << std::endl; 
	file2.close();

    //Ex. 3
    // Initialize Matrix A with a1 - a10 as vector space
    Eigen::VectorXd a1(4), a2(4), a3(4), a4(4), a5(4), a6(4), a7(4), a8(4), a9(4), a10(4);
    a1 << 4., 1., 2., 4.;
    a2 << 1., 2., 6., 2.;
    a3 << 2., 6., 9., 8.;
    a4 << 5., 3., 8., 6.;
    a5 << 8., 2., 4., 8.;
    a6 << 6., 7., 11., 12.;
    a7 << 3., 8., 15., 10.; 
    a8 << -2., 5., 7., 4.;
    a9 << 12,  3., 6., 12.;
    a10 << 3., -1., -4., 2.; 

    Eigen::MatrixXd A(4, 10);
    A << a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;

    // find orthonormal basis of a1-a10 with singular value decomposition
    Eigen::JacobiSVD<Eigen::MatrixXd> svd2(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd W2 = svd2.singularValues().asDiagonal();
    Eigen::MatrixXd U2 = svd2.matrixU();
    Eigen::MatrixXd V2 = svd2.matrixV().transpose();

    std::ofstream file3;
    file3.open("output/basis.txt", std::ofstream::trunc);
    // basis vectors are the columns of U with w_i != 0
    for (int i = 0; i< 4; i++){
        if(W2(i,i) >= 1e-10){
            file3 << "Value W_" << i << " is not zero. The corresponding vector is\n" << U2.col(i).format(CSVFormat) << "\n\n";
        }
    }
    file3 << "\nDie U-Matrix" << U2.format(CSVFormat) << std::endl;
    file3 << "\nDie W-Matrix" << W2.format(CSVFormat) << std::endl;
    file3.close();

    return 0;
}