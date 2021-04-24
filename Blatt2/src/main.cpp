#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <string>
#include <fstream>
#include <service.cpp>

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
    file.open("bin/data.txt", std::ofstream::out | std::ofstream::trunc);
	file << "# Picture with k:\n" << matA.format(CSVFormat) << std::endl; 
	file.close();

 // Compute U, W, V matrix 
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matA, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd W = svd.singularValues().asDiagonal();
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV().transpose();
    int range[3]={10,20,50};
	for(int i : range){
        Eigen::MatrixXd W_new = W(Eigen::seq(0,i), Eigen::seq(0,i));
        Eigen::MatrixXd U_new = U(Eigen::all, Eigen::seq(0,i));
        Eigen::MatrixXd V_new = V(Eigen::seq(0,i), Eigen::all);
        Eigen::MatrixXd mat_new = U_new * W_new* V_new;
        std::ofstream file;
        file.open("bin/data"+std::to_string(i)+".txt", std::ofstream::out | std::ofstream::trunc);
        file << "# Picture with k:"+std::to_string(i)+"\n" << mat_new.format(CSVFormat) << std::endl; 
        file.close();
    }
// 	int quality[3] = {10,20,50};
// 	for(int i : quality){
// 	Eigen::MatrixXd W_new = W(Eigen::seq(0,i), Eigen::seq(0,i));
//     Eigen::MatrixXd U_new = U(Eigen::all, Eigen::seq(0,i));
//     Eigen::MatrixXd V_new = V(Eigen::seq(0,i), Eigen::all);
//     Eigen::MatrixXd mat_new = U_new * W_new* V_new;
//     std::ofstream file;
//     file.open("bin/data"+std::to_string(i)+".txt", std::ofstream::out | std::ofstream::trunc);
// 	file << "# Picture with k:"+std::to_string(i)+"\n" << mat_new.format(CSVFormat) << std::endl; 
// 	file.close();
// 	}
 
	// std::ofstream output;
	// output.open("bin/Aufgabe1.txt", std::ofstream::out | std::ofstream::trunc);
	// output << "Matrix matA:\n" << matA.format(CSVFormat) << std::endl;
	// output.close();

    return 0;
}