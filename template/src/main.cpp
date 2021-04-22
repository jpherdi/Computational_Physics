#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n");

int main()
{
	// Initialization of vectors and matrices
    Eigen::VectorXd v(3);
	Eigen::MatrixXd matA(3,3);

	v << 0, 1, 2;
	matA << v, v, v;
 
	std::ofstream output;
	output.open("bin/Aufgabe1.txt", std::ofstream::out | std::ofstream::trunc);
	output << "Matrix matA:\n" << matA.format(CSVFormat) << std::endl;
	output.close();

    return 0;
}