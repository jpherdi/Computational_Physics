#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <fstream>

const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n");

// int saveMatrix(Eigen::MatrixXd &mat, std::string filename, std::string header)
// {
// 	std::ofstream output;
// 	output.open(filename, std::fstream::trunc);
// 	if (output.is_open())
// 	{
// 		output << "# "+header << "\n" << mat.format(CSVFormat);
// 		output.close();
// 		return 1;
// 	}
// 	return 0;
// }

// int saveVector(Eigen::VectorXd &mat, std::string filename, std::string header)
// {
// 	std::ofstream output;
// 	output.open(filename, std::fstream::trunc);
// 	if (output.is_open())
// 	{
// 		output << "# "+header << "\n" << mat.format(CSVFormat);
// 		output.close();
// 		return 1;
// 	}
// 	return 0;
// }

int main()
{
	// Initialization of vectors and matrices
    Eigen::VectorXd a_1(3), a_2(3), a_3(3), x(3), y(3);
	Eigen::MatrixXd matA(3,3);

	// Definition of vectors and matrices
    a_1 << 1./2., sqrt(3.)/2., 0.;
    a_2 << -1./2., sqrt(3.)/2., 0.;
    a_3 << 0., 0., 1.;
    matA << a_1, a_2, a_3;

    // b), c) 
    x << 2., 0., 2.;
    y << 1., 2.*sqrt(3.), 3;

	// Use PartialPivLU as it was told in the exercise. 
	Eigen::PartialPivLU<Eigen::MatrixXd> matPivLU(matA);
	Eigen::VectorXd x_new = matPivLU.solve(x);
	Eigen::VectorXd y_new = matPivLU.solve(y);

	// Produce U, L, P matrices and check if everything works as expected. 
	Eigen::MatrixXd U = matPivLU.matrixLU().triangularView<Eigen::Upper>();
    Eigen::MatrixXd L = matPivLU.matrixLU().triangularView<Eigen::UnitLower>();
    Eigen::MatrixXd P = matPivLU.permutationP();
    Eigen::MatrixXd P_inv = P.inverse();
    Eigen::MatrixXd A_new = P_inv * L * U;

	// d) Do the same with the reversed (rev) basis vectors in the matrix.  
	Eigen::MatrixXd matA_rev(3,3);
    matA_rev << a_3, a_2, a_1;

	// Everything is done as before. 
	Eigen::PartialPivLU<Eigen::MatrixXd> matPivLU_rev(matA_rev);
	Eigen::VectorXd x_new_rev = matPivLU_rev.solve(x);
	Eigen::VectorXd y_new_rev = matPivLU_rev.solve(y);

	Eigen::MatrixXd U_rev = matPivLU_rev.matrixLU().triangularView<Eigen::Upper>();
    Eigen::MatrixXd L_rev = matPivLU_rev.matrixLU().triangularView<Eigen::UnitLower>();
    Eigen::MatrixXd P_rev = matPivLU_rev.permutationP();
    Eigen::MatrixXd P_inv_rev = P_rev.inverse();
    Eigen::MatrixXd A_new_rev = P_inv_rev * L_rev * U_rev;


	// Aufgabe 2 

	// Initialization of the vector and the matrix.
	Eigen::VectorXd x2(10), y2(10);
	Eigen::MatrixXd A(10, 2);
	
	// x2 and y2 are the input vectors. A is the matrix resulting from the original minimization problem.
	x2 <<  0., 2.5, -6.3, 4., -3.2, 5.3, 10.1, 9.5, -5.4, 12.7;
	A << x2, Eigen::VectorXd::Ones(10);
	y2 << 4., 4.3, -3.9, 6.5, 0.7, 8.6, 13., 9.9, -3.6, 15.1;

	// Transform A to a pseudo quadratic matrix as explained in the lecture and transform y2 as well. 
	Eigen::MatrixXd mat = A.transpose() * A;
	Eigen::VectorXd b = A.transpose() * y2;

	// Use the same calculation as in Exercise 1 to solve the problem.
	Eigen::PartialPivLU<Eigen::MatrixXd> matPivLU2(mat);
	Eigen::Vector2d alpha = matPivLU2.solve(b);

	Eigen::VectorXd result = A * alpha;
	// Save everything in a file to copy the results into an output file and to use it as input for a visualization script in python.

	Eigen::MatrixXd xy(10, 2);
	xy << x2, y2;
	std::ofstream file;
	file.open("bin/python_results.txt", std::ofstream::out | std::ofstream::trunc);
	file << "# x, y\n" << xy.format(CSVFormat) << std::endl; 
	file << "#m, n\n " << alpha.transpose().format(CSVFormat);
	file.close();


	// std::ofstream file2;
	// file2.open("bin/python_result.txt", std::ofstream::out | std::ofstream::trunc);
	
	// file2.close();
    return 0;
}