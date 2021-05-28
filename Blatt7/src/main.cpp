#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <fstream>

template <typename T>
void print(T x)
{
    std::ofstream output;
    output.open("output/output.txt", std::ios_base::app); // std::ofstream::trunc);
    std::cout << x << std::endl;
    output.close();
}

int main()
{

    return 0;
}