#include <iostream>
#include <Eigen/Dense>
#include <math.h>
#include <string>
#include <fstream>

double numericIntegral(double delta_x, Eigen::VectorXd f, Eigen::VectorXd alpha)
{
    double result;
    for (int i = 0; i < alpha.size(); i++)
    {
        result += f(i) * alpha(i);
    }
    return delta_x * result;
}

double trapezRegel(Eigen::VectorXd integrand(Eigen::VectorXd), double a, double b, double h)
{
    double result = 0;
    double old_result = 0;
    int iterations = 0;

    do
    {
        h = h / 2.;
        old_result = result;
        result = 0;
        // Trapezregel
        Eigen::VectorXd alpha(2);
        alpha << 0.5, 0.5;
        Eigen::VectorXd x = Eigen::VectorXd::LinSpaced((b - a) / h + 1, a, b);
        Eigen::VectorXd f = integrand(x);

        for (int i = 0; i < f.size() - 1; i++)
        {
            result += numericIntegral(h, f(Eigen::seq(i, i + 1)), alpha);
        }
        iterations++;
    } while ((iterations < 10) & ((result - old_result) > 1e-1));

    std::cout << iterations << std::endl;
    return result;
}

Eigen::VectorXd quadraticFunction(Eigen::VectorXd x)
{
    return x.array() * x.array();
}

double integral(double a, double b)
{
    return 1. / 3. * ((b * b * b) - (a * a * a));
}

int main()
{
    std::ofstream output;
    output.open("output/result_ex1.txt", std::ofstream::trunc);

    output << trapezRegel(quadraticFunction, 1., 4., .5) << std::endl;

    output << integral(1., 4.) << std::endl;

    output.close();
    return 0;
}