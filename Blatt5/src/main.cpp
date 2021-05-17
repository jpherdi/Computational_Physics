#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <fstream>
#include <bits/stdc++.h>

template <typename T>
void print(T x)
{
    std::ofstream output;
    output.open("output/output.txt", std::ios_base::app); // std::ofstream::trunc);
    std::cout << x << std::endl;
    output.close();
}

double trapezRule(Eigen::VectorXd function(Eigen::VectorXd), double x_min, double x_max, double h)
{
    double result = 100;
    double old_result = 0;
    int iterator = 0;

    do
    {
        h = h / 2;
        iterator++;
        Eigen::VectorXd deltaX = Eigen::VectorXd::LinSpaced((x_max - x_min) / h + 1, x_min, x_max);
        Eigen::VectorXd f = function(deltaX);
        Eigen::VectorXd result_vec = f * h;
        old_result = result;
        result = 0.5 * (result_vec(0) + result_vec(result_vec.size() - 1)) + result_vec(Eigen::seq(1, result_vec.size() - 2)).sum();

    } while ((abs(result - old_result) > 1e-6) & (iterator <= 30));
    print("Iterationen:");
    print(iterator);
    return result;
}

Eigen::VectorXd f1a(Eigen::VectorXd x)
{
    return Eigen::exp(x.array()) / x.array();
}

Eigen::VectorXd f1b(Eigen::VectorXd x)
{
    return Eigen::exp(-x.array()) / Eigen::sqrt(x.array());
}

Eigen::VectorXd f1c(Eigen::VectorXd x)
{
    return Eigen::sin(x.array()) / x.array(); //sin(x) / x;
}

double f3a(double x)
{
    return (x * x) - 2;
}

double bisection(double function(double), double a, double b, double c)
{
    if ((a < b) & (b < c) & (function(b) < function(a)) & (function(b) < function(c)))
    {
        int iteration = 0;
        do
        {
            iteration++;

            if (abs(b - a) > abs(c - b))
            {
                double a_new = a;
                double b_new = (a + b) / 2;
                double c_new = b;
                if (function(b_new) < function(c_new))
                {
                    a = a_new;
                    b = b_new;
                    c = c_new;
                }
                else
                {
                    a = b_new;
                    b = c_new;
                    c = c;
                }
            }
            else
            {
                double a_new = b;
                double b_new = (b + c) / 2;
                double c_new = c;
                if (function(b_new) < function(a_new))
                {
                    a = a_new;
                    b = b_new;
                    c = c_new;
                }
                else
                {
                    a = a;
                    b = a_new;
                    c = b_new;
                }
            }
        } while ((abs(a - c) > 1e-9) & (iteration < 100));
        print("Iterationen: ");
        print(iteration);
        print("a - c = ");
        print(abs(a - c));
    }
    else
    {
        print("First requirement of bisection is not fullfiled.");
    }
    return a;
}

double fDot(double function(double), double x0)
{
    double deltaX = (x0 + 1) / 1000;
    // Eigen::VectorXd deltaX = Eigen::VectorXd::LinSpaced((x(x.size()-1) -x(0) ) / diff_x + 1, x(0), x(x.size()-1));
    return (function(x0 + deltaX) - function(x0 - deltaX)) / (2 * deltaX);
}

double fDDot(double function(double), double x0)
{
    double deltaX = (x0 + 1) / 1000;
    // Eigen::VectorXd deltaX = Eigen::VectorXd::LinSpaced((x(x.size()-1) -x(0) ) / diff_x + 1, x(0), x(x.size()-1));
    return (function(x0 + deltaX) - 2 * function(x0) + function(x0 - deltaX)) / (deltaX * deltaX);
}

double newton(double function(double), double x0)
{
    double x_new = x0;
    double x_old = 0;
    double iteration = 0;
    do
    {
        iteration++;
        x_old = x_new;
        x_new = x_old - fDot(function, x_old) / fDDot(function, x_old);
    } while ((abs(x_new - x_old) > 1e-9) & (iteration < 100));
    print("Iterationen: ");
    print(iteration);
    print("x_new - x_old = ");
    print(abs(x_new - x_old));
    return x_new;
}

int main()
{
    for (int m = 3; m < 5; m++)
    {
        int j = pow(2, m);
        int N = 10;
        int l = N;
        std::complex<double> complex_(0.0, 1.0);
        Eigen::VectorXcd F_j(j);
        Eigen::VectorXcd f_l(l);
        Eigen::MatrixXcd omega(j, l);

        for (int i = 0; i < l; i++)
        {
            f_l(i) = std::sqrt(1 + i);
        }

        for (int i = 0; i < l; i++)
        {
            for (int k = 0; k < j; k++)
            {
                omega(k, i) = std::exp(k * i / N * 2 * M_PI * complex_);
            }
        }

        F_j = omega * f_l;
        Eigen::VectorXd F_j_real = F_j.real();
        print(F_j_real);
    }

    // // Ex. 1
    // // a)

    // print("Exercise a): ");
    // double zero = 1e-6;
    // print("Zero is:");
    // print(zero);
    // double i1 = trapezRule(f1a, -1, -zero, 0.5) + trapezRule(f1a, zero, 1, 0.5);
    // print("Result I1: ");
    // print(i1);

    // // b) see WolframAlpha for x_max
    // print("Exercise b): ");
    // double infty = 100.;
    // print("Infinity is:");
    // print(infty);
    // double epsilon = 0.001;
    // print("Epsilon is:");
    // print(epsilon);
    // zero = 1e-11;
    // print("Zero is:");
    // print(zero);
    // double i2 = trapezRule(f1b, zero, epsilon, 0.0005) + trapezRule(f1b, epsilon, infty, 0.5);
    // print("Result I2: ");
    // print(i2);

    // // c)
    // print("Exercise c): ");
    // infty = 1000000.;
    // print("Infinity is:");
    // print(infty);
    // zero = 1e-6;
    // print("Zero is:");
    // print(zero);
    // double i3 = trapezRule(f1c, -infty, -zero, 0.5) + trapezRule(f1c, zero, infty, 0.5);
    // print("Result I3: ");
    // print(i3);

    // // Ex. 3
    // // a)
    // double a = -0.5;
    // double b = -0.1;
    // double c = 2.0;
    // double result_bisec = bisection(f3a, a, b, c);

    // // b)
    // double x0 = 1;
    // double result_newton = newton(f3a, x0);
    // print("Bisec - Newton:");
    // print(abs(result_bisec - result_newton));

    return 0;
}