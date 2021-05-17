#include <iostream>
#include <cmath>
#include <vector>

//Integrand i)
double f(double x)
{
    return exp(-x) / x;
}

//Integrad ii)
double g(double x)
{
    if (x == 0)
    {
        return 1;
    }
    return x * sin(1 / x);
}

//Trapezregel i)
double trapez(double a, double b, int n)
{
    int i;
    double h, s;

    h = (b - a) / n;
    s = 0.5 * (f(a) + f(b));
    for (i = 1; i <= n - 1; i++)
        s = s + f(a + i * h);

    return s * h;
}

//Simpsonregel i)
double simpson(double a, double b, int n)
{
    int i;
    double h, s;

    h = (b - a) / n;
    s = f(a) + f(b) + 4 * f(a + h / 2);
    for (i = 1; i <= n - 1; i++)
        s = s + 2 * f(a + i * h) + 4 * f(a + i * h + h / 2);

    return s * h / 6;
}

//Mittelpunktsregel i)
double mittelpunkt(double a, double b, int n)
{
    int i;
    double h, s;

    h = (b - a) / n;
    for (i = 1; i <= n; i++)
        s = s + f(a + i * h + h / 2);

    return h * s;
}

//Formeln noch mal für ii) (also mit g statt f)

//Trapezregel ii)
double trapez_ii(double a, double b, int n)
{
    int i;
    double h, s;

    h = (b - a) / n;
    s = 0.5 * (g(a) + g(b));
    for (i = 1; i <= n - 1; i++)
        s = s + g(a + i * h);

    return s * h;
}

//Simpsonregel ii)
double simpson_ii(double a, double b, int n)
{
    int i;
    double h, s;

    h = (b - a) / n;
    s = g(a) + g(b) + 4 * g(a + h / 2);
    for (i = 1; i <= n - 1; i++)
        s = s + 2 * g(a + i * h) + 4 * g(a + i * h + h / 2);

    return s * h / 6;
}

//Mittelpunktsregel ii)
double mittelpunkt_ii(double a, double b, int n)
{
    int i;
    double h, s;

    h = (b - a) / n;
    for (i = 1; i <= n; i++)
        s = s + g(a + i * h + h / 2);

    return h * s;
}

int main()
{
    //Integrationsgrenzen
    int a1 = 1;
    int b1 = 100;
    int a2 = 0;
    int b2 = 1;

    int n = 99; //Anzahl Intervalle

    int i = 0; //Counter

    //Ergebnisvektor
    std::vector<double> ergebnisse;

    ergebnisse.push_back(1.0);
    ergebnisse.push_back(2.0);

    //i)

    //a) Trapezregel für i)
    while (abs(ergebnisse.back() - ergebnisse[ergebnisse.size() - 2]) > pow(10, -4))
    {
        ergebnisse.push_back(trapez(a1, b1, n));
        i++;
        n = 2 * n;
    }
    std::cout << "Trapezregel für i)" << std::endl;
    std::cout << "Counter = " << i << std::endl;
    std::cout << "n = " << n << std::endl; //daraus lässt sich h bestimmen

    //zurücksetzen
    i = 0;
    n = 99;
    ergebnisse.clear();
    ergebnisse.push_back(1.0);
    ergebnisse.push_back(2.0);

    //b) Mittelpunktsregel für i)
    while (abs(ergebnisse.back() - ergebnisse[ergebnisse.size() - 2]) > pow(10, -4))
    {
        ergebnisse.push_back(mittelpunkt(a1, b1, n));
        i++;
        n = 2 * n;
    }
    std::cout << "Mittelpunktsregel für i)" << std::endl;
    std::cout << "Counter = " << i << std::endl;
    std::cout << "n = " << n << std::endl; //daraus lässt sich h bestimmen

    //zurücksetzen
    i = 0;
    n = 99;
    ergebnisse.clear();
    ergebnisse.push_back(1.0);
    ergebnisse.push_back(2.0);

    //c) Simpsonregel für i)
    while (abs(ergebnisse.back() - ergebnisse[ergebnisse.size() - 2]) > pow(10, -4))
    {
        ergebnisse.push_back(simpson(a1, b1, n));
        i++;
        n = 2 * n;
    }
    std::cout << "Simpsonregel für i)" << std::endl;
    std::cout << "Counter = " << i << std::endl;
    std::cout << "n = " << n << std::endl; //daraus lässt sich h bestimmen

    //zurücksetzen
    i = 0;
    n = 99;
    ergebnisse.clear();
    ergebnisse.push_back(1.0);
    ergebnisse.push_back(2.0);

    //ii)
    n = 10;

    //a)Trapezregel für ii)
    while (abs(ergebnisse.back() - ergebnisse[ergebnisse.size() - 2]) > pow(10, -4))
    {
        ergebnisse.push_back(trapez_ii(a2, b2, n));
        i++;
        n = 2 * n;
    }
    std::cout << "Trapezregel für ii)" << std::endl;
    std::cout << "Counter = " << i << std::endl;
    std::cout << "n = " << n << std::endl; //daraus lässt sich h bestimmen

    //zurücksetzen
    i = 0;
    n = 10;
    ergebnisse.clear();
    ergebnisse.push_back(1.0);
    ergebnisse.push_back(2.0);

    //b) Mittelpunktsregel für ii)
    while (abs(ergebnisse.back() - ergebnisse[ergebnisse.size() - 2]) > pow(10, -4))
    {
        ergebnisse.push_back(mittelpunkt_ii(a2, b2, n));
        i++;
        n = 2 * n;
    }
    std::cout << "Mittelpunktsregel für ii)" << std::endl;
    std::cout << "Counter = " << i << std::endl;
    std::cout << "n = " << n << std::endl; //daraus lässt sich h bestimmen

    //zurücksetzen
    i = 0;
    n = 10;
    ergebnisse.clear();
    ergebnisse.push_back(1.0);
    ergebnisse.push_back(2.0);

    //c) Simpsonregel für ii)
    while (abs(ergebnisse.back() - ergebnisse[ergebnisse.size() - 2]) > pow(10, -4))
    {
        ergebnisse.push_back(simpson_ii(a2, b2, n));
        i++;
        n = 2 * n;
    }
    std::cout << "Simpsonregel für ii)" << std::endl;
    std::cout << "Counter = " << i << std::endl;
    std::cout << "n = " << n << std::endl; //daraus lässt sich h bestimmen

    //zurücksetzen
    i = 0;
    n = 10;
    ergebnisse.clear();
    ergebnisse.push_back(1.0);
    ergebnisse.push_back(2.0);

    return 0;
}