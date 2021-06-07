#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

// Print Funktion zur Vereinfachung a lá python
template <typename T>
void print(T x) {
    std::cout << x << std::endl;
}

// Bisection Verfahren Analog zu Blatt 5 nur mit Vektoren
Eigen::Vector2d bisection(double function(Eigen::Vector2d), Eigen::Vector2d x0, Eigen::Vector2d gradient) {
    // print("Start bisection method.");
    double a, b, c;

    a = 0;
    b = 10;
    c = 100;
    int iterator = 0;
    while ((function(x0 + a * gradient) < function(x0 + b * gradient)) & (iterator < 100)) {
        a--;
        iterator++;
    }
    iterator = 0;
    // print("Find a. Iterations:");
    // print(iterator);
    while ((function(x0 + c * gradient) < function(x0 + b * gradient)) & (iterator < 100)) {
        c++;
        iterator++;
    }
    // print("Find c. Iterations:");
    // print(iterator);

    if ((a < b) & (b < c) & (function(x0 + b * gradient) < function(x0 + a * gradient)) & (function(x0 + b * gradient) < function(x0 + c * gradient))) {
        int iteration = 0;
        do {
            iteration++;
            if (abs(b - a) > abs(c - b)) {
                double a_new = a;
                double b_new = (a + b) / 2;
                double c_new = b;
                if (function(x0 + b_new * gradient) < function(x0 + c_new * gradient)) {
                    a = a_new;
                    b = b_new;
                    c = c_new;
                } else {
                    a = b_new;
                    b = c_new;
                    c = c;
                }
            } else {
                double a_new = b;
                double b_new = (b + c) / 2;
                double c_new = c;
                if (function(x0 + b_new * gradient) < function(x0 + a_new * gradient)) {
                    a = a_new;
                    b = b_new;
                    c = c_new;
                } else {
                    a = a;
                    b = a_new;
                    c = b_new;
                }
            }
        } while ((abs(a - c) > 1e-9) & (iteration < 100));
        // print("Iterations: ");
        // print(iteration);
        // print("a - c = ");
        // print(abs(a - c));
    } else {
        print("First requirement of bisection is not fullfiled.");
    }
    Eigen::Vector2d x_min = x0 + a * gradient;
    // print("x_min after Bisection.");
    // print(x_min);
    return x_min;
}

// BFGS Algorithmus
Eigen::Vector2d calc_BFGS(double function(Eigen::Vector2d), Eigen::Vector2d gradient(Eigen::Vector2d), Eigen::Vector2d x0, Eigen::Matrix2d C0, double epsilon, std::string filename) {
    // Definiere alle nötigen Objekte, x_k = x_k1, x_(k-1)= x_k0, b analog
    Eigen::Vector2d x_k0, x_k1, b_k0, b_k1, x_min, s, y;
    // Neue C Matrix, C1
    Eigen::Matrix2d C1;
    // Iterator um Schritte zu zählen
    double rho, iterator;
    iterator = 0;
    // Analog zu Skript. Bestimme b_0, mit bisections-Verfahren x1 und damit b1
    x_k0 = x0;
    b_k0 = gradient(x_k0);
    x_k1 = bisection(function, x_k0, b_k0);
    b_k1 = gradient(x_k1);
    // Speichere x_k in Datei
    std::ofstream output;
    output.open(filename, std::ofstream::trunc);  // std::ofstream::trunc);
    // Starte Iteration: Bestimmte s, y Vektoren und dann C1, anschließend nach Skript das neue x_k
    do {
        s = x_k1 - x_k0;
        y = b_k1 - b_k0;
        rho = 1. / (s.transpose() * y);
        C1 = (Eigen::Matrix2d::Identity() - rho * s * y.transpose()) * C0 * (Eigen::Matrix2d::Identity() - rho * y * s.transpose()) + rho * s * s.transpose();
        x_k0 = x_k1;
        C0 = C1;
        b_k0 = b_k1;

        x_k1 = x_k0 - C0 * b_k0;
        b_k1 = gradient(x_k1);
        iterator++;
        output << std::setprecision(8) << x_k1(0) << " " << x_k1(1) << std::endl;
        // Solange die Norm des Gradienten größer als epsilon ist oder bis eine gewisse Anzahl an Iterationen durchgelaufen ist.
    } while ((b_k1.norm() > epsilon) & (iterator < 1000));
    output.close();
    print("Iterations of BFGS:");
    print(iterator);
    x_min = x_k1;
    // Gib Minimum (1,1) zurück
    return x_min;
}

// 1. C0 Matrix Variante: Inverse Hesse-Matrix
Eigen::Matrix2d inv_HesseC_a(Eigen::Vector2d x) {
    Eigen::Matrix2d Hesse_C;
    Hesse_C(0, 0) = 2 - 400 * (x(1) - pow(x(0), 2.)) + 800 * pow(x(0), 2.);  // f_x1_x1
    Hesse_C(0, 1) = -400 * x(0);                                             // f_x1_x2
    Hesse_C(1, 0) = -400 * x(0);                                             // f_x2_x1
    Hesse_C(1, 1) = 200 * x(1);                                              // f_x2_x2
    return Hesse_C.inverse();
}

// 2. C0 Matrix Variante: inverse diag. Hesse-Matrix
Eigen::Matrix2d inv_HesseC_b(Eigen::Vector2d x) {
    Eigen::Matrix2d Hesse_C;
    Hesse_C(0, 0) = 2 - 400 * (x(1) - pow(x(0), 2.)) + 800 * pow(x(0), 2.);  // f_x1_x1
    Hesse_C(0, 1) = 0;                                                       // f_x1_x2
    Hesse_C(1, 0) = 0;                                                       // f_x2_x1
    Hesse_C(1, 1) = 200 * x(1);                                              //
    return Hesse_C.inverse();
}

// 3. C0 Matrix Variante: f(x_0) * identity
Eigen::Matrix2d inv_HesseC_c(Eigen::Vector2d x, double function(Eigen::Vector2d)) {
    Eigen::Matrix2d identity = Eigen::Matrix2d::Identity(2, 2);
    return function(x) * identity;
}

// Rosenbrock Funktion zu Aufgabe 1
double function_1(Eigen::Vector2d x) {
    return pow(1 - x(0), 2.) + 100 * pow(x(1) - pow(x(0), 2.), 2);
}

// Gradient zur Aufgabe 1
Eigen::Vector2d gradient_1(Eigen::Vector2d x) {
    Eigen::Vector2d g(2);
    g(0) = -2 * (1 - x(0)) - 400 * x(0) * (x(1) - pow(x(0), 2.));
    g(1) = 200 * (x(1) - pow(x(0), 2.));
    return g;
}

// double calc_RungeKutta(double y_prime(double, double), double t, double y_n, double h){
//     double k_1 = h* y_prime(t, y_n);
//     double k_2 = h* y_prime(t+ h/2., y_n + 0.5*k_1);
//     double k_3 = h* y_prime(t+h/2., y_n + 0.5*k_2);
//     double k_4 = h* y_prime(t+h, y_n + k_3);
//     double y_n_1 = y_n + 1./6. * (k_1 + 2* k_2 + 2*k_3 + k_4);
//     return y_n_1;
// }

int main() {
    // Nr. 1 Definiere Startvektor, Gradient und die C0 Matrix mit drei verschiedenen Varianten
    Eigen::Vector2d x(2);
    x << -1., -1.;
    Eigen::Vector2d grad = gradient_1(x);
    Eigen::Matrix2d C0_a = inv_HesseC_a(x);
    Eigen::Matrix2d C0_b = inv_HesseC_b(x);
    Eigen::Matrix2d C0_c = inv_HesseC_c(x, function_1);

    // Berechne die Ergebnisse mit dem oben definierten BFGS-Algorithmus für die drei verschiedenen Methoden
    print("1 a)");
    Eigen::Vector2d x_min_a = calc_BFGS(function_1, gradient_1, x, C0_a, 1e-5, "output/output_a.txt");
    print("Result of BFGS:");
    print(x_min_a);
    print("");

    print("1 b)");
    Eigen::Vector2d x_min_b = calc_BFGS(function_1, gradient_1, x, C0_b, 1e-5, "output/output_b.txt");
    print("Result of BFGS:");
    print(x_min_b);
    print("");

    print("1 c)");
    Eigen::Vector2d x_min_c = calc_BFGS(function_1, gradient_1, x, C0_c, 1e-5, "output/output_c.txt");
    print("Result of BFGS:");
    print(x_min_c);
    print("");

    return 0;
}