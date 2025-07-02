#include <iostream>
#include <fstream>
#include <cmath>

// dy/dt = log(abs(y)) - cos(t)
double dydt(double y, double t) {
    return std::log(std::abs(y)) - std::cos(t);
}

double next_step(double y, double h, double t) {
    double k_1 = h * dydt(y, t);
    double k_2 = h * dydt(y + k_1/2.0, t + h/2.0);
    double k_3 = h * dydt(y + k_2/2.0, t + h/2.0);
    double k_4 = h * dydt(y + k_3, t + h);

    return y + ((k_1 + 2.0 * k_2 + 2.0 * k_3 + k_4) / 6.0);
}

int main() {
    std::ofstream data;
    data.open("example1.txt");
    std::cout << "Hello" << std::endl;
    // initial values
    double y = 1;
    double h = 0.1;
    for (double t = 0; t <= 3; t = t + 0.01) {
        data << t << "," << y << std::endl;
        double dvt = dydt(y, t);
        y = next_step(y, h, dvt);
    }
    data.close();
}