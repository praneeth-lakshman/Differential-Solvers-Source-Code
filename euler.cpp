#include <iostream>
#include <fstream>
#include <cmath>

double dydt(double y, double t);
double next_step(double y, double h, double dydt);

int main() {
    std::ofstream data;
    data.open("example.txt");
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

// differential equation: dy/dt = sin(y) - cos(t), iv of y(0) = 1

double dydt(double y, double t) {
    return std::log(std::abs(y)) - std::cos(t);
}

// h is the step index
double next_step(double y, double h, double dydt) {
    return y + h * dydt;
}