#include <cmath>
#include <fstream>
#include <iostream>

double dydt(double y, double t);
double next_step(double y, double h, double dydt);

int main() {
  std::ofstream data;
  data.open("example.txt");
  std::cout << "Hello" << std::endl;
  // initial values
  double y = 0;
  double h = 0.1;
  for (double t = 0; t <= 5.0; t = t + h) {
    data << t << "," << y << std::endl;
    double dvt = dydt(y, t);
    y = next_step(y, h, dvt);
  }
  data.close();
}

// differential equation: dy/dt = sin(y) - cos(t), iv of y(0) = 1

double dydt(double y, double t) {
  return std::pow(y-1, 2) * std::pow(t-1, 2);
}

// h is the step index
double next_step(double y, double h, double dydt) { return y + h * dydt; }