#include <cmath>
#include <fstream>
#include <iostream>

// dy/dt = log(abs(y)) - cos(t)
double dydt(double y, double t) { return std::pow(y-1, 2) * std::pow(t-1, 2);}

double next_step(double y, double h, double t) {
  double k_1 = h * dydt(y, t);
  double k_2 = h * dydt(y + k_1 / 2.0, t + h / 2.0);
  double k_3 = h * dydt(y + k_2 / 2.0, t + h / 2.0);
  double k_4 = h * dydt(y + k_3, t + h);

  return y + ((k_1 + 2.0 * k_2 + 2.0 * k_3 + k_4) / 6.0);
}

int main() {
  std::ofstream data;
  data.open("example1.txt");
  std::cout << "Hello" << std::endl;
  // initial values
  double y = 0;
  double h = 0.01;
  for (double t = 0; t <= 5.0; t = t + h) {
    data << t << "," << y << std::endl;
    y = next_step(y, h, t);
  }
  data.close();
}