#include <cmath>
#include <functional>
#include <iostream>

double finite_difference(double y, std::function<double(double, double)> f,
                         double t, double h = 0.001) {
  return (f(y + h, t) - f(y, t)) / h;
}

double find_zero(std::function<double(double, double)> f, double t, double y,
                 double e = 1e-6) {
  std::cout << e << std::endl;
  double y_next;
  bool flag = true;
  int max_iter = 100;
  int iter = 0;
  while (flag && iter++ < max_iter) {
    std::cout << y << std::endl;
    y_next = y - f(y, t) / finite_difference(y, f, t);
    std::cout << y_next << std::endl;
    if ((std::abs(y_next - y)) < e) {
      flag = false;
      break;
    }
    y = y_next;
  }
  return y;
}

double test(double y, double t) { return std::pow(y, 2) - y; }

double BEuler_next_step(std::function<double(double, double)> f, double t,
                        double y, double h) {
  auto G = [=](double t, double y_next) {
    return y_next - y - h * f(y_next, t + h);
  };
  return find_zero(G, t + h, y);
}
int main() { std::cout << find_zero(test, 4, 4) << std::endl; }
