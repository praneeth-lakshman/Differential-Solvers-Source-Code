#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

class Values {
 public:
  double y;
  double h;
  bool accepted;
};

const double a11 = 0;
const double a21 = 1.0 / 5;
const double a31 = 3.0 / 40, a32 = 9.0 / 40;
const double a41 = 44.0 / 45, a42 = -56.0 / 15, a43 = 32.0 / 9;
const double a51 = 19372.0 / 6561, a52 = -25360.0 / 2187, a53 = 64448.0 / 6561,
             a54 = -212.0 / 729;
const double a61 = 9017.0 / 3168, a62 = -355.0 / 33, a63 = 46732.0 / 5247,
             a64 = 49.0 / 176, a65 = -5103.0 / 18656;
const double a71 = 35.0 / 384, a72 = 0.0, a73 = 500.0 / 1113, a74 = 125.0 / 192,
             a75 = -2187.0 / 6784, a76 = 11.0 / 84;

const double c1 = 0;
const double c2 = 1.0 / 5;
const double c3 = 3.0 / 10;
const double c4 = 4.0 / 5;
const double c5 = 8.0 / 9;
const double c6 = 1.0;
const double c7 = 1.0;

const double b11 = 35.0 / 184;
const double b12 = 0;
const double b13 = 500.0 / 1113;
const double b14 = 125.0 / 192;
const double b15 = -2187.0 / 6784;
const double b16 = 11.0 / 84;
const double b17 = 0;
const double b21 = 5179.0 / 57600;
const double b22 = 0.0;
const double b23 = 7571.0 / 16695;
const double b24 = 393.0 / 640;
const double b25 = -92097.0 / 339200;
const double b26 = 187.0 / 2100;
const double b27 = 1.0 / 40;

const double eps = 0.00005;

double f(double y, double t) { return std::pow(y - 1, 2) * std::pow(t - 1, 2); }

Values next_step(double y, double h, double t) {
  double k1 = h * f(y, t + h * c1);
  double k2 = h * f(y + a21 * k1, t + h * c2);
  double k3 = h * f(y + a31 * k1 + a32 * k2, t + h * c3);
  double k4 = h * f(y + a41 * k1 + a42 * k2 + a43 * k3, t + h * c4);
  double k5 = h * f(y + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4, t + h * c5);
  double k6 = h * f(y + a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5,
                    t + h * c6);
  double k7 =
      h * f(y + a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6,
            t + h * c7);

  double rk4_y = y + b21 * k1 + b22 * k2 + b23 * k3 + b24 * k4 + b25 * k5 +
                 b26 * k6 + b27 * k7;
  double rk5_y = y + b11 * k1 + b12 * k2 + b13 * k3 + b14 * k4 + b15 * k5 +
                 b16 * k6 + b17 * k7;
  double e = std::abs(rk4_y - rk5_y);
  double factor = std::pow(eps / (e + 1e-10), 0.25);
  double h_new = h * std::clamp(0.84 * factor, 0.1, 5.0);
  Values result;
  result.h = h_new;
  result.accepted = (e <= eps);
  if (result.accepted) {
    result.y = rk5_y;
  } else {
    result.y = y;
  }
  return result;
}

int main() {
  std::ofstream data("example3.txt");
  std::cout << "Starting adaptive DP Method..." << std::endl;

  double y = 0;
  double h = 0.1;
  double t = 0.0;

  while (t <= 5.0) {
    Values vals = next_step(y, h, t);

    if (vals.h < 1e-6) {
      std::cerr << "Step size too small at t = " << t << std::endl;
      break;
    }

    if (vals.accepted) {
      data << t << "," << y << std::endl;
      y = vals.y;
      t += h;
    }

    h = vals.h;
  }

  data.close();
  std::cout << "Integration complete. Data saved to example3.txt\n";
  return 0;
}