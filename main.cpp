#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <pybind11/pybind11.h>
#include <functional>

namespace py = pybind11;

class Euler {
    public:
    Euler() = default;
    double next_step(double y, double h, py::function f, double t) { return y + h * f(y, t).cast<double>(); }
};

class RK4 {
    public:
    double next_step(double y, double h, py::function f, double t) {
    double k_1 = h * f(y, t).cast<double>();
    double k_2 = h * f(y + k_1 / 2.0, t + h / 2.0).cast<double>();
    double k_3 = h * f(y + k_2 / 2.0, t + h / 2.0).cast<double>();
    double k_4 = h * f(y + k_3, t + h).cast<double>();

    return y + ((k_1 + 2.0 * k_2 + 2.0 * k_3 + k_4) / 6.0);
    }
};

class RKF {
    public:
    class Values {
        public:
        double y;
        double h;
        bool accepted;
    };
    const double a11 = 0;
    const double a21 = 1.0 / 4;
    const double a31 = 3.0 / 32;
    const double a32 = 9.0 / 32;
    const double a41 = 1932.0 / 2197;
    const double a42 = -7200.0 / 2197;
    const double a43 = 7296.0 / 2197;
    const double a51 = 439.0 / 216;
    const double a52 = -8.0;
    const double a53 = 3680.0 / 513;
    const double a54 = -845.0 / 4104;
    const double a61 = -8.0 / 27;
    const double a62 = 2.0;
    const double a63 = -3544.0 / 2565;
    const double a64 = 1859.0 / 4104;
    const double a65 = -11.0 / 40;

    const double c1 = 0;
    const double c2 = 1.0 / 4;
    const double c3 = 3.0 / 8;
    const double c4 = 12.0 / 13;
    const double c5 = 5.0;
    const double c6 = 1.0 / 2;

    const double b11 = 16.0 / 135;
    const double b12 = 0;
    const double b13 = 6656.0 / 12825;
    const double b14 = 28561.0 / 56430;
    const double b15 = -9.0 / 50;
    const double b16 = 2.0 / 55;
    const double b21 = 25.0 / 216;
    const double b22 = 0.0;
    const double b23 = 1408.0 / 2565;
    const double b24 = 2197.0 / 4104;
    const double b25 = -1.0 / 5;
    const double b26 = 0.0;
    double epsilon;
    RKF(double e=0.005) {
        epsilon = e;
    }
    Values next_step(double y, double h, py::function f, double t) {
          double k1 = h * f(y, t + h * c1).cast<double>();
  double k2 = h * f(y + a21 * k1, t + h * c2).cast<double>();
  double k3 = h * f(y + a31 * k1 + a32 * k2, t + h * c3).cast<double>();
  double k4 = h * f(y + a41 * k1 + a42 * k2 + a43 * k3, t + h * c4).cast<double>();
  double k5 = h * f(y + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4, t + h * c5).cast<double>();
  double k6 = h * f(y + a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5,
                    t + h * c6).cast<double>();

  double rk4_y =
      y + b21 * k1 + b22 * k2 + b23 * k3 + b24 * k4 + b25 * k5 + b26 * k6;
  double rk5_y =
      y + b11 * k1 + b12 * k2 + b13 * k3 + b14 * k4 + b15 * k5 + b16 * k6;
  double e = std::abs(rk4_y - rk5_y);
  double factor = std::pow(epsilon / (e + 1e-10), 0.25);
  double h_new = h * std::clamp(0.84 * factor, 0.1, 5.0);
  Values result;
  result.h = h_new;
  result.accepted = (e <= epsilon);
  if (result.accepted) {
    result.y = rk5_y;
  } else {
    result.y = y;
  }
  return result;
    }
};

class DormandPrince {
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

const double b11 = 35.0/184;
const double b12 = 0;
const double b13 = 500.0/1113;
const double b14 = 125.0/192;
const double b15 = - 2187.0/6784;
const double b16 = 11.0/84;
const double b17 = 0;
const double b21 = 5179.0/57600;
const double b22 = 0.0;
const double b23 = 7571.0 / 16695;
const double b24 = 393.0 / 640;
const double b25 = -92097.0/339200;
const double b26 = 187.0 /2100;
const double b27 = 1.0 / 40;

double eps;
DormandPrince(double e = 0.005) {
    eps = e;
}
Values next_step(double y, double h,  double (*f)(double y, double t), double t) {
  double k1 = h * f(y, t + h * c1);
  double k2 = h * f(y + a21 * k1, t + h * c2);
  double k3 = h * f(y + a31 * k1 + a32 * k2, t + h * c3);
  double k4 = h * f(y + a41 * k1 + a42 * k2 + a43 * k3, t + h * c4);
  double k5 = h * f(y + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4, t + h * c5);
  double k6 = h * f(y + a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5,
                    t + h * c6);
  double k7 = h * f(y + a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6, t + h * c7);

  double rk4_y =
      y + b21 * k1 + b22 * k2 + b23 * k3 + b24 * k4 + b25 * k5 + b26 * k6 + b27 * k7;
  double rk5_y =
      y + b11 * k1 + b12 * k2 + b13 * k3 + b14 * k4 + b15 * k5 + b16 * k6 + b17 * k7;
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
};

double dydt(double y, double t) {
  return std::pow(y-1, 2) * std::pow(t-1, 2);
}

PYBIND11_MODULE(ODESolvers, handle) {
    handle.doc() = "My First Library!";
    py::class_<Euler>(handle, "Euler").def(py::init<>()).def("next_step", &Euler::next_step);
    py::class_<RK4>(handle, "RK4").def(py::init<>()).def("next_step", &RK4::next_step);
    py::class_<RKF::Values>(handle, "RKFValues")
        .def_readwrite("y", &RKF::Values::y)
        .def_readwrite("h", &RKF::Values::h)
        .def_readwrite("accepted", &RKF::Values::accepted);

    py::class_<RKF>(handle, "RKF")
        .def(py::init<double>(), py::arg("epsilon") = 0.005)
        .def("next_step", &RKF::next_step);
}