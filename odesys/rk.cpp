#include <iostream>
#include <vector>
#include <valarray>
#include <utility>
#include <functional>
#include <fstream>
#include <iomanip>
#include "rk.h"

std::valarray<double> rhs(double x, std::valarray<double> y) {
  return {-49*y[0] + 125*y[1], 20*y[0] - 49*y[1]};
}

int main() {
  double A = -5, B = -2;
  std::valarray<double> y0(2);
  y0[0] = A;
  y0[1] = B;
  const double x0 = 0, D = 20;
  rk4 Solver(x0, D, y0, rhs);
  Solver.solve_precision(10e-7);
  Solver.store_last();
  Solver.gen_table();
  Solver.store_error();
  return 0;
}
