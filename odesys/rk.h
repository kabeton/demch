#include <valarray>
#include <functional>
#include <fstream>
#include <iostream>
#include <vector>
#include <float.h>
#include <iomanip>

typedef std::function<std::valarray<double>(double, std::valarray<double>)> func;

double norm(std::valarray<double> el) {
  return abs(el).max();
}

double sol_dist(std::vector<std::valarray<double>> s1, std::vector<std::valarray<double>> s2) {
  double t = 0, tm = -DBL_MAX;
  for(int i = 0; i < s2.size(); i++) {
    t = norm(s1[i*2] - s2[i]);
    if(t > tm) tm = t;
  }
  return tm;
}

class solver
{
private:
  std::valarray<double> y0;
  std::vector<std::vector<std::valarray<double>>> sol_seq;
  std::vector<double> h_seq;
protected:
  int dims;
  double x0, xl;
  func rhs;
public:
  solver(double _x0, double _xl, std::valarray<double> _y0, func _rhs) {
    dims = _y0.size();
    x0 = _x0;
    xl = _xl;
    y0 = _y0;
    rhs = _rhs;
  };

  virtual std::valarray<double> step(int i, double h, std::valarray<double> pr) = 0;

  void solve_step(double h) {
    int iters = (xl - x0) / h;
    std::vector<std::valarray<double>> sol;
    sol.reserve((xl - x0) / h + 1); // с этими +-1 хз, на всякий пока сделаю так
    sol.push_back(y0);
    std::valarray<double> rp(dims);

    for(int i = 1; i < iters; i++) {
      rp = step(i, h, sol[i - 1]);
      sol.push_back(rp);
    }

    h_seq.push_back(h);
    sol_seq.push_back(sol);
  }

  void solve_precision(double eps) {
    double h0 = (xl - x0) / 10;
    solve_step(h0);
    double h = h0 / 2;
    solve_step(h);
    double local_err = eps + 1;
    while(local_err > eps) { 
      h /= 2;
      solve_step(h);
      local_err = sol_dist(sol_seq[sol_seq.size()-1], sol_seq[sol_seq.size() - 2]);
    }
  }

  void store_last() {
    std::ofstream out("trail.dat");
    for(int i = 0; i < (xl - x0) / h_seq.back(); i++) {
      out << x0 + i*h_seq.back() << " " << sol_seq.back()[i][0] << " " << sol_seq.back()[i][1] << std::endl;
    }
  }

  void gen_table() {
    std::ofstream table("table");
    table << std::setprecision(9);
    int sols = sol_seq.size();
    for(int i = 0; i < sols; i++) {
      table << "--------------------------------------------------------------------------" << std::endl;
      int k = 1;
      for(int j = 1*pow(2, i); j < sol_seq[i].size() - 1; j+= pow(2, i)) {
        double c0, c1, m0, m1;
        c0 = sol_seq[i][j][0];
        c1 = sol_seq[i][j][1];
        m0 = sol_seq.back()[k*pow(2,sols-1)][0];
        m1 = sol_seq.back()[k*pow(2, sols-1)][1];
        table << j << " | " << x0 + h_seq[i]*j << " | " << c0 << " | " << m0 << " | " << fabs(c0 - m0) << " | "
         << c1 << " | " << m1 << " | " << fabs(c1 - m1) << std::endl;
        k++;
      }
    }
  }

  ~solver(){};
};

class rk4: public solver
{
public:
  rk4(double _x0, double _xl, std::valarray<double> _y0, func _rhs): solver(_x0, _xl, _y0, _rhs){};

  std::valarray<double> step(int i, double h, std::valarray<double> pr) {
    std::valarray<double> f1(dims), f2(dims), f3(dims), f4(dims);
    f1 = rhs(x0 + i*h, pr);
    f2 = rhs(x0 + i*h + h/2, pr + h*f1/2);
    f3 = rhs(x0 + i*h + h/2, pr + h*f2/2);
    f4 = rhs(x0 + i*h + h, pr + h*f3);
    return pr + h/6*(f1 + 2*f2 + 2*f3 + f4);
  }
};

class imp_euler1: public solver
{
public:
  imp_euler1(double _x0, double _xl, std::valarray<double> _y0, func _rhs): solver(_x0, _xl, _y0, _rhs){};

  std::valarray<double> step(int i, double h, std::valarray<double> pr) {
    return pr; 
  }
};