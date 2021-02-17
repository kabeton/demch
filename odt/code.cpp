#include "solver.h"
#include <utility>
#include <vector>
#include <functional>
#include <cstdlib>
#include <fstream>

void draw(grid y) {
  std::ofstream ofi("data.dat");

  for(int i = 0; i < y.x.size(); i++) {
    ofi << y.x[i].first << " " << y.x[i].second << std::endl;
  }

  FILE *gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set grid x y\n");
  fprintf(gp, "show grid\n");
  fprintf(gp, "plot 'data.dat' with lines\n");
  fprintf(gp, "pause mouse close\n");
  pclose(gp);
}

double f(double x, double y) {
  return (y/x)*(y/x) + 1.5*y/x - 1/x;
}

int main() 
{
  double x0 = 1, xn = 11, y0 = 0.5, err = 0.0001;
  std::function<double(double, double)> eq;
  eq = f;
  rk4 Solver(x0, xn, y0, eq);
  grid y = Solver.eps_sol(err);
  //draw(y);
  std::vector<grid> res = Solver.eps_sol_saved(err);
  grid trace = res.back();
  std::vector<double> Dist;
  std::cout << "Results for n= " << trace.x.size() << std::endl;
  for(int i = 0; i < trace.x.size(); i+=32) {
    std::cout << "x = " << trace.x[i].first << "  " << "y = " << trace.x[i].second << std::endl;
  }
  for(int i = 0; i < res.size() - 1; i++) {
    Dist.push_back(dist(trace, res[i]));
    std::cout << "n = " << res[i].x.size() << " distance: " <<  Dist[i] << std::endl;
  }
  return 0;
}
