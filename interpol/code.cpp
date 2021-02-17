#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <cmath>
#include "poly.h"
#include "spline.h"

Poly newton_interpol(std::vector<double> x, std::vector<double> y) {
  int n = x.size();
  double ds[n][n];
  for (int k = 0; k < n; k++) {
    ds[0][k] = y[k];
  }

  for (int i = 1; i < n; i++) {
    for(int k = 0; k < n - i; k++) {
      ds[i][k] = (ds[i-1][k+1] - ds[i-1][k]) / (x[k+i] - x[k]);
    }
  }
  
  Poly inter, temp{1};
  for(int i = 0; i < n; i++) {
    inter = inter + ds[i][0] * temp;
    temp = temp * Poly{{-1 * x[i], 1}};
  }
  
  return inter;
}

template <class T>
void draw(std::vector<double> x, std::vector<double> y, T newt) {
  std::ofstream ofi("data.dat");
  std::ofstream giv("poly.dat");

  for(int i = 0; i < x.size(); i++) {
    giv << x[i] << " " << y[i] << std::endl;
  }
  for(double i = x[0]; i <= x[x.size() - 1]; i += 0.00001) {
    ofi << i << " " << newt.value(i) << std::endl;
  }

  FILE *gp = popen("gnuplot -persist", "w");
  fprintf(gp, "set grid x y\n");
  fprintf(gp, "show grid\n");
  fprintf(gp, "plot 'data.dat' with lines, 'poly.dat'\n");
  fprintf(gp, "pause mouse close\n");
  pclose(gp);
}

int main() 
{
  std::vector<double> x = {0.17453, 0.5236, 0.87267, 1.22173, 1.5708, 1.91986};
  std::vector<double> y = {38e-6, 0.00052, 0.00254, 0.01013, 0.03636, 0.11699};

  Poly newt = newton_interpol(x, y);
  std::cout << newt;

  Spline s(newt, x);
  std::cout << s;

  std::cout << newt(1.55) << std::endl;
  std::cout << s.value(1.55) << std::endl;

  std::vector<double> newx = {-4, -2, -1, 1, 2, 4};
  std::vector<double> newy = {-256, -8, 2, 4, 8, -128};

  draw(x, y, newt);

  return 0;
}
