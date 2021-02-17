#include <iostream>
#include <math.h>
#include <vector>
#include <cstdio>

std::vector<double> linspace(double start, double end, int n) {
  std::vector<double> res;
  res.reserve(n);
  double d = end - start;
  for(int i = 0; i < n; i++) {
    res.push_back(start + i * d / n);
  }
  return res;
}

std::vector<double> count_roots(std::vector<double> grid, std::vector<double> poly, std::vector<double> k) {
  std::vector<double> res;
  res.push_back(0.0);
  for (int i = 1; i < grid.size(); i++) {
    if (poly[i] * poly[i-1] < 0) {
      res[0]++;
      res.push_back(grid[i-1]);
      res.push_back(grid[i]);
    }
  }
  return res;
}

double poly(double c, std::vector<double> &k) {
  return k[0]*pow(c, 14) + k[1]*pow(c, 9) + k[2]*pow(c, 8) + k[3]*pow(c, 7) + k[4]*c*c + k[5]*c + k[6];
}

int main() {
  double g0 = 5.0 / 3;
  double r0 = 0.00001;
  double p0 = 3848;
  double u0 = 0, u3 = 0;
  double g3 = 7.0 / 5;
  double c3 = 25324.8;
  double p3 = 3.04 * 1000000000;
  
  double a0 = (g0 + 1) / (g0 - 1);
  double n = 2 * g3 / (g3 - 1);
  double mu = 0;
  double r3 = g3 * p3 / c3 / c3;
  double nu = 2 / (g3 - 1) * sqrt(r0 / r3 * p3 / p0 * 
      g3 * (g0 - 1) / 2);
  double X = p3 / p0;
  
  std::vector<double> k;
  k.reserve(7);
  k[0] = X * X;
  k[1] = -a0 * nu * nu * X;
  k[2] = 2 * a0 * nu * (mu + nu) * X;
  k[3] = -(2 + (mu + nu) * (mu + nu) * a0) * X;
  k[4] = -nu * nu;
  k[5] = 2 * nu * (mu + nu);
  k[6] = 1 - (mu + nu) * (mu + nu); 

  //for (int i = 0; i < 7; i++) {
  //  std::cout << k[i] << std::endl;
  //}
  //std::cout << std::endl;

  double zl, zr;
  zl = abs(k[6]) / (abs(k[6]) + abs(k[0]));
  zr = 1 + abs(k[2]) / abs(k[0]);

  std::cout << zl << " <= |z| <= " << zr << std::endl;
  std::vector<double> grid = linspace(zl, zr, 100);
  std::vector<double> poly;
  poly.reserve(grid.size());
  for (int i = 0; i < grid.size(); i++) {
    poly.push_back(k[0]*pow(grid[i], 14) + k[1]*pow(grid[i], 9) + k[2]*pow(grid[i], 8) +
        k[3]*pow(grid[i], 7) + k[4]*grid[i]*grid[i] + k[5]*grid[i] + k[6]);                  
  }

  std::vector<double> roots = count_roots(grid, poly, k);
  
  std::cout << roots[0] << " possible roots: " << std::endl;
  for(int i = 0; i < roots[0]; i++) {
    std::cout << roots[i*2 + 1] << " " << roots[2 + i*2] << std::endl;
  }
  
  double eps = 0.00001;
  double left = roots[1], right = roots[2], len = right - left;
  
  while(len > eps) {
    double c = (left + right) / 2;
    if((k[0]*pow(left, 14) + k[1]*pow(left, 9) + k[2]*pow(left, 8) + k[3]*pow(left, 7) + k[4]*left*left + k[5]*left + k[6]) * 
        (k[0]*pow(c, 14) + k[1]*pow(c, 9) + k[2]*pow(c, 8) + k[3]*pow(c, 7) + k[4]*c*c + k[5]*c + k[6]) > 0) {
      left = c;
    } else {
      right = c;
    }
    len = right - left;
  }
  double c = (left + right) / 2;
  printf("%g\n", c);
  double u1 = u3 + 2*c3/(g3 - 1)*(1 - c);
  printf("%g\n", u1);
  double d0 = ((p3*pow(c, n) - p0)/r0 - u0*u0 + u0*u1)/(u1-u0);
  printf("%g\n", d0);


  return 0;
}

