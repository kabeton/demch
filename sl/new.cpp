#include <cstdio>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <array>
#include <utility>
#include "atridiag.h"
#include "logic.h"

double p(double x) 
{
  return x*x/(1 + x);
}

double q(double x) 
{
  return 2 - x;
}

int main() 
{
  double eps = 1e-4;
  int L = 10;
  double h = 1/(double)L;
  double l1 = 0.1, l2 = 100.1;

  std::array<double, 4> lams;
  std::array<double, 4> lt = {4.9348022, 44.4132198, 123.370055, 241.805307};
  model_container sample{'b', h, l1, L + 1, q, p};
  lams = get_l(L, h, q, p, sample);
  printf("l1 = %f\nl2 = %f\nl3 = %f\nl4 = %f\n", lams[0], lams[1], lams[2], lams[3]);

  L *= 2;
  h = 1/(double)L;

  std::array<double, 4> lams2;
  lams2 = get_l(L, h, q, p, sample);
  printf("l1 = %f\nl2 = %f\nl3 = %f\nl4 = %f\n", lams2[0], lams2[1], lams2[2], lams2[3]);

  while (dist(lams, lams2) > 1e-4 && L < 2000000) {
    lams = lams2;
    L *= 2;
    h = 1/(double)L;
    lams2 = get_l(L, h, q, p, sample);
    printf("\nL = %d\n", L);
    printf("l1 = %f\nl2 = %f\nl3 = %f\nl4 = %f\n", lams2[0], lams2[1], lams2[2], lams2[3]);
    //printf("dist = %f\n", dist(lams, lams2));
  }
  return 0;
}