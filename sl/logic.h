#include "atridiag.h"
#include <cstdio>
#include <cmath>
#include <utility>
#include <array>
#include <float.h>

template<typename container>
double pinpoint(double l1, double l2, double h, int L,
                atridiag<container> &matl1, 
                atridiag<container> &matl2) {
  double d1 = matl1.det(), d2 = matl2.det();
  double ln = (l1 + l2)/2;
  double lprev = 0;

  while(fabs(ln - lprev) > 1e-8) {
    //printf("l = %f\n", (l1 + l2)/2);
    //printf("d1 = %f d2 = %f d1-d2 = %f\n", d1, d2, fabs(d1 - d2));
    double l3 = (l1 + l2)/2;
    atridiag<container> matl3 = matl1;
    matl3.update(l3);
    double d3 = matl3.det(); 
    if(d3*d1 > 0) {
      l1 = l3;
      d1 = d3;
    } else {
      l2 = l3;
      d2 = d3;
    }
    lprev = l3;
    ln = (l1 + l2)/2;
  }
  return (l1 + l2)/2;
}

template<typename container>
void locate_next(double &l1, double &l2,
                atridiag<container> &m1, 
                atridiag<container> &m2) 
{
  while(m1.det()*m2.det() > 0) {
    //printf("l1 = %f l2 = %f d1d2 = %f\n", l1, l2, m1.det()*m2.det());
    l1 += 10; 
    l2 += 10;
    m1.update(l1);
    m2.update(l2);
  }
}

template<typename container>
std::array<double, 4> get_l(int L, double h, std::function<double(double)> q, std::function<double(double)> p, 
                                  container sample_container) 
{
  double l1 = 0.1, l2 = 10.1;
  std::array<double, 4> ans;

  atridiag<container> matl1(L, 
                           container{'a', h, l1, L - 1, q, p},
                           container{'b', h, l1, L + 1, q, p},
                           container{'c', h, l1, L - 1, q, p});

  atridiag<container> matl2(L, 
                           container{'a', h, l2, L - 1, q, p},
                           container{'b', h, l2, L + 1, q, p},
                           container{'c', h, l2, L - 1, q, p});


  int counter = 0;
  while(counter < 4) {
    locate_next(l1, l2, matl1, matl2);
    double l = pinpoint(l1, l2, h, L, matl1, matl2);
    ans[counter] = l;
    l1 += 10;
    l2 += 10;
    matl1.update(l1);
    matl2.update(l2);
    counter++;
  }
  return ans;
}

double dist(std::array<double, 4> a1, std::array<double, 4> a2) {
    double t = 0, tm = -DBL_MAX;
    tm = fabs(a1[0] - a2[0]);
  for(int i = 0; i < 4; i++) {
    t = fabs(a1[i] - a2[i]);
    printf("dist%d = %f\n", i, t);
    if(t > tm) tm = t;
  }
  printf("dist estimation = %f\n\n", tm);
  return tm;
}

