#include <eigen3/Eigen/Dense>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <climits>
#include <iomanip>

Eigen::ArrayXd tridiag_solve(Eigen::ArrayXd &a, Eigen::ArrayXd &b, Eigen::ArrayXd &c, Eigen::ArrayXd &d) {
  int n = d.size();
  Eigen::ArrayXd x(n);
  for(int i = 1; i < n; i++) {
    double w = a[i]/b[i - 1];
    b[i] = b[i] - w*c[i - 1];
    d[i] = d[i] - w*d[i - 1];
  }
  x[n - 1] = d[n - 1]/b[n - 1];
  for(int i = n - 2; i >= 0; i--) {
    x[i] = (d[i] - c[i]*x[i + 1])/b[i];
  } 
  return x;
}

int main(int argc, char *argv[]) {
  int N = atoi(argv[1]);
  int L = 150, M = 150;
  double hr = 1./L, hf = M_PI/2/M, tau = 1./N;
  double eps = 1e-4;
  std::cout << "L = M = " << L << std::endl;
  std::cout << "N = " << N << std::endl;

  std::vector<Eigen::ArrayXXd> u(N + 1, Eigen::ArrayXXd::Zero(L + 1, M + 1));
  Eigen::ArrayXXd uan(L + 1, M + 1);

  int it = N;
  //get analytic solution
  for(int l = 0; l <= L; l++) {
    for(int m = 0; m <= M; m++) {
      uan(l, m) = (l*hr)*(l*hr)*sin(m*hf)*sin(m*hf)/(7 - 6*(tau*(it)));
    }
  }

  //boundary conditions
  //first layer
  for(int l = 0; l <= L; l++) {
    for(int m = 0; m <= M; m++) {
      u[0](l, m) = (l*hr)*(l*hr)*sin(m*hf)*sin(m*hf)/7;
    }
  }
  //edges
  for(int n = 0; n <= N; n++) {
    for(int k = 0; k <= M; k++) {
      u[n](0, k) = 0;
      u[n](L, k) = sin(k*hf)*sin(k*hf)/(7 - 6*(tau*n));
      u[n](k, M) = (k*hr)*(k*hr)/(7 - 6*(tau*n));
      u[n](k, 0) = 0;
    }
  }

  //cycle over time
  
  for(int n = 0; n < N; n++) {
  /*
    if(n % 10 == 0) {
      std::cout << n << " iteration over n" << std::endl;
    }
  */
    Eigen::ArrayXXd uk = u[n];
    Eigen::ArrayXXd ukk = uk;
    do {
      uk = ukk;
      //first subsystem
      Eigen::ArrayXXd ut = Eigen::ArrayXXd::Zero(L + 1, M + 1);
      for(int m = 1; m < M; m++) {
        Eigen::ArrayXd a(L + 1), b(L + 1), c(L + 1), d(L + 1);
        for(int l = 1; l < L; l++) {
          a[l] = -(hr*(l + 1./2))*(uk(l + 1, m) + uk(l, m))*tau/2/(l*hr)/hr/hr;
          c[l] = -(hr*(l - 1./2))*(uk(l, m) + uk(l - 1, m))*tau/2/(l*hr)/hr/hr;
          b[l] = 1 - a[l] - c[l];
          d[l] = u[n](l, m);
        }
        a[0] = a[L] = 0; c[0] = c[L] = 0;
        b[0] = 1; b[L] = 1; d[0] = u[n](0, m); d[L] = u[n](L, m);
        ut.col(m) = tridiag_solve(c, b, a, d);
      }
      ut.row(0) = u[n].row(0);
      ut.row(L) = u[n].row(L);
      ut.col(0) = u[n].col(0);
      ut.col(M) = u[n].col(M);
      //std::cout << ut << std::endl;
      //std::cout << "first solved" << std::endl;

      //second subsystem
      for(int l = 1; l < L; l++) {
        Eigen::ArrayXd a(L + 1), b(L + 1), c(L + 1), d(L + 1);
        for(int m = 1; m < M; m++) {
          a[m] = -(uk(l, m + 1) + uk(l, m))*tau/2/(hr*l)/(hr*l)/hf/hf;
          c[m] = -(uk(l, m) + uk(l, m - 1))*tau/2/(hr*l)/(hr*l)/hf/hf;
          b[m] = 1 - a[m] - c[m];
          d[m] = ut(l, m);
        }
        a[0] = a[M] = 0; c[0] = c[M] = 0;
        b[0] = 1; b[M] = 1; d[0] = ut(l, 0); d[M] = ut(l, M);
        ukk.row(l) = tridiag_solve(c, b, a, d);
      }
      ukk.row(0) = u[n + 1].row(0);
      ukk.row(L) = u[n + 1].row(L);
      ukk.col(0) = u[n + 1].col(0);
      ukk.col(M) = u[n + 1].col(M);
      //std::cout << "second solved" << std::endl;
      //std::cout << ukk << std::endl;
      //std::cout << ((ukk - uk)/ukk).abs().maxCoeff() << std::endl;

    } while((((ukk - uk)/ukk).abs() >= eps).any());
    u[n + 1] = ukk;
  }


  int step = L/5;
  std::cout << std::setw(9) << "(r,phi)" << " : " \
            << std::setw(9) << "u" << " | " << std::setw(9) << "u_an" \
            << " | " << std::setw(9) << "error" << std::endl;
  double md = __DBL_MIN__;
  for(int l = 0; l <= L; l += step) {
    for(int m = 0; m <= M; m += step) {
      double d = fabs(u[it](l, m) - uan(l, m));
      md = (d > md) ? d : md;
      std::cout << std::fixed << std::setprecision(1) << "(" << l*hr << "," << m*hf << ") : " << \
                   std::setprecision(7) << std::setw(9) << u[it](l, m) << " | " << std::setw(9) << uan(l, m) << " | " << \
                   std::setw(9) << d << std::endl;

    }
  }
  std::cout << "Maximum error value: " << md << std::endl << std::endl << std::endl;

  return 0;
}