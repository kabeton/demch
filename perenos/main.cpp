#include <iostream>
#include <eigen3/Eigen/Dense>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cfloat>

int main(int argc, char *argv[]) {
  int NX = atoi(argv[1]);
  double h = 1. / NX, tau = h/6;
  int NT = (int)(1/tau);

  Eigen::ArrayXXd u(NX + 1, NT + 1);
  Eigen::ArrayXXd uan(NX + 1, NT + 1);

  //analitic solution
  for(int i = 0; i <= NT; i++) {
    for(int j = 0; j <= NX; j++) {
      uan(j, i) = j*h - i*tau + exp(i*tau);
    }
  }

  //initial conditions
  for(int i = 0; i <= NX; i++) {
    u(i, 0) = i*h + 1;
  }
  for(int i = 0; i <= NT; i++) {
    u(0, i) = exp(tau*i) - tau*i;
    u(1, i) = exp(tau*i) - tau*i + h;
    u(2, i) = exp(tau*i) - tau*i + 2*h;
  }

  //scheme realisation
  for(int n = 0; n <= NT - 1; n++) {
    for(int l = 3; l <= NX; l++) {
      u(l, n + 1) = u(l, n) + tau/6/h*(2*u(l - 3, n) - 9*u(l - 2, n) + 18*u(l - 1, n) - 11*u(l, n)) + \
                    tau*tau/2/h/h*(-u(l - 3, n) + 4*u(l - 2, n) - 5*u(l - 1, n) + 2*u(l, n)) - \
                    tau*tau*tau/6/h/h/h*(-u(l - 3, n) + 3*u(l - 2, n) - 3*u(l - 1, n) + u(l, n)) + \
                    tau*exp(tau*n)*(1 + tau/2 + tau*tau/6);
    }
  }

  //compare scheme and analitic solution
  int xst = NX / 10;
  double dmax = -DBL_MAX;
  std::cout << "t = " << (NT)*tau << std::endl;
  std::cout << std::setw(3) << " x " << std::setw(10) << " u computed " << std::setw(10) << " u analitic " \
            << std::setw(10) << " error" << std::endl;
  
  for(int i = 0; i <= 10; i++) {
    double diff = fabs(u(i*xst, NT) - uan(i*xst, NT));
    if(diff > dmax) dmax = diff;
  }

  int mord = ((-1*(int)log10(dmax) + 1) >= 0) ? (-1*(int)log10(dmax) + 1) : 8;
  std::cout << mord << std::endl;
  std::cout << std::setprecision(mord);
  for(int i = 0; i <= 10; i++) {
    double diff = fabs(u(i*xst, NT) - uan(i*xst, NT));
    std::cout << std::setw(3) << i*xst*h << " " << std::setw(mord + 2) << u(i*xst, NT) << " " << std::setw(mord + 2) \
              << uan(i*xst, NT) << " " << std::setw(mord + 2) << diff << std::endl;
  }

  std::cout << "max diff: " << dmax << std::endl;

  return 0;
}