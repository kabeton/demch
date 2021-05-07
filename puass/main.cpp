#include <eigen3/Eigen/Dense>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <climits>
#include <cfloat>

std::vector<int> cheb(int N) 
{
  std::vector<int> nseq;
  nseq.push_back(N);
  while(nseq.back() != 1) {
    if(nseq.back() % 2 == 0) {
      nseq.push_back(nseq.back() / 2);
    } else {
      nseq.push_back(nseq.back() - 1);
    }
  }

  std::vector<int> vp = {1};
  std::vector<int> vl;
  for(auto it = nseq.rbegin() + 1; it != nseq.rend(); it++) {
    auto n = it;
    if(it != nseq.rend() - 1) n++;
    if(*it % 2 == 0) {
      vl.resize(*it);
      if(*n == *it * 2) {
        for(int i = 0; i < vp.size(); i++) {
          vl[2*i] = vp[i];
          vl[2*i + 1] = 4 * *it/2 - vp[i];
        }
      } else {
        for(int i = 0; i < vp.size(); i++) {
          vl[2*i] = vp[i];
          vl[2*i + 1] = 4 * *it/2 + 2 - vp[i];
        }
      }
    } else {
      vl.resize(*it);
      for(int i = 0; i < vp.size(); i++) {
        vl[i] = vp[i];
      }
      vl.back() = 2 * (*it - 1)/2 + 1;
    }
    vp = vl;
  }
  return vl;
}

Eigen::ArrayXXd get_ansol(int L, int M) {
  Eigen::ArrayXXd ans(L + 1, M + 1);
  double hx = 1./L, hy = 1./M;
  for(int l = 0; l <= L; l++) {
    for(int m = 0; m <= M; m++) {
      ans(l, m) = sin(3*M_PI*l*hx)*sin(4*M_PI*m*hy);
    }
  }
  return ans;
}

Eigen::ArrayXXd get_sol(std::vector<int> chseq, int L, int M) {
  int N = chseq.size();
  double hx = 1./L, hy = 1./M;
  Eigen::ArrayXXd up = Eigen::ArrayXXd::Ones(L + 1, M + 1);
  for(int i = 0; i <= L; i++) {
    up(i, 0) = 0;
    up(i, M) = 0;
    up(0, i) = 0;
    up(L, i) = 0;
  }

  Eigen::ArrayXXd u = up;
  double maxm = 4*((cos(M_PI/2/L)*cos(M_PI/2/L))/hx/hx + (cos(M_PI/2/M)*cos(M_PI/2/M))/hy/hy);
  double minm = 4*((sin(M_PI/2/L)*sin(M_PI/2/L))/hx/hx + (sin(M_PI/2/M)*sin(M_PI/2/M))/hy/hy);

  for(int n = 0; n < N; n++) {
    double tn = 2/((maxm + minm) + (maxm - minm)*cos(M_PI*(2*chseq[n] - 1)/2/N));

    for(int l = 1; l < L; l++) {
      for(int m = 1; m < M; m++) {
        u(l, m) = up(l, m) + tn*((up(l + 1, m) - 2*up(l, m) + up(l - 1, m))/hx/hx + \
                                 (up(l, m + 1) - 2*up(l, m) + up(l, m - 1))/hy/hy + \
                                  25*M_PI*M_PI*sin(3*M_PI*l*hx)*sin(4*M_PI*m*hy));
      }
    }
    up = u;
  }

  return u;
}

double norm(Eigen::ArrayXXd f, Eigen::ArrayXXd s) {
  int step = (int)((f.rows() - 1)/5);
  double ma = -DBL_MAX;
  std::cout << std::setw(9) << "(x,y)" << " : " \
            << std::setw(9) << "u" << " | " << std::setw(9) << "u_an" \
            << " | " << std::setw(9) << "norm" << std::endl;
  std::cout << "---------------------------------------------" << std::endl;
  for(int l = 0; l < f.rows(); l += step) {
    for(int m = 0; m < f.cols(); m += step) {
      double d = fabs(f(l, m) - s(l, m));
      std::cout << "(" << std::fixed << std::setprecision(1) << l*(1./(f.rows() - 1)) << "," << m*(1./(s.rows() - 1)) << ") : " \
                << std::setprecision(5) << std::setw(9) << f(l, m) << " | " << std::setw(9) << s(l, m) \
                << " | " << std::setw(9) <<  d << std::endl;
      if(d > ma) ma = d;
    }
  }
  return ma;
}

int main(int argc, char *argv[]) {
  int L = 150, M = 150;
  double hx = 1./L, hy = 1./M;
  double eps=1e-4;
  double maxm = 4*((cos(M_PI/2/L)*cos(M_PI/2/L))/hx/hx + (cos(M_PI/2/M)*cos(M_PI/2/M))/hy/hy);
  double minm = 4*((sin(M_PI/2/L)*sin(M_PI/2/L))/hx/hx + (sin(M_PI/2/M)*sin(M_PI/2/M))/hy/hy);

  //волшебная формула с оценкой из учебника демча
  int N = (int)(log(2/eps)/log((sqrt(maxm) + sqrt(minm))/(sqrt(maxm) - sqrt(minm)))) + 1;
  if(argc == 2) N = atoi(argv[1]);
  std::cout << "N = " << N << std::endl;

  std::vector<int> iseq = cheb(N);

  Eigen::ArrayXXd u(L + 1, M + 1);
  Eigen::ArrayXXd uan(L + 1, M + 1);

  uan = get_ansol(L, M);
  u = get_sol(iseq, L, M);

  double dist = norm(u, uan);

  std::cout << "max distance: " <<  dist << std::endl;

  return 0;
}