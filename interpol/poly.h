#ifndef POLY_H
#define POLY_H 1

class Poly {
  int deg;
  std::vector<double> a;
public:
  Poly() {
    this->deg = 0;
    this->a = {0};
  }

  Poly(double n) {
    this->deg = 0;
    this->a = std::vector<double>(1, n);
  }

  Poly(std::vector<double> const &vec) {
    this->deg = vec.size() - 1;
    this->a = vec;
  }

  int getdeg() const {
    return deg;
  }
  
  double getdegc(int n) const {
    if(n <= deg) {
      return a[n];
    } else {
      std::cout << "error: " << n << std::endl;
      return -1;
    }
   }
   
  double value(double x) const {
    double val = a[0];
    for(int i = 1; i <= deg; i++) {
      double pow = 1;
      for(int j = 0; j < i; j++) {
        pow *= x;
      }
      val += a[i] * pow;
    }
    return val;
  }

  Poly derivative() const {
    if(deg == 0) return Poly();
    std::vector<double> newv(deg);
    for(int i = 0; i < deg; i++) {
      newv[i] = a[i + 1] * (i + 1);
    }
    return Poly(newv);
  }

  Poly operator+(Poly const &b) const {
    if(this->deg == b.getdeg()) {
      std::vector<double> newv(this->deg+1, 0);
      for(int i = 0; i <= this->deg; i++) {
        newv[i] += b.getdegc(i) + this->a[i];
      }
      return Poly(newv);
    } else {
      int n = std::max(this->deg, b.getdeg());
      std::vector<double> newv(n+1, 0);
      for(int i = 0; i <= this->deg; i++) {
        newv[i] += this->a[i];
      }
      for(int i = 0; i <= b.getdeg(); i++) {
        newv[i] += b.getdegc(i);
      }
      return Poly(newv);
    }
  }

  Poly operator=(Poly const &b) {
    this->deg = b.getdeg();
    std::vector<double> newv(b.getdeg()+1);
    for(int i = 0; i < newv.size(); i++) {
      newv[i] = b.getdegc(i);
    }
    this->a = newv;
    return *this;
  }

  Poly operator*(double const &b) const {
    std::vector<double> newv(this->deg+1);
    for(int i = 0; i <= this->deg; i++) {
      newv[i] = this->a[i] * b;
    }
    return Poly{newv};
  }

  friend Poly operator*(const double &b, const Poly &p) {
    return p * b;
  }

  Poly operator*(Poly const &b) const {
    int s = this->deg + b.getdeg();
    std::vector<double> newv(s+1, 0);
    for(int i = 0; i != this->deg + 1; i++) {
      for(int j = 0; j != b.getdeg() + 1; j++) {
        newv[i + j] += this->a[i] * b.getdegc(j);
      }
    }
    return Poly{newv};
  }

  double operator()(double x) {
    return value(x);
  }

  ~Poly() {};
};


std::ostream& operator<<(std::ostream &os, const Poly x) {
  for(int i = 0; i <= x.getdeg(); i++) {
    os << x.getdegc(i) << " " << "x^" << i << " ";
  }
  os << std::endl;
  return os;
}

#endif
