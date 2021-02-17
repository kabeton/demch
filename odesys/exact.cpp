#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

int main() {
  std::ofstream exact("exact.dat");
  for(double x = 0.0; x <= 1.0; x += 0.001) {
    exact << x << " " << (1/10. - 1/4.)*(5)*exp(-99*x) + (1/10. + 1/4.)*5*exp(x) << " " << (1/10. - 1/4.)*(-2)*exp(-99*x) + (1/10. + 1/4.)*2*exp(x) << std::endl;
  }
  return 0;
}