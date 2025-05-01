#include "w_kernels.hpp"
#include <iostream>
using namespace Wkernels;

int main()
{
    Mat4 u{}, l{}, s{};          // Fill these with real data!
    double eps = 0.7;
    double phi = 0.5;            //   ϕ  in radians
    double kap = 0.3;            //   φ  in radians

    u[h('0')][h('0')][h('+')][h('+')] = {0.0123, 0.0000};   // real example
    l[h('+')][h('-')][h('0')][h('+')] = {0.0,    -0.0071};  // imaginary
    
    auto wUU = UU(u, eps, phi, kap);
    std::cout << "W_UU^LL = " << wUU.LL << '\n';
}
