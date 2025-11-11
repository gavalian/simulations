#include "generator.h"

void generate2(){
  // define a reaction to decay to proton and rho
  // where rho decays to two pions
  candidate cr(10.6,0.77,0.13957,0.13957);
  generator gen(&cr);
  // set range for generation, Q2 min, Q2 max, Xb min, Xb max
  gen.setRange(1.5,2.5,0.05,1.0);

  for(int j = 0; j < 12000; j++) {
    gen.generate();
    cr.react.show();
  }
}
