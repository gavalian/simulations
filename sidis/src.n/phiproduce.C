#include "generator.h"

void phiproduce(){
  // define a reaction to decay to proton and rho
  // where rho decays to two pions
  //sim::candidate cr(10.6,0.77,0.13957,0.13957);
  sim::candidate cr(10.6,1.02,0.49368,0.49368);
  sim::generator gen(&cr);
  sim::event     event;
  
  // set range for generation, Q2 min, Q2 max, Xb min, Xb max
  gen.setRange(1.5,2.5,0.05,1.0);

  for(int j = 0; j < 24; j++) {
    gen.generate();
    //cr.react.show();
    cr.react.getEvent(event);
    event.show();
  }
}
