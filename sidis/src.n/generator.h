#include "reaction.h"

#ifndef __GENERATOR__
#define __GENERATOR__

class generator {
 private:

  static constexpr double Mp = 0.93827;
  
  candidate *cand;
  TRandom3   rand;
  
  double q2_min, q2_max;
  double xb_min, xb_max;
  double E;
 public:
  generator(candidate *c){ cand = c; E = cand->react.E();}
  
  void setRange(double __q2min, double __q2max, double __xbmin, double __xbmax){
    q2_min = __q2min; q2_max = __q2max;
    xb_min = __xbmin; xb_max = __xbmax;
  }

  int  is_valid(double q2, double xb){
    double upper = 2.0*Mp*E*xb/(1+Mp*xb/E);
    double lower = 2.73*xb/(1.0-xb);
    if(q2>1.0&&q2>lower&&q2<upper) return 1;
      return 0;
  }
  
  void generate(){
    double q2 = rand.Uniform(q2_min,q2_max);
    double xb = rand.Uniform(xb_min,xb_max);
    int counter = 0;
    while(is_valid(q2,xb)==0){
      q2 = rand.Uniform(q2_min,q2_max);
      xb = rand.Uniform(xb_min,xb_max);
      counter++;
    }
    cand->react.generate(q2,xb);
    printf("random (Q2,xb) : %4d %8.5f %8.5f\n",counter,q2,xb);
  }
};

#endif
