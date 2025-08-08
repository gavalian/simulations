// RhoLundIO.h ─────────────────────────────────────────────────────────
// Helper for printing a RhoEvent in LUND format (CLAS convention)

#ifndef RHO_CROSS_H
#define RHO_CROSS_H

#include "RhoEvent.h"
#include <iomanip>
#include <ostream>
#include <vector>

#include <iomanip>
#include <vector>
#include "w_kernels.hpp"

struct Candidate {
     RhoEvent  ev;
     double weight;
 };
 
 struct PhysicsInput {
     Wkernels::Mat4 u, l, s;
     double dsigmaT_dt, dsigmaL_dt;   // [nb / GeV²]
     double Pl = 0.0, SL = 0.0, ST = 0.0;   // polarisations
 };

class Cross {

  TRandom3 rng;
  double beamEnergy = 10.6;
  static constexpr double alem  = 1.0 / 137.035999084;

public:

  Cross(double beam = 10.6){ beamEnergy = beam; rng.SetSeed(0); }

  void generate(Candidate &cand, double Q2_min = 2.0, double Q2_max=2.05, double xB_min = 0.3, double xB_max=0.31){
     double q2 = rng.Uniform(Q2_min,Q2_max);
     double xb = rng.Uniform(xB_min,xB_max);
     double pl = rng.Uniform(-1.0,1.0);
     
     cand.ev.Q2  = q2;
     cand.ev.xB  = xb;
     cand.ev.pol = pl;
     
     cand.ev.beamE = beamEnergy;
     cand.ev.BuildScatteredElectron();
     cand.ev.BuildVirtualPhoton();
     cand.ev.BuildRecoilProton();
     cand.ev.BuildRecoilPions();
  }

   double dsigma_3fold(double xB, double Q2, double y,
                            double eps, double dsigmaT_dt, double dsigmaL_dt)
 {
     double pref = alem/(2.0*M_PI);
     double kin  = (y*y)/(1.0-eps)*(1.0-xB)/xB/Q2;
     return pref * kin * (dsigmaT_dt + eps*dsigmaL_dt);
 }

  double dsigma_7fold(double WUU,double WLU,double WUL,
                            double WLL,double WUT,double WLT,
                            double Pl,double SL,double ST,
                            double sigma3)
 {
     double S = WUU + Pl*WLU + SL*WUL + Pl*SL*WLL
                      + ST*WUT + Pl*ST*WLT;
     return sigma3 * S / (4.0*M_PI*M_PI);
 }


  void updateWeight(Candidate &cand){
    
    double t = (cand.ev.p_out-cand.ev.p_in).M2();
    PhysicsInput ph = loadPhysicsInputs(cand.ev.xB,cand.ev.Q2,t);
    auto WUU = Wkernels::UU(ph.u,cand.ev.eps, cand.ev.phi,
			    cand.ev.phiPi);
    
    auto WLU = Wkernels::LU(ph.u,cand.ev.eps, cand.ev.phi,
			    cand.ev.phiPi);
    
    double cosTheta = std::cos(cand.ev.thetaPi);
    double sinTheta = std::sin(cand.ev.thetaPi);

    double σ3 = dsigma_3fold(cand.ev.xB,cand.ev.Q2,cand.ev.y,cand.ev.eps,ph.dsigmaT_dt,ph.dsigmaL_dt);
    double W_LU = cosTheta*cosTheta*WLU.LL+std::sqrt(2)*cosTheta*sinTheta*WLU.LT+sinTheta*sinTheta*WLU.TT;
    double W_UU = cosTheta*cosTheta*WUU.LL+std::sqrt(2)*cosTheta*sinTheta*WUU.LT+sinTheta*sinTheta*WUU.TT;


    ph.Pl = cand.ev.pol<0?-1:1;
    
    double w  = dsigma_7fold(W_UU,W_LU,0,0,0,0,
			     ph.Pl,ph.SL,ph.ST, σ3);

    //double wcheck = ph.Pl*W_LU;
    /*std::cerr << WLU.LL << " " << WLU.LT << "  " << WLU.TT << " W LU = "
      << W_LU << " sigma3 = " << σ3 << " w =  " << w << "  pol = " << ph.Pl << std::endl;     */
    //std::cout << "  w uu "  << W_UU <<  " weight = " << w << "  sigma 3 " <<  σ3 <<  std::endl;
    cand.weight = w;
    //cand.weight = σ3;
    //cand.weight = 1.0;
    if(w<0) {
      std::cerr << " Error: the weight of the event for q2 = " << cand.ev.Q2 << " / xb = "
		<< cand.ev.xB << " POL = " << ph.Pl << "  σ = " << w
		<<  " , σ3 = " << σ3 << std::endl;
      cand.weight = 0.0;
    }
  }
  
  PhysicsInput loadPhysicsInputs(double xB, double Q2,
                                double t /*GeV²*/)
  {
    PhysicsInput ph{};
    
    /* ---------------------------------------------------------------
     *  TODO:
     *  1.  Set the u,l,s helicity matrices for the chosen (xB,Q²,t)
     *      or read them from a file.
     *  2.  Provide dσT/dt and dσL/dt (nb/GeV²).
     * --------------------------------------------------------------*/
    ph.dsigmaT_dt = 1.0;   // placeholder
    ph.dsigmaL_dt = 1.0;   // placeholder
    
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('+')][Wkernels::h('+')] = {1.00, 0.0000};   // real example
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')] = {1.00, 0.0000};   // real example
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('+')] = {0.25, 0.2000};   // real example
    ph.u[Wkernels::h('0')][Wkernels::h('0')][Wkernels::h('-')][Wkernels::h('+')] = {0.00, 0.0000};
    
    return ph;
 }
};


#endif /* RHO_LUND_IO_H */
