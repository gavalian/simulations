#include "TwoBodyDecay.h"


#ifndef __RECATION__
#define __REACTION__

class reaction {

  static constexpr double Mp = 0.93827;

  double Mv = 1.02;
  
  double rq2, rxb, rbeam, rpol;
  
  TLorentzVector e_in,e_out;
  TLorentzVector p_in,p_out;
  TLorentzVector q,vec;
  TLorentzVector decayOne, decayTwo;

  double eps,  y,  enu;
  double decayTheta, decayPhi;
  double prodTheta, prodPhi;
  double decayDaughterMassOne,decayDaughterMassTwo;
    
  TRandom3 rng;
  
public:
  reaction(){}
  reaction(double beame_, double massv_){
    rbeam = beame_; Mv = massv_;
  }
  void set(double beame_, double massv_){
    rbeam = beame_; Mv = massv_;
  }

  void setDecay(double dm1, double dm2){
    decayDaughterMassOne=dm1; decayDaughterMassTwo = dm2;
  }
  
  void generate(double q2_, double xb_){
    rq2 = q2_; rxb = xb_;
    produce__(rbeam,rq2,rxb);
    decay();
  }
  
  void produce__(double beame_, double q2_, double xb_){
    rng.SetSeed(0);
    p_in.SetPxPyPzE(0,0,0, Mp);
    // incoming e⁻ along z
    e_in.SetPxPyPzE(0,0, beame_, sqrt(beame_*beame_+0.0005*0.0005));
    e_out = scatteredElectron(beame_,q2_,xb_);
    q = e_in-e_out;
    scatteredProton(Mv);
  }

  double Q2(){return rq2;}
  double xB(){return rxb;}
  double pol(){return rpol;}
  
  void printVector(const char *name,TLorentzVector& vec){
    printf("%s: %9.5f %9.5f %9.5f [%9.5f]\n",name,vec.Px(),vec.Py(),vec.Pz(),vec.M());
  }
  void show(){
    printf("header: production [%9.5f, %9.5f], decay [%9.5f %9.5f]\n",prodTheta,prodPhi,
	   decayTheta,decayPhi);
    
    printf("e  in: %9.5f %9.5f %9.5f [%9.5f]\n",e_in.Px(),e_in.Py(),e_in.Pz(),e_in.M());
    printf("e out: %9.5f %9.5f %9.5f [%9.5f]\n",e_out.Px(),e_out.Py(),e_out.Pz(),e_out.M());
    printf("p out: %9.5f %9.5f %9.5f [%9.5f]\n",p_out.Px(),p_out.Py(),p_out.Pz(),p_out.M());
    printf("decay: %9.5f %9.5f %9.5f [%9.5f]\n",vec.Px(),vec.Py(),vec.Pz(),vec.M());
    printVector("dgt 1", decayOne);
    printVector("dgt 2", decayTwo);
  }

  void scatteredProton(double mass){
    TLorentzVector  cm = q + p_in;
        //std::cout << " mass = " << cm.M2() << "  q2 " << q.M2() << std::endl;
    TVector3       ref = e_out.Vect();
    TwoBodyDecay decay (cm , Mp, mass);
    decay.Generate(ref);
    p_out = decay.daughter1;
    //std::cout << " d1 = " << decay.daughter1.M() << " d2 " << decay.daughter2.M() << std::endl;
    vec = decay.daughter2;
    prodTheta = decay.theta;
    prodPhi = decay.phi;
  }
  
  TLorentzVector scatteredElectron( double beamE,double Q2, double xB){
    double nu = Q2 / (2.0 * Mp * xB);
    enu = nu;
    // 2) scattered‐electron energy E' = E_beam − ν
    double Eprime = beamE - nu;
    y  = nu/beamE;
    double gam     = gamma_bj(xB,Q2);
    eps   = epsilon(y,gam);
    // 3) scattering angle θ_e from Q2 = 2 E E' (1 - cos θ)
    double cosTh = 1.0 - Q2 / (2.0 * beamE * Eprime);
    //printf("cosTh = %f, eprime = %f\n",cosTh, Eprime);
    // guard against small numerical slip:
    if (cosTh > +1) cosTh = +1;
    if (cosTh < -1) cosTh = -1;
    double sinTh = std::sqrt(1.0 - cosTh*cosTh);
    //printf("cosTh = %f\n",cosTh);
    // 4) build four‐vector (px, py, pz, E)
    //    assume scattering in the x–z plane: py = 0
    double px = Eprime * sinTh;
    double pz = Eprime * cosTh;
    
    return TLorentzVector(px, 0.0, pz, sqrt(Eprime*Eprime+0.0005*0.0005));
  }

  double gamma_bj(double xB, double Q2)
  { return 2.0 * xB * Mp / std::sqrt(Q2); }
  
  double epsilon(double y, double γ)
  {
    double y2γ2 = y*y*γ*γ;
    return (1.0 - y - 0.25*y2γ2) /
      (1.0 - y + 0.5*y*y + 0.25*y2γ2);
  }

  void decay(){
     double Mrho = vec.M();
     double pMag = sqrt((Mv*Mv - 4*decayDaughterMassOne*decayDaughterMassOne))/2.;
     double cos_t = rng.Uniform(-1,1);
     double sin_t = sqrt(1 - cos_t*cos_t);
     double ph    = rng.Uniform(0, 2*M_PI);
     TLorentzVector p1(pMag*sin_t*cos(ph),
		       pMag*sin_t*sin(ph),
		       pMag*cos_t,
		       sqrt(pMag*pMag + decayDaughterMassOne*decayDaughterMassOne));
     TLorentzVector p2 = p1;
     p2.SetVect(-p1.Vect());
     TVector3 b = vec.BoostVector();
     decayOne = p1; decayOne.Boost(b);
     decayTwo = p2; decayTwo.Boost(b);
     decayPhi = ph;
     decayTheta = acos(cos_t);
     
  }
};

class candidate {
public:
  reaction react;
  double  weight;
  candidate(double beam, double massv, double m1, double m2){
    react.set(beam,massv);weight=0.0;
    react.setDecay(m1,m2);
  }
  
};

#endif
