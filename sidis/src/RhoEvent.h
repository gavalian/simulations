// RhoEvent.h
#ifndef RHO_EVENT_H
#define RHO_EVENT_H

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TwoBodyDecay.h"
#include <cmath>
#include <iostream>

class RhoEvent {

public:
    //--- inputs
  double Q2, xB, beamE, pol;
  static constexpr double Mp = 0.93827;
  //--- stored four‑vectors
    TLorentzVector e_in, e_out;
    TLorentzVector p_in, p_out;
    TLorentzVector q,   v;        // virtual γ and ρ
    TLorentzVector piPlus, piMinus;

    double eps , y, enu;
    //--- decay angles in ρ rest frame
    double theta;   // polar
    double phi;     // azimuthal

    //--- decay angles in ρ rest frame
    double thetaPi;   // polar
    double phiPi;     // azimuthal

    // ctor: set kinematics, beamE defaults to 10.6 GeV
    RhoEvent(double beamE_ = 10.6):Q2(2.0),xB(0.3),beamE(beamE_){
        rng.SetSeed(0);
        // target proton at rest
        p_in.SetPxPyPzE(0,0,0, Mp);
        // incoming e⁻ along z
        e_in.SetPxPyPzE(0,0, beamE, beamE);
    }
    RhoEvent(double Q2_, double xB_, double beamE_ = 10.6)
      : Q2(Q2_), xB(xB_), beamE(beamE_), theta(0), phi(0)
    {
        rng.SetSeed(0);
        // target proton at rest
        p_in.SetPxPyPzE(0,0,0, 0.938);
        // incoming e⁻ along z
        e_in.SetPxPyPzE(0,0, beamE, beamE);
    }

    // build full event + decay + angles
    void Generate() {
        BuildVirtualPhoton();
        BuildScatteredElectron();
        BuildRecoilProton();
        BuildRhoMeson();
        DecayRho();
        //ComputeDecayAngles();
    }

    // print everything
    void Print() const {
        TLorentzVector t = p_out-p_in;
        std::cout << "--- RhoEvent ---\n";
        std::cout << "e_in:    "; e_in.Print();
        std::cout << "e_out:   "; e_out.Print();
        std::cout << "p_in:    "; p_in.Print();
        std::cout << "p_out:   "; p_out.Print();
        std::cout << "q:       "; q.Print();
        std::cout << "rho(v):  "; v.Print();
        std::cout << "pi+:     "; piPlus.Print();
        std::cout << "pi-:     "; piMinus.Print();
        std::cout << " Q2  = " << Q2 << std::endl;
        std::cout << " xB  = " << xB << std::endl;
        std::cout << " y   = " << y << std::endl;
        std::cout << " eps = " << eps << std::endl;
        std::cout << " -t  = " << t.M2() << std::endl;
        std::cout << "theta:   " << theta*180/M_PI << "°\n";
        std::cout << "phi:     " << phi*180/M_PI   << "°\n";
        std::cout << "theta: pi  " << thetaPi*180/M_PI << "°\n";
        std::cout << "phi: pi    " << phiPi*180/M_PI   << "°\n";
    }

public:
    TRandom3 rng;

    // 1) virtual photon: q = e_in - e_out
    void BuildVirtualPhoton() {
        q = e_in-e_out;
    }

    TLorentzVector ScatteredElectron(double Q2, double xB, double beamE) {
    // 1) energy transfer ν = Q2 / (2 Mp xB)
        double nu = Q2 / (2.0 * Mp * xB);
        enu = nu;
    // 2) scattered‐electron energy E' = E_beam − ν
     double Eprime = beamE - nu;
       y     = nu/beamE;
     double γ     = gamma_bj(xB,Q2);
      eps   = epsilon(y,γ);
    // 3) scattering angle θ_e from Q2 = 2 E E' (1 - cos θ)
        double cosTh = 1.0 - Q2 / (2.0 * beamE * Eprime);
    // guard against small numerical slip:
        if (cosTh > +1) cosTh = +1;
        if (cosTh < -1) cosTh = -1;
        double sinTh = std::sqrt(1.0 - cosTh*cosTh);

        // 4) build four‐vector (px, py, pz, E)
        //    assume scattering in the x–z plane: py = 0
        double px = Eprime * sinTh;
        double pz = Eprime * cosTh;

        return TLorentzVector(px, 0.0, pz, Eprime);
    }

    // 2) scattered electron (simple elastic kinematics)
    void BuildScatteredElectron() {
        e_out = ScatteredElectron(Q2,xB,beamE);
    }

    // 3) recoil proton from t
    void BuildRecoilProton() {
        TLorentzVector  cm = q + p_in;
        //std::cout << " mass = " << cm.M2() << "  q2 " << q.M2() << std::endl;
        TVector3       ref = e_out.Vect();
        TwoBodyDecay decay (cm , Mp, 0.77526);
        decay.Generate(ref);
        p_out = decay.daughter1;
        //std::cout << " d1 = " << decay.daughter1.M() << " d2 " << decay.daughter2.M() << std::endl;
        v = decay.daughter2;
        theta = decay.theta;
        phi = decay.phi;
    }
    void BuildRecoilPions() {

        TLorentzVector  vrho = v;
        TVector3         ref = q.Vect();
        TwoBodyDecay decay (vrho , 0.139570, 0.139570);
        decay.Generate(ref);
        piPlus = decay.daughter1;
        //std::cout << " d1 = " << decay.daughter1.M() << " d2 " << decay.daughter2.M() << std::endl;
        piMinus = decay.daughter2;
        thetaPi = decay.theta;
        phiPi = decay.phi;
    }
    // 4) ρ meson from 4‑momentum conservation
    void BuildRhoMeson() {
        //v = q + p_in - e_out - p_out;
        v = q + p_in - p_out;
    }

    // 5) isotropic 2‑body decay in ρ rest frame
    void DecayRho() {
        const double Mpi = 0.13957;
        double Mrho = v.M();
        double pMag = sqrt((Mrho*Mrho - 4*Mpi*Mpi))/2.;
        double cos_t = rng.Uniform(-1,1);
        double sin_t = sqrt(1 - cos_t*cos_t);
        double ph    = rng.Uniform(0, 2*M_PI);

        // momenta in ρ‑rest
        TLorentzVector p1(pMag*sin_t*cos(ph),
                          pMag*sin_t*sin(ph),
                          pMag*cos_t,
                          sqrt(pMag*pMag + Mpi*Mpi));
        TLorentzVector p2 = p1;
        p2.SetVect(-p1.Vect());  // back‑to‑back

        // boost back to lab
        TVector3 b = v.BoostVector();
        piPlus  = p1;  piPlus.Boost(b);
        piMinus = p2;  piMinus.Boost(b);
    }

    // 6) compute θ,φ in ρ rest frame with helicity axes

     double gamma_bj(double xB, double Q2)
    { return 2.0 * xB * Mp / std::sqrt(Q2); }

  double epsilon(double y, double γ)
 {
     double y2γ2 = y*y*γ*γ;
     return (1.0 - y - 0.25*y2γ2) /
            (1.0 - y + 0.5*y*y + 0.25*y2γ2);
 }
};

#endif // RHO_EVENT_H
