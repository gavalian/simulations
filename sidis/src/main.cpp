// ----
// compile g++ -O2 -std=c++17 main.cpp -I. `root-config --cflags --glibs` -o main
#include "RhoEvent.h"
#include "RhoLundIO.h"
#include "Cross.h"

int main(int argc,char* argv[]) {

     const long     Nrequest = (argc>1) ? std::stol(argv[1]) : 200;
     const double   Ebeam    = (argc>2) ? std::stod(argv[2]) : 10.6;
    // choose your kinematics:
    double Q2 = 2.0;      // GeV^2
    double xB = 0.3;      
    //double t  = -0.35;    // GeV^2
    Cross cross;
    Candidate cand;
TRandom3 R(0);
    cross.generate(cand);
    //cand.ev.Print();

    cross.updateWeight(cand);


    const int Nscan = 50000;
     double wMax = 0.0;
     for(int i=0;i<Nscan;++i) {
         //double w = generateCandidate(R,Ebeam).weight;
         cross.generate(cand);
         cross.updateWeight(cand);
         double w = cand.weight;
         if(w>wMax) wMax = w;
     }
     wMax *= 1.2;                                // safety margin

     std::cerr << "Cross section max value = " << wMax << std::endl;
     //std::cout << " Max Cross = " << wMax << std::endl;
     //return 1;
     std::cout << std::fixed << std::setprecision(6);
     long Nkept=0, Nattempt=0;
     while(Nkept < Nrequest) {
        //auto cand = generateCandidate(R,Ebeam,0.299,0.301,2.49,2.51,-0.51,-0.49);
         //auto cand = generateCandidate(R,Ebeam);
         cross.generate(cand);
         cross.updateWeight(cand);
         ++Nattempt;
         if( R.Uniform(0.,wMax) <= cand.weight ) {
             ++Nkept;
             printLundCLAS12(cand.ev,12);
             //const auto& e = cand.ev;
             //write_lund(e);
             /* -------  write LUND  --------------------------------- *
              * You can of course switch to ROOT/hipo here.           */
             /*std::cout << 6 << ' ' << 1 << '\n'   // #part, event# (dummy)
                       << " 1  -1  "  << e.l.Px()  << ' ' << e.l.Py()  << ' ' << e.l.Pz()  << ' ' << e.l.E()  << '\n'
                       << " 1  -1  "  << e.lp.Px() << ' ' << e.lp.Py() << ' ' << e.lp.Pz() << ' ' << e.lp.E() << '\n'
                       << " 1   1  "  << e.k1.Px() << ' ' << e.k1.Py() << ' ' << e.k1.Pz() << ' ' << e.k1.E() << '\n'
                       << " 1  -1  "  << e.k2.Px() << ' ' << e.k2.Py() << ' ' << e.k2.Pz() << ' ' << e.k2.E() << '\n'
                       << " 1 2212  "<< e.pp.Px() << ' ' << e.pp.Py() << ' ' << e.pp.Pz() << ' ' << e.pp.E() << '\n'
                       << " 1 2212  0 0 0 " << Mp << '\n';  // target
                       */
         }
     }
     std::cerr << "Generated " << Nkept << " events in "
               << Nattempt << " attempts  (ε = "
               << 100.*Nkept/Nattempt << " %)\n";

    /*RhoEvent ev(Q2, xB);
    for(int i = 0; i < 200000; i++){
    ev.BuildScatteredElectron();
    ev.BuildVirtualPhoton();
    ev.BuildRecoilProton();
    ev.BuildRecoilPions();
    //ev.BuildRhoMeson();
    TLorentzVector t = ev.p_in - ev.p_out;

    //ev.Generate();
    //ev.Print();

    //std::cout << " t = " << t.M2() << "  " << ev.v.M() << "  " << ev.p_out.M() << std::endl;

        printLundCLAS12(ev,12);
    }*/
    return 0;
}
