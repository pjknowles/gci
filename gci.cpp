#include "gci.h"
#include "gciFile.h"
#include "gciHamiltonian.h"
#include "gciDeterminant.h"
#include "gciWavefunction.h"
#include "gciStringSet.h"
#include "gciExcitationSet.h"
#include "gciOrbitalSpace.h"
#include "FCIdump.h"
#include <iostream>
using namespace gci;

std::vector<double> gci::RSPT(std::vector<gci::Hamiltonian*>& hamiltonians , State &prototype, int maxOrder)
{
  std::vector<double> e(maxOrder+1,(double)0);
  if (hamiltonians.size() < 1) throw "not enough hamiltonians";
//  xout << "H0: " << *hamiltonians[0] << std::endl;
//  xout << "H1: " << *hamiltonians[1] << std::endl;
  Wavefunction w(hamiltonians[0],prototype.nelec,prototype.symmetry,prototype.ms2);
  Wavefunction g(w);
  g.diagonalHamiltonian(*hamiltonians[0]);
  size_t reference = g.minloc();
  e[0]=g.at(reference);
  g-=e[0];g.set(reference,(double)1);
//  xout << "Møller-Plesset denominators: " << g.str(2) << std::endl;
  gci::File h0file; g.put(h0file);
  w.set((double)0); w.set(reference, (double) 1);
  gci::File wfile; w.put(wfile,0);
  gci::File gfile;
  for (int k=0; k < (int) hamiltonians.size(); k++) {
    g.set((double)0);
    g.hamiltonianOnWavefunction(*hamiltonians[k],w);
    g.put(gfile,k);
  }
  for (int n=1; n < maxOrder; n++) {
    // construct  |n> = -(H0-E0)^{-1} ( -sum_k^{n-1} E_{n-k} |k> + sum_{k=n-h}^{n-1} H|k>) where h is the maximum order of hamiltonian
//    xout <<std::endl<<std::endl<<"MAIN ITERATION n="<<n<<std::endl;
    g.set((double)0);
    for (int k=0; k<n; k++) {
      w.get(wfile,k);
      if (n-k < (int) hamiltonians.size())
        g.hamiltonianOnWavefunction(*hamiltonians[n-k], w);
//      if (n-k < (int) hamiltonians.size())
//        xout << "g after H.w: " << g.str(2) <<std::endl;
      if (n == 1) e[1]=g.at(reference);
//        xout << "k, E:"<<k<<" "<<e[n-k]<<", g before -E.w: " << g.str(2) <<std::endl;
//        xout <<"w="<<w.str(2)<<std::endl;
      g -= e[n-k] * w;
//        xout << "k, E:"<<k<<" "<<e[n-k]<<", g after -E.w: " << g.str(2) <<std::endl;
    }
    w = -g;
    g.get(h0file);
//        xout <<std::endl<< "Perturbed wavefunction before precondition: " << w.str(2) <<std::endl;
    w.set(reference,(double)0);
    w /= g;
//    xout <<std::endl<< "Perturbed wavefunction, order="<<n<<": " << w.str(2) <<std::endl;
    w.put(wfile,n);
    for (int k=1; k < (int) hamiltonians.size(); k++) {
      if (n+k > maxOrder) break;
      g.get(gfile,k);
//      xout <<"gfile "<<g.str(2)<<std::endl;
//      xout <<"contribution from n="<<n<<", k="<<k<<" to E("<<n+k<<")="<<g*w<<std::endl;
      e[n+k]+=g*w;
    }
  }
  return e;
}

void gci::HamiltonianPrint(Hamiltonian &hamiltonian, State &prototype, int verbosity)
{
  Wavefunction w(&hamiltonian,prototype.nelec,prototype.symmetry,prototype.ms2);
  Wavefunction g(w);
  xout << std::endl << "Full Hamiltonian matrix"<<std::endl;
  if (verbosity >= 0) {
    for (size_t i=0; i < w.size(); i++)
    {
      w.set((double)0); w.set(i,(double)1);
      g.set((double)0);
      g.hamiltonianOnWavefunction(hamiltonian,w);
      for (size_t j=0; j < w.size(); j++)
        xout <<g.at(j)<< " ";
      xout <<std::endl;
    }
  }
}


#ifndef MOLPRO
//int main(int argc, char *argv[])
int main()
{
  try {
    FCIdump dump("FCIDUMP");
    //        OrbitalSpace os("FCIDUMP");
    //        xout <<"Orbital space:" << os << std::endl;
    Hamiltonian hh(&dump);
//    xout << "Hamiltonian: " <<hh.str()<<std::endl;
    Wavefunction w(&dump);
    //    OrbitalSpace os = hh;
    //    xout << "Orbital space: " << os.str(1) <<std::endl;
    //    exit(0);
    //    State ss(&dump);
    //    ss.orbitalSpace=&os;

    //    Determinant d1(&ss);

    //    d1.create(3);
    //    d1.create(-1);


    //    String s1(&ss);
    //    int phase;
    //    unsigned int orbital;
    //    orbital=1;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<std::endl;
    //    orbital=3;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<std::endl;
    //    orbital=2;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<std::endl;
    //    s1.first(3);
    ////    orbital=5;phase=s1.create(orbital); xout << "Add orbital " << orbital << "; phase=" <<phase <<"; string=" <<s1.printable() <<std::endl;
    //    int n=0;
    //    for (bool i=true; i && n<20; n++,i=s1.next()) {
    //        xout << "Advance string=" <<s1.printable() <<std::endl;
    //    }
    //    xout <<"Total number of string="<<n<<std::endl;

    //    xout << "Scan through constructing determinants:" << std::endl;
    //    d1.first(); while(d1.next()) {
    //        xout << " Determinant " <<d1.printable() <<std::endl;
    //    }
    //    xout <<"done scanning through determinants"<<std::endl;

    //    xout << "Wavefunction after constructor:"<<w.str(2)<<std::endl
    //         <<"...end of Wavefunction after constructor."<<std::endl<<std::endl;
    //    w.set((double)0.12345);
    //    xout << "Wavefunction after assign:"<<w.str(2)<<std::endl
    //         <<"...end of Wavefunction after assign."<<std::endl<<std::endl;
    //    //    w.buildStrings();
    //    //    xout << "Wavefunction after buildStrings:"<<w.str(1)<<std::endl;
    //    Wavefunction w2=w;
    //    xout << "Copied wavefunction:"<<w2.str(2)<<std::endl
    //         <<"...end of copied wavefunction."<<std::endl<<std::endl;
    //    xout << "Original wavefunction after copy:"<<w.str(2)<<std::endl
    //         <<"...end of original wavefunction."<<std::endl<<std::endl;
    //    w.set((double)1);
    //    xout << "Original wavefunction after original changed:"<<w.str(2)<<std::endl
    //         <<"...end of original wavefunction."<<std::endl<<std::endl;
    //    xout << "Copied wavefunction after original changed:"<<w2.str(2)<<std::endl
    //         <<"...end of copied wavefunction."<<std::endl<<std::endl;

    //    xout << "w.w=" << w*w << std::endl;
    //    xout << "w2.w2=" << w2*w2 << std::endl;
    //    xout << "Original wavefunction after w2.w2:"<<w.str(2)<<std::endl
    //         <<"...end of original wavefunction."<<std::endl<<std::endl;

    //    //Wavefunction w3;
    //    xout << "w:"<<w.str(2)<<std::endl <<"...end w."<<std::endl<<std::endl;
    //    w2 = w;
    //    xout << "after w2=w, w:"<<w.str(2)<<std::endl <<"...end w."<<std::endl<<std::endl;
    //    xout << "after w2=w, w2:"<<w2.str(2)<<std::endl <<"...end w2."<<std::endl<<std::endl;
    //    Wavefunction w3 = w2+(double)98*w-(w*((double)99));
    //    xout << "back from w3=..." <<std::endl;
    //    xout << "w:"<<w.str(2)<<std::endl <<"...end w."<<std::endl<<std::endl;
    //    xout << "w2:"<<w2.str(2)<<std::endl <<"...end w2."<<std::endl<<std::endl;
    //    xout << "w3:"<<w3.str(2)<<std::endl <<"...end w3."<<std::endl<<std::endl;

    //    for (unsigned int syma=0; syma<8; syma++) {
    //      //        unsigned int symb = syma ^ w.symmetry;
    //      std::vector<ExcitationSet> seta;
    //      seta = w.alphaStrings[syma].allExcitations(w.alphaStrings[syma],1,1);
    //      xout << "Excitations from alpha strings of symmetry " << syma+1 <<std::endl;
    //      for (std::vector<ExcitationSet>::iterator a=seta.begin(); a!=seta.end(); a++)
    //        xout <<"ExcitationSet: " <<a->str()<<std::endl;
    //      xout << "Alpha occupation numbers"<<std::endl;
    //      std::vector<double> on = w.alphaStrings[syma].occupationNumbers();
    //      for (size_t i=0; i < w.alphaStrings[syma].size(); i++) {
    //        for (size_t j=0; j < w.orbitalSpace->total(); j++)
    //          xout << " " << on[i+j*w.alphaStrings[syma].size()];
    //        xout << std::endl;
    //      }
    //    }


    //    w.diagonalHamiltonian(hh);
    //    xout << "Diagonal elements: " << w.str(2) << std::endl;
    //    size_t i = w.minloc();
    //    Determinant d = w.determinantAt(i);
    //    xout << "Lowest determinant " << d <<" with energy "<<w.at(i)<<std::endl;

    //    Hamiltonian fh = hh.FockHamiltonian(d);
    //    xout << "Fock hamiltonian: " << fh.str(3) << std::endl;

    //    w2.set((double)0); w2.set(w.minloc(), (double) 1);

    //    xout << "trial wavefunction: " << w2.str(2) <<std::endl;

    //    w3.hamiltonianOnWavefunction(fh, w2);
    //    xout << "action of Fock hamiltonian on trial wavefunction: " << w3.str(2) <<std::endl;

    //    w3.hamiltonianOnWavefunction(hh, w2);
    //    xout << "action of hamiltonian on trial wavefunction: " << w3.str(2) <<std::endl;

    //    xout << "Start looking for annihilation spaces using w=" << w.str(5) << std::endl;
    //    for (unsigned int syma=0; syma<8; syma++) {
    //        xout << "w.alphaStrings[0]" << w.alphaStrings[0].str(1) << std::endl;
    //        StringSet ka(w.alphaStrings,1,0,syma);
    //        xout << "N-1 alpha StringSet: " << ka.str(1) << std::endl;
    //    }
    //    xout << "Hamiltonian: " <<hh.str(3)<<std::endl;

    //    File ff;
    //    std::vector<double> v(3,99.0);
    //    ff.write(v,0);
    //    std::vector<double> v2(3);
    //    ff.read(v2,0);
    //    xout <<"vector read " <<v2[0]<<" "<<v2[1]<<" "<<v2[2]<<std::endl;

    State prototype(&dump);
    std::vector<gci::Hamiltonian*> hamiltonians;
    w.diagonalHamiltonian(hh);
    Hamiltonian fh = hh.FockHamiltonian(w.determinantAt(w.minloc()));
    hamiltonians.push_back(&fh);
    Hamiltonian h1(hh); h1-=fh;
    hamiltonians.push_back(&h1);
    std::vector<double> emp = gci::RSPT(hamiltonians, prototype,6);
    xout <<"MP energies" ;
    for (int i=0; i<(int)emp.size(); i++)
      xout <<" "<<emp[i];
    xout <<std::endl;
    xout <<"MP total energies" ;
    double totalEnergy=0;
    for (int i=0; i<(int)emp.size(); i++)
      xout <<" "<<(totalEnergy+=emp[i]);
    xout <<std::endl;

    gci::HamiltonianPrint(hh,prototype);
    return 0;
  }
  catch (const std::string& ex) {
    xout << "uncaught exception: " << ex << std::endl;
    throw "Error";
  }
  catch (char const* ex) {
    xout << "uncaught exception: " << ex << std::endl;
    throw "Error";
  }
  //    catch (...) {
  //        xout << "uncaught exception: "  << std::endl;
  //       throw "Error";
  //    }
}
#endif
