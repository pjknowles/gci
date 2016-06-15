#include "gciRun.h"
#include "gciWavefunction.h"
//#ifndef MOLPRO
#include "gciMolpro.h"
//#endif
#include "FCIdump.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
using namespace gci;

Run::Run(std::string fcidump)
{
  globalFCIdump = new FCIdump(fcidump);
}

  Profiler gci::profiler("GCI");
std::vector<double> Run::run()
{
  profiler.reset("GCI");
  xout <<"PROGRAM * GCI (General Configuration Interaction)     Author: Peter Knowles, 2014" << std::endl;
#ifdef MOLPRO
   parallel_rank=(int64_t)GET_IPROCS_CXX();
   parallel_size=(int64_t)GET_NPROCS_CXX();
#elif defined(GCI_PARALLEL)
  PPIDD_Rank(&parallel_rank);
  PPIDD_Size(&parallel_size);
  xout << "Parallel run of "<<parallel_size<<" processes"<< std::endl;
#else
  parallel_rank=0;
  parallel_size=1;
#endif
  std::vector<double>energies;
  std::string method = parameter("METHOD",std::vector<std::string>(1,"")).at(0);
  xout << "method="<<method<<std::endl;
  if (method == "MBPT" || method == "MOLLER") method="RSPT";
  xout << "METHOD="<<method<<std::endl;

  profiler.start("load Hamiltonian");
  Hamiltonian hh(globalFCIdump);
  std::vector<double> rot(hh.total(0,0));
  double theta=-0*std::acos(-1.0)/10;
  for (std::vector<double>::iterator i=rot.begin(); i!=rot.end(); i++) *i=(double)0;
  for (unsigned int i=0; i<hh.basisSize; i++) rot[hh.pairIndex(i+1,i+1,0)]=(double)1;
  rot[hh.pairIndex(1,1,0)]=rot[hh.pairIndex(2,2,0)]=std::cos(theta); rot[hh.pairIndex(2,1,0)]=std::sin(theta); rot[hh.pairIndex(1,2,0)]=-std::sin(theta);
  xout << "rotation matrix"; for (std::vector<double>::const_iterator i=rot.begin(); i!=rot.end(); i++) xout <<" "<< *i; xout <<std::endl;
  hh.rotate(&rot);
  profiler.stop("load Hamiltonian");
  size_t referenceLocation;
  Determinant referenceDeterminant;
  State prototype;
  profiler.start("find reference");
  { // so that w goes out of scope
    Wavefunction w(globalFCIdump);
    w.diagonalHamiltonian(hh);
    referenceLocation = w.minloc();
    referenceDeterminant = w.determinantAt(referenceLocation);
    xout.precision(8);
    xout <<std::fixed;
    xout << "Lowest energy determinant " << referenceDeterminant <<" with energy "<<w.at(referenceLocation)<<std::endl;
    prototype = State(&hh,w.nelec,w.symmetry,w.ms2);
  }
  profiler.stop("find reference");
  if (parameter("EXPLICIT1").at(0)==0 && method != "RSPT") hh.constructBraKet(
        referenceDeterminant.nelec+referenceDeterminant.ms2,
        referenceDeterminant.nelec-referenceDeterminant.ms2
        );

  if (method == "RSPT") {
    xout << "Rayleigh-Schroedinger perturbation theory with the Fock hamiltonian" << std::endl;
    double scs_opposite = parameter("SCS_OPPOSITE",std::vector<double>(1,(double)1)).at(0);
    double scs_same = parameter("SCS_SAME",std::vector<double>(1,(double)1)).at(0);
    xout << "First-order hamiltonian contains " << scs_opposite<<" of opposite-spin and "<< scs_same <<" of same spin"<<std::endl;
    xout << "Second-order hamiltonian contains " << 1-scs_opposite<<" of opposite-spin and "<< 1-scs_same <<" of same spin"<<std::endl;
    Hamiltonian h0 = hh.FockHamiltonian(referenceDeterminant);
//    xout <<"h0.spinUnrestricted="<<h0.spinUnrestricted<<std::endl;
//    xout <<"h0="<<h0<<std::endl;
    Hamiltonian ssh = hh.sameSpinHamiltonian(referenceDeterminant);
//    xout <<"ssh.spinUnrestricted="<<ssh.spinUnrestricted<<std::endl;
//    xout <<"ssh="<<ssh<<std::endl;
    Hamiltonian osh(hh,true); osh -= ssh; osh-=h0;
//    xout <<"osh.spinUnrestricted="<<osh.spinUnrestricted<<std::endl;
//    xout <<"osh="<<osh<<std::endl;
    Hamiltonian h1 = osh*scs_opposite + ssh*scs_same;
//    xout <<"h1.spinUnrestricted="<<h1.spinUnrestricted<<std::endl;
//    xout <<"h1="<<h1<<std::endl;
    Hamiltonian h2(hh,true); // spinUnrestricted not yet optimised
    h2-=h1; h2-=h0;
//    xout <<"h2.spinUnrestricted="<<h2.spinUnrestricted<<std::endl;
//    xout <<"h2="<<h2<<std::endl;
    std::vector<gci::Hamiltonian*> hamiltonians;
    hamiltonians.push_back(&h0);
    hamiltonians.push_back(&h1);
    if (scs_opposite != (double) 1 || scs_same != (double) 1) hamiltonians.push_back(&h2);
//    xout << "hamiltonians.size()" << hamiltonians.size() << std::endl;
    std::vector<double> emp = RSPT(hamiltonians, prototype);
    xout <<std::fixed << std::setprecision(8);
    xout <<"MP energies" ; for (int i=0; i<(int)emp.size(); i++) xout <<" "<<emp[i]; xout <<std::endl;
    xout <<"MP total energies" ; double totalEnergy=0; for (int i=0; i<(int)emp.size(); i++) xout <<" "<<(emp[i]=totalEnergy+=emp[i]); xout <<std::endl;
    energies.resize(1); energies[0] = totalEnergy;
#ifdef MOLPRO
    itf::SetVariables( "ENERGY_MP", &(emp.at(1)), (unsigned int) emp.size()-1, (unsigned int) 0, "" );
#endif
  } else if (method == "DAVIDSON") {
    energies = Davidson(hh, prototype);
#ifdef MOLPRO
//    itf::SetVariables( "ENERGY_METHOD", &(emp.at(1)), (unsigned int) emp.size()-1, (unsigned int) 0, "" );
#endif
  } else if (method=="SD") {
    energies.resize(1);energies[0] = SteepestDescent(hh, prototype);
  } else if (method=="HAMILTONIAN")
     HamiltonianMatrixPrint(hh,prototype);
  else if (method=="PROFILETEST") {
    double a=1.234;
    for (int i=0; i<100000000; i++) a=(a+1/std::sqrt(a));
    energies.resize(1);energies[0]=a;
  }
  else {
    xout << "Unknown method in GCI, " << method << std::endl;
  }
  xout <<profiler.str(parameter("PROFILER",std::vector<int>(1,-1)).at(0)) <<std::endl;
  return energies;
}

void Run::addParameter(const std::string& key, const std::vector<std::string>& values, const bool echo)
{
  globalFCIdump->addParameter(key,values);
  if (echo) {
    xout << "Run::addParameter "<<key<<" = ";
    for (std::vector<std::string>::const_iterator v =values.begin(); v!=values.end(); v++)
      xout <<*v<<",";
    xout <<std::endl;
  }
}

void Run::addParameter(const std::string& key, const std::string& value, const bool echo)
{
  std::vector<std::string> values(1,value);
  addParameter(key,values,echo);
}

#ifdef MOLPRO
#include "gciMolpro.h"
using namespace itf;
#endif

#include <cmath>
double Run::SteepestDescent(const Hamiltonian &hamiltonian, const State &prototype, double energyThreshold, int maxIterations)
{
  profiler.start("SteepestDescent");
  profiler.start("SteepestDescent preamble");
//  xout << "on entry to Run::SteepestDescent energyThreshold="<<energyThreshold<<std::endl;
  if (maxIterations < 0)
    maxIterations = parameter("MAXIT",std::vector<int>(1,1000)).at(0);
  xout << "MAXIT="<<maxIterations<<std::endl;
  if (energyThreshold <= (double)0)
    energyThreshold = parameter("TOL",std::vector<double>(1,(double)1e-8)).at(0);
  xout << "after parameter in Run::SteepestDescent energyThreshold="<<energyThreshold<<std::endl;
  Wavefunction w(prototype);
  Wavefunction g(w);
  g.diagonalHamiltonian(hamiltonian);
  size_t reference = g.minloc();
  double e0=g.at(reference);
//  g -= (e0-(double)1e-10);
//    xout << "Diagonal H: " << g.str(2) << std::endl;
  gci::File h0file; h0file.name="H0"; g.put(h0file);
  gci::File wfile; wfile.name="Wavefunction vectors";
  w.set((double)0); w.set(reference, (double) 1);
  double elast=e0+1;
  double e;
  profiler.stop("SteepestDescent preamble");
  for (int n=0; n < maxIterations; n++) {
    profiler.start("SteepestDescent Hc");
    g.set((double)0); g.hamiltonianOnWavefunction(hamiltonian, w);
    profiler.stop("SteepestDescent Hc");
    e=g*w;
    g.axpy(-e,w);
    xout << "Iteration "<<n<<", energy:"<< std::fixed; xout.precision(8) ;xout << e <<"; "<<std::endl;

    bool olddistw=w.distributed; w.distributed=true;
    bool olddistg=g.distributed; g.distributed=true;
    // at this point we have the residual in g, and wavefunction in w
    w.put(wfile);
    w.get(h0file); w -= (e-(double)1e-8); // get preconditioner
    // form update
    double discarded; double ePredicted = g.update(w,discarded);
    w.get(wfile); // retrieve original wavefunction
    w+=g;
    double norm2=w*w;
    gsum(&norm2,1);
    w.distributed=olddistw;
    g.distributed=olddistg;

    double econv=std::fabs(e-elast);

    elast=e;
    // normalise
    w *= ((double)1/std::sqrt(norm2));
//    xout << "norm2="<<norm2<<", econv="<<econv<<" "<<energyThreshold<<std::endl;
    if (econv < energyThreshold) break;
  }
  profiler.stop("SteepestDescent");
  return e;
}

std::vector<double> Run::Davidson(const Hamiltonian& hamiltonian,
                                  const State &prototype,
                                  double energyThreshold, int nState, int maxIterations)
{
  profiler.start("Davidson");
  profiler.start("Davidson preamble");
  xout << "on entry to Run::Davidson energyThreshold="<<energyThreshold<<std::endl;
  if (nState < 0)
    nState = parameter("NSTATE",std::vector<int>(1,1)).at(0);
  xout << "nState "<<nState<<std::endl;
  if (maxIterations < 0)
    maxIterations = parameter("MAXIT",std::vector<int>(1,1000)).at(0);
//  xout << "MAXIT="<<maxIterations<<std::endl;
  if (energyThreshold <= (double)0)
    energyThreshold = parameter("TOL",std::vector<double>(1,(double)1e-8)).at(0);
  xout << "after parameter in Run::Davidson energyThreshold="<<energyThreshold<<std::endl;
  double compressionK = parameter("COMPRESSIONK",std::vector<double>(1,2)).at(0);
  int compressionL = parameter("COMPRESSIONL",std::vector<int>(1,1)).at(0);
  bool compressive = compressionK != 2; // whether to use compressive sampling penalty
  if (compressive)
    xout << "Compressive sampling algorithm, k="<<compressionK<<", l="<<compressionL<<", epsilon="<<energyThreshold<<std::endl;
  Wavefunction w(prototype);
  Wavefunction g(w);
  g.diagonalHamiltonian(hamiltonian);
  size_t reference = g.minloc();
  double e0=g.at(reference);
  g -= (e0-(double)1e-10);
  std::vector<double> e;
  //  xout << "Denominators: " << g.str(2) << std::endl;
  gci::File h0file; h0file.name="H0"; g.put(h0file);
  gci::File wfile; wfile.name="Wavefunction vectors";
  gci::File gfile; gfile.name="Action vectors";
  w.set((double)0); w.set(reference, (double) 1);
  std::vector<double> reducedHamiltonian;
  std::vector<double> elast(nState,e0+1);
  profiler.stop("Davidson preamble");
  for (int n=0; n < maxIterations; n++) {
//    xout <<" start of iteration "<<n<<std::endl;
    w.put(wfile,n);
//    xout << "w="<<w.str(2)<<std::endl;
    g.set((double)0);
    profiler.start("Davidson Hc");
    g.hamiltonianOnWavefunction(hamiltonian, w);
//    xout << "g="<<g.str(2)<<std::endl;
    profiler.stop("Davidson Hc");
    g.put(gfile,n);
    reducedHamiltonian.resize((size_t)(n+1)*(n+1));
    profiler.start("Davidson build rH");
    {
    for (int i=n-1; i>-1; i--)
      for (int j=n-1; j>-1; j--)
        reducedHamiltonian[j+i*(n+1)] = reducedHamiltonian[j+i*n];
    bool olddistw=w.distributed; w.distributed=true;
    bool olddistg=g.distributed; g.distributed=true;
    for (int i=0; i<n+1; i++) {
      g.get(gfile,i);
      reducedHamiltonian[i+n*(n+1)] = g * w;
    }
    gsum(&reducedHamiltonian[n*(n+1)],n+1);
    for (int i=0; i<n+1; i++)
      reducedHamiltonian[n+i*(n+1)] = reducedHamiltonian[i+n*(n+1)];
    w.distributed=olddistw; g.distributed=olddistg;
    }
    profiler.stop("Davidson build rH");
//     { xout << "Reduced hamiltonian:"<<std::endl; for (int i=0; i < n+1; i++) { for (int j=0; j < n+1; j++) xout <<" "<<reducedHamiltonian[j+(n+1)*i]; xout << std::endl; } }
    std::vector<double> eigenvectors(reducedHamiltonian);
    std::vector<double> eigenvalues(n+1);
    Diagonalize(&eigenvectors[0], &eigenvalues[0], (unsigned int)(n+1), (unsigned int)(n+1));
    e.resize((nState > n+1 ? n+1 : nState));
    e.assign(eigenvalues.begin(),eigenvalues.begin()+e.size());
    xout << "Iteration "<<n<<", energies:";
    xout << std::fixed; xout.precision(8);
    for (int i=0; i < (int)e.size(); i++) xout <<" "<<eigenvalues[i];
    xout <<"; ";
    // xout << std::endl << "Eigenvectors:"<<std::endl; for (int i=0; i < nState; i++) { for (int j=0; j < n+1; j++) xout <<" "<<eigenvectors[j+(n+1)*i]; xout << std::endl; }
    int track=0; double tracktest=0;
    for (int i=0; i < (int)e.size(); i++) {
      if (std::fabs(eigenvectors[n+1+i*(n+1)]) > tracktest) {
        track=i; tracktest=std::fabs(eigenvectors[n+1+i*(n+1)]);
      }
    }

    double energy = eigenvalues[track];
    std::vector<double> alpha;
    for (int i=0; i <= n; i++)
        alpha.push_back(eigenvectors[i+track*(n+1)]);
    std::vector<double> hamiltonianInverse;  // to contain (H-E)^{-1} in subspace
    for (int i=0; i<=n; i++)
      for (int j=0; j<=n; j++) {
        double hinvij=0;
        for (int k=0; k<=n; k++)
          if (k != track)
            hinvij += eigenvectors[i+k*(n+1)]*eigenvectors[j+k*(n+1)]/(eigenvalues[k]-energy);
//        xout << "hinv["<<i<<","<<j<<"]="<<hinvij<<std::endl;
        hamiltonianInverse.push_back(hinvij);
      }

    bool olddistw=w.distributed; w.distributed=true;
    bool olddistg=g.distributed; g.distributed=true;
    // determine penalty magnitude, first solving for d.alpha/d.mu
    w.set((double)0);
    for (int i=0; i <= n; i++) {
//      xout << "alpha "<<alpha[i]<<std::endl;
      g.get(wfile,i);
      w.axpy(alpha[i] , g);
    } // w contains current wavefunction
    double l2norm = w.norm((double)2);
    double lknorm = w.norm((double)compressionK);
//    xout << "l2norm="<<l2norm<<" "<<w*w<<std::endl;
//    xout << "lknorm="<<lknorm<<std::endl;
    double factor = pow(lknorm,compressionL) * pow(l2norm,-compressionK*compressionL*(double)0.5);
    if (compressionL*(2-compressionK) < 0) factor=-factor;
    xout << "factor="<<factor<<std::endl;
    double Pkl=(factor -
                ((compressionL*(2-compressionK) > 0) ? 1 : -1)) / (compressionK*compressionL);
    xout << "Pkl from eigenvectors = " <<Pkl<<std::endl;
    // construct dP/dmu in g
    g.set((double)0);
    g.addAbsPower(w,(double)compressionK-2,factor/lknorm);
    g.axpy(-factor/l2norm,w);
    // project dP/dmu onto subspace
    std::vector<double> dPdmu(n+1);
    std::vector<double> dalphadmu(n+1);
    for (size_t i=0; i<=(size_t)n; i++){
      w.get(wfile,i);
      dPdmu[i]=g*w;
      //      xout << "dPdmu[] "<<dPdmu[i]<<std::endl;
    }
    double d2Edmu2=(double)0;
    for (size_t i=0; i<=(size_t)n; i++){
        dalphadmu[i]=(double)0;
      for (size_t j=0; j<=(size_t)n; j++){
        dalphadmu[i]-= hamiltonianInverse[j+i*(n+1)]*dPdmu[j];
	//if (!i) xout << "dPdmu["<<j<<"]="<<dPdmu[j]<<std::endl;
      }
      d2Edmu2-=dalphadmu[i]*dPdmu[i];
      //xout << "dalphadmum["<<i<<"]="<<dalphadmu[i]<<std::endl;
    }
    double mu = d2Edmu2 == (double)0 ? (double) 0 : sqrt(2*energyThreshold/d2Edmu2) * parameter("PENALTY_SCALE",std::vector<double>(1,1)).at(0);
    xout << "d2Edmu2="<< d2Edmu2<<", mu="<<mu<<std::endl;

    // penalised equation solver here

    profiler.start("Davidson residual");
    if (compressive) {
      g.set((double)0);
      for (int i=0; i <= n; i++) {
        w.get(wfile,i);
        g.axpy(alpha[i] , w);
      } // g contains the current wavefunction
      w.set((double)0);
      w.addAbsPower(g,(double)compressionK-2,mu*factor/(2*lknorm));
      w.axpy(-mu*factor/(2*l2norm),g);
    }
    else // !compressive
      w.set((double)0);
    for (int i=0; i <= n; i++) {
      g.get(wfile,i);
//      w += energy*alpha[i] * g;
      w.axpy(energy*alpha[i] , g);
      g.get(gfile,i);
//      w -= alpha[i] * g;
      w.axpy( -alpha[i] , g);
    }
    // at this point we have the residual
    g.get(h0file);
    //g -= (energy-e0); // Davidson
    // form update
    w.put(h0file,1); // save a copy
//    xout <<"residual "<<w.str(2)<<std::endl;
    double etruncate = parameter("ETRUNCATE",std::vector<double>(1,-1)).at(0);
//    xout << "energyThreshold="<<energyThreshold<<std::endl;
    if (etruncate < 0) etruncate = energyThreshold;
//    xout << "etruncate="<<etruncate<<std::endl;
    double discarded;
    double ePredicted = w.update(g,discarded,etruncate);
//    xout << "discarded="<<discarded<<std::endl;
    for (double etrunc=etruncate*.3; etrunc > 1e-50 && discarded > etruncate; etrunc*=0.3) {
      w.get(h0file,1); // retrieve original
      ePredicted = w.update(g,discarded,etrunc);
      xout << "etrunc="<<etrunc<<", discarded="<<discarded<<std::endl;
    }
//    xout <<"new expansion vector before orthogonalisation "<<w.str(2)<<std::endl;
    // orthogonalize to previous expansion vectors
    for (int i=0; i <= n; i++) {
      g.get(wfile,i);
      double factor = -(g*w)/(g*g);
      gsum(&factor,1);
//      w += factor*g;
      w.axpy(factor,g);
    }
    profiler.stop("Davidson residual");
//    xout <<"new expansion vector before normalisation "<<w.str(2)<<std::endl;
    double norm2=w*w;
    gsum(&norm2,1);
    w.distributed=olddistw;
    g.distributed=olddistg;

    double econv=0;
    for (int i=0; i<(int)e.size(); i++) {
      //if (i != track)
	econv+=std::fabs(e[i]-elast[i]);
      //else {
	//econv += std::fabs(ePredicted);
	//e[i] += ePredicted;
	xout <<"econv="<<econv << " (ePredicted="<<ePredicted<<", e-elast="<<e[i]-elast[i]<<std::endl;
	//}
    }

    elast=e;
    // normalise
    w *= ((double)1/std::sqrt(norm2));
    xout << "norm2="<<norm2<<", econv="<<econv<<" "<<energyThreshold<<std::endl;
    if (norm2 >(double) 1e-30 && econv > energyThreshold) continue;

    { profiler.start("Histogram");
      w.set((double)0);
      bool olddist=w.distributed; w.distributed=true;
      for (int i=0; i <= n; i++) {
        g.get(wfile,i);
        w.axpy(alpha[i] , g);
      }
      w.gather();
      w.distributed=olddist;
      double histmin=1e-14,histmax=1.1;
      size_t nhist=25;
      double ratio=std::pow(histmin/histmax,1/((double)nhist));
      std::vector<double> edges(nhist);
      edges[0]=histmax*ratio;
      for (size_t i=1;i<nhist;i++)
        edges[i]=edges[i-1]*ratio;
      std::vector<std::size_t> cumulative = w.histogram(edges);
      while (cumulative.size()>0 && *cumulative.end()==*(cumulative.end()-1) ) cumulative.pop_back();
      std::vector<double> fcumulative(cumulative.size());
      for (size_t i=0;i<nhist;i++) {
        fcumulative[i]=((double)cumulative[i])/(double)w.size();
        if (fcumulative[i]<1e-8) continue;
        xout << "Histogram: "<<fcumulative[i]*100<<"% > "<<edges[i]<<std::endl;
        if (fcumulative[i]>1-1e-8 || (edges[i] < 1e-4 && i>0 && fcumulative[i]==fcumulative[i-1])) break;
      }
#ifdef MOLPRO
      // put the histogram to Molpro variable space
      SetVariables ("HISTOGRAM_X",&edges[0],(unsigned int)nhist, 0, "");
      SetVariables ("HISTOGRAM_Y",&fcumulative[0],(unsigned int)nhist, 0, "");
#endif
      profiler.stop("Histogram"); }
    break;
  }
  profiler.stop("Davidson");
    return e;
}

std::vector<double> Run::RSPT(const std::vector<gci::Hamiltonian*>& hamiltonians,
                              const State &prototype,
            int maxOrder,
                              double energyThreshold,
            int maxIterations)
{
  if (maxOrder < 0)
    maxOrder = parameter("MAXORDER",std::vector<int>(1,1000)).at(0);
  if (maxIterations < 0)
    maxIterations = parameter("MAXIT",std::vector<int>(1,1000)).at(0);
  if (energyThreshold < (double)0)
    energyThreshold = parameter("TOL",std::vector<double>(1,(double)1e-8)).at(0);
  std::vector<double> e(maxOrder+1,(double)0);
//  for (int k=0; k<(int)hamiltonians.size(); k++)
//    HamiltonianMatrixPrint(*hamiltonians[k],prototype);
//  return e;
  if (hamiltonians.size() < 1) throw "not enough hamiltonians";
//  for (int k=0; k<(int)hamiltonians.size(); k++) xout << "H("<<k<<"): " << *hamiltonians[k] << std::endl;
  Wavefunction w(prototype);
  xout <<"RSPT wavefunction size="<<w.size()<<std::endl;
  Wavefunction g(w);
  g.diagonalHamiltonian(*hamiltonians[0]);
  size_t reference = g.minloc();
  e[0]=g.at(reference);
  g-=e[0];g.set(reference,(double)1);
//  xout << "Møller-Plesset denominators: " << g.str(2) << std::endl;
  gci::File h0file; h0file.name="H0"; g.put(h0file);
  w.set((double)0); w.set(reference, (double) 1);
  gci::File wfile; wfile.name="Wavefunction vectors"; w.put(wfile,0);
  gci::File gfile; gfile.name="Action vectors";
  for (int k=0; k < (int) hamiltonians.size(); k++) {
    g.set((double)0);
//    xout << "hamiltonian about to be applied to reference: "<< *hamiltonians[k] <<std::endl;
    g.hamiltonianOnWavefunction(*hamiltonians[k],w);
//    xout << "hamiltonian on reference: " << g.str(2) << std::endl;
    g.put(gfile,k);
  }
  int nmax = maxOrder < maxIterations ? maxOrder : maxIterations+1;
  e.resize(nmax+1);
  for (int n=1; n < maxOrder && n <= maxIterations; n++) {
    // construct  |n> = -(H0-E0)^{-1} ( -sum_k^{n-1} E_{n-k} |k> + sum_{k=n-h}^{n-1} H|k>) where h is the maximum order of hamiltonian
//    xout <<std::endl<<std::endl<<"MAIN ITERATION n="<<n<<std::endl;
    g.set((double)0);
//            xout <<std::endl<< "g after set 0: " << g.str(2) <<std::endl;
    for (int k=n; k>0; k--) {
      w.getAll(wfile,n-k);
      if (k < (int) hamiltonians.size()) {
//            xout <<"k="<<k<< " g before H.w: " << g.str(2) <<std::endl;
        g.hamiltonianOnWavefunction(*hamiltonians[k], w);
//            xout << "g after H.w: " << g.str(2) <<std::endl;
        if (n == k) e[n]+=g.at(reference);
//        if (n == k) xout << "k, E:"<<k<<" "<<e[k]<<std::endl;
      }
//        xout << "k, E:"<<k<<" "<<e[k]<<", g before -E.w: " << g.str(2) <<std::endl;
//        xout <<"w="<<w.str(2)<<std::endl;
      g -= e[k] * w;
//        xout << "k, E:"<<k<<" "<<e[k]<<", g after -E.w: " << g.str(2) <<std::endl;
    }
      {
    bool olddistw=w.distributed; w.distributed=true;
    bool olddistg=g.distributed; g.distributed=true;
    w = -g;
    g.get(h0file);
//    xout <<std::endl<< "Perturbed wavefunction before precondition: " << w.str(2) <<std::endl;
    w.set(reference,(double)0);
    w /= g;
//     xout <<std::endl<< "Perturbed wavefunction, order="<<n<<": " << w.str(2) <<std::endl;
    w.put(wfile,n);
    for (int k=1; k < (int) hamiltonians.size(); k++) {
      if (n+k > maxOrder) break;
      g.get(gfile,k);
//      xout <<"gfile "<<g.str(2)<<std::endl;
//      xout <<"contribution from n="<<n<<", k="<<k<<" to E("<<n+k<<")="<<g*w<<std::endl;
      e[n+k]+=g*w;
    }
    gsum(&e[n+1],(size_t)(hamiltonians.size()-1));
    w.distributed=olddistw; g.distributed=olddistg;
  }
    xout << "n="<<n<<", E(n+1)="<<e[n+1]<<std::endl;
    if ((e[n+1] < 0 ? -e[n+1] : e[n+1]) < energyThreshold && e[n+1] != (double)0) {e.resize(n+2);break;}
  }
  return e;
}

#include <cmath>
void Run::HamiltonianMatrixPrint(Hamiltonian &hamiltonian, const State &prototype, int verbosity)
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
        xout << (std::abs(g.at(j))> 1e-7 ? g.at(j) : 0) << " ";
      xout <<std::endl;
    }
  }
}

std::vector<std::string> Run::parameter(std::string key, std::vector<std::string> def)
{
//  xout <<"string parameter request "<<key<<std::endl;
#ifdef MOLPRO
  std::string r = GetOptionS(key.c_str(),"GCI");
  if (r != std::string("")) return std::vector<std::string>(1,r);
#endif
  if (globalFCIdump != NULL) return globalFCIdump->parameter(key,def);
  return def;
}

std::vector<int> Run::parameter(std::string key, std::vector<int> def)
{
//  xout <<"integer parameter request "<<key<<std::endl;
#ifdef MOLPRO
  FORTINT r = GetOptionI(key.c_str(),"GCI");
//  xout << "GetOptionI gives "<<r<<std::endl;
  if (r != (FORTINT) -1) return std::vector<int>(1,(int) r);
#endif
  if (globalFCIdump != NULL) return globalFCIdump->parameter(key,def);
//  xout <<"dropped through"<<std::endl;
  return def;
}

std::vector<double> Run::parameter(std::string key, std::vector<double> def)
{
#ifdef MOLPRO
  FORTDBL r = GetOptionF(key.c_str(),"GCI");
  if (r != (FORTDBL) -1) return std::vector<double>(1,(double) r);
#endif
  if (globalFCIdump != NULL) return globalFCIdump->parameter(key,def);
  return def;
}
