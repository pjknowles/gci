#include "gciRun.h"
//#ifndef MOLPRO
#include "gciMolpro.h"
//#endif
#include <iomanip>
#include <algorithm>
#include "IterativeSolver.h"

using namespace gci;



using namespace LinearAlgebra;
using ParameterVectorSet = LinearAlgebra::vectorSet<double>;
using scalar = double;

static double _lastEnergy;
static double _mu;
static double _residual_q;
static bool parallel_stringset;
struct residual {
 protected:
  const gci::Operator &m_hamiltonian;
  const bool m_subtract_Energy;
  const gci::Operator *m_Q;
 public:
  residual(const gci::Operator &hamiltonian, bool subtract_Energy, gci::Operator *Q = nullptr) : m_hamiltonian(
      hamiltonian), m_subtract_Energy(subtract_Energy), m_Q(Q) {}
  void operator()(const ParameterVectorSet &psx, ParameterVectorSet &outputs, bool append = false) const {
    for (size_t k = 0; k < psx.size(); k++) {
      const std::shared_ptr<Wavefunction> x = std::static_pointer_cast<Wavefunction>(psx[k]);
      std::shared_ptr<Wavefunction> g = std::static_pointer_cast<Wavefunction>(outputs[k]);
//        profiler->start("density");
//        SMat natorb=x->naturalOrbitals();
      //    activeHamiltonian->rotate(&natorb);
//        profiler->stop("density");
      profiler->start("Hc");
      if (not append)
        g->zero();
//        xout << "x in residual "<<x->str(2)<<std::endl;
//        xout << "g "<<g->str(2)<<std::endl;
//        xout <<"g->buffer"<<g->data()<<std::endl;
//        xout << "activeHamiltonian "<<activeHamiltonian->str(2)<<std::endl;
      // HERE!!
//          g->operatorOnWavefunction(*activeHamiltonian, *x, parallel_stringset);
//        xout << "g "<<g->str(2)<<std::endl;
//        g->zero();
      g->operatorOnWavefunction(m_hamiltonian, *x, parallel_stringset);
//        xout << "g "<<g->str(2)<<std::endl;
      profiler->stop("Hc");
//        xout << "g=Hc "<<g->str(2)<<std::endl;
      if (m_subtract_Energy) {
        double cc = x->dot(*x);
        double cg = x->dot(*g);
        _lastEnergy = cg / cc;
        double epsilon = cg / cc;
        if (m_Q != nullptr) {
//                xout << "@ m_Q in _residual"<<std::endl<<*m_Q<<std::endl;
          Wavefunction m(*g);
          m.zero();
          m.operatorOnWavefunction(*m_Q, *x);
//                xout << "m "<<m.str(2)<<std::endl;
          double cm = x->dot(m);
          double gm = g->dot(m);
          _mu = cm == 0 ? 0 : (cg * cm - cc * gm) / (cm * cm - cm * cc);
          epsilon = (cg - cm * _mu + cc * _mu * _residual_q) / (cc);
          g->axpy(-_mu, m);
//                xout << "cm="<<cm<<std::endl;
//                xout << "gm="<<gm<<std::endl;
//                xout << "mu="<<_mu<<std::endl;
//                xout << "epsilon="<<epsilon<<", cg/cc="<<cg/cc<<std::endl;
//                xout << "residual after subtracting m "<<g->str(2)<<std::endl;
          // FIXME idempotency constraint to follow
          _lastEnergy = epsilon - _mu * _residual_q;
        }
//            xout << "_lastEnergy "<<_lastEnergy<<std::endl;
        g->axpy(-_lastEnergy, x);
      }
//        xout << "final residual "<<g->str(2)<<std::endl;
    }
  }
};

struct Presidual {
 private:
  const gci::Operator &m_hamiltonian;
 public:
  std::vector<size_t> pvec;
  Presidual(const gci::Operator &hamiltonian) : m_hamiltonian(hamiltonian) {}
  void operator()(const std::vector<std::vector<double> > &Pcoeff, ParameterVectorSet &outputs) const {
    for (size_t k = 0; k < Pcoeff.size(); k++) {
//      xout << "k "<<k<<", pvec.size() "<<pvec.size()<<", Pcoeff[k].size() "<<Pcoeff[k].size()<<std::endl;
      assert(pvec.size() == Pcoeff[k].size());
      std::shared_ptr<Wavefunction> g = std::static_pointer_cast<Wavefunction>(outputs[k]);
      Wavefunction w(*g);
      w.m_sparse = true;
      for (auto i = 0; i < pvec.size(); i++)
        w.buffer_sparse.insert({pvec[i], Pcoeff[k][i]});
//      xout << "Presidual: wavefunction:" << w.buffer_sparse.size() << "\n" << w << std::endl;
      auto prof = profiler->push("HcP");
      g->operatorOnWavefunction(m_hamiltonian, w);
//      xout << "g[1708407]: " << g->at(1708407) << std::endl;
//      xout << "g[1708429]: " << g->at(1708429) << std::endl;
//        xout << "g=Hc "<<g->str(2)<<std::endl;
    }
  }
};
static std::vector<gci::Operator> _IPT_Fock;
static std::vector<double> _IPT_Epsilon;
static std::vector<double> _IPT_eta;
static std::vector<Wavefunction> _IPT_c;
static std::unique_ptr<Wavefunction> _IPT_b0m;
static std::unique_ptr<gci::Operator> _IPT_Q;
struct meanfield_residual : residual {
 public:
  meanfield_residual(const gci::Operator &hamiltonian, bool subtract_Energy, gci::Operator *Q = nullptr)
      : residual(hamiltonian, subtract_Energy, Q) {}
  void operator()(const ParameterVectorSet &psx,
                  ParameterVectorSet &outputs,
                  const std::vector<double> &shift = std::vector<double>(),
                  bool append = false) const {
    for (size_t k = 0; k < psx.size(); k++) {
      const std::shared_ptr<Wavefunction> x = std::static_pointer_cast<Wavefunction>(psx[k]);
      std::shared_ptr<Wavefunction> g = std::static_pointer_cast<Wavefunction>(outputs[k]);
      auto p = profiler->push("Mean field residual");
//        xout << "_meanfield_residual: append"<<append<<std::endl;
//        xout << "_meanfield_residual: b0m"<<_IPT_b0m->values()<<std::endl;
//        xout <<&g->buffer[0]<<" "<<&_IPT_b0m->buffer[0]<<std::endl;
//        xout << "_meanfield_residual: b0m"<<_IPT_b0m->values()<<std::endl;
//        xout << "_meanfield_residual: x"<<x->values()<<std::endl;
      if (not append)
        g->set(0);
      *g -= *_IPT_b0m;
//        xout << "_meanfield_residual: g after b0m"<<g->values()<<std::endl;
      g->operatorOnWavefunction(_IPT_Fock[0], *x, parallel_stringset);
//        xout << "_meanfield_residual: g after fock"<<g->values()<<std::endl;
      g->axpy(-_IPT_Epsilon[0], x);
      gci::Wavefunction Qc(*g);
      Qc.set(0);
      Qc.operatorOnWavefunction(*_IPT_Q, *x);
//        xout << "_meanfield_residual: Qc "<<Qc.values()<<std::endl;
      g->axpy(-_IPT_eta[0], Qc);
      g->operatorOnWavefunction(m_hamiltonian.fock(x->density(1, false, true, &_IPT_c[0]) * 2, false),
                                _IPT_c[0],
                                parallel_stringset);
      g->operatorOnWavefunction(m_hamiltonian.fock(
          _IPT_c[0].density(1, false, true, &_IPT_c[0]) * (2 * (-(*x) * _IPT_c[0])),
          false),
                                _IPT_c[0],
                                parallel_stringset);
      auto dd = _IPT_c[0] * *g;
      g->axpy(-dd, _IPT_c[0]);
//        xout << "final residual "<<g->str(2)<<std::endl;
//        xout << "_meanfield_residual: c00"<<_IPT_c[0].values()<<std::endl;

      if (false) {
        // shift away ref and Koopmans to avoid singularities
        for (int k = 0; k < 2; k++)
          g->axpy(_IPT_c[k] * *x, _IPT_c[k]);
      }

//        xout << "_meanfield_residual: g"<<g->values()<<std::endl;
    }
  }
};

struct updater {
  updater(const Wavefunction &diagonals, bool subtractDiagonal)
      : m_diagonals(diagonals), m_subtractDiagonal(subtractDiagonal) {}
 private:
  const Wavefunction &m_diagonals;
  const bool m_subtractDiagonal;
 public:
  void operator()(ParameterVectorSet &psc,
                  const ParameterVectorSet &psg,
                  std::vector<double> shift = std::vector<double>(),
                  bool append = false) const {
    std::vector<double> shifts = shift;
    for (size_t state = 0; state < psc.size(); state++) {
      if (m_subtractDiagonal)
        shifts[state] -= m_diagonals.at(m_diagonals.minloc(state + 1));
      std::shared_ptr<Wavefunction> cw = std::static_pointer_cast<Wavefunction>(psc[state]);
      std::shared_ptr<const Wavefunction> gw = std::static_pointer_cast<const Wavefunction>(psg[state]);
      if (shift[state] == 0) {
        cw->times(gw.get(), &m_diagonals);
      } else {
        shifts[state] += 2 * std::numeric_limits<scalar>::epsilon()
            * std::fmax(1, std::fabs(m_diagonals.at(m_diagonals.minloc(state + 1)))); // to guard against zero
//                xout << "initial gw  in preconditioner"<<gw->str(2)<<std::endl;
//                xout << "initial cw  in preconditioner"<<cw->str(2)<<std::endl;
//                xout << "diag  in preconditioner"<<diag->str(2)<<std::endl;
//                xout << "append "<<append<<std::endl;
        cw->divide(gw.get(), &m_diagonals, shifts[state], append, true);
//                xout << "cw after divide in preconditioner"<<cw->str(2)<<std::endl;

//            if (m_Q != nullptr) {
        if (false) { //FIXME needs reimplementing
          //FIXME this is fragile to the case that cw does not have any component in Q
          // but this has to be dealt with by providing an appropriate trial function
          // however we have the diagonals right here so we can do it.
          Wavefunction m(*cw);
          m.zero();
//                xout << "cw  in preconditioner"<<cw->str(2)<<std::endl;

//                m.operatorOnWavefunction(*m_Q,*cw); // FIXME reinstate this!

//                xout << "m  in preconditioner"<<m.str(2)<<std::endl;
          double cm = cw->dot(m);
//                if (cm==0) throw std::runtime_error("IPT wavefunction has no component in Q");
          if (std::fabs(cm) < 1e-12) {
            // generate an ion trial vector
            xout << "generating ion trial vector" << std::endl;
            Wavefunction d(m_diagonals);
            d.zero();
//                    xout << "diag"<<std::endl<<diag->str(2)<<std::endl;

//                    d.operatorOnWavefunction(*m_Q,m_diagonals); //  FIXME reinstate this!

//                    xout << "d"<<std::endl<<d.str(2)<<std::endl;
            m.set(d.minloc(state + 1), 1);
//                    xout << "m"<<std::endl<<m.str(2)<<std::endl;
            double lambda = std::sqrt(_residual_q / (1 - _residual_q));
            cw->axpy(lambda, m);
            m.axpy(lambda - 1, m);
//                    xout << "cw after initial generation"<<std::endl<<cw->str(2)<<std::endl;
//                    xout << "m after initial generation"<<std::endl<<m.str(2)<<std::endl;
            cm = cw->dot(m);
          }
          double cc = cw->dot(cw);
          double lambda = -1 + std::sqrt(_residual_q * (cc - cm) / ((1 - _residual_q) * cm));
//                xout << "cc="<<cc<<std::endl;
//                xout << "cm="<<cm<<std::endl;
          cw->axpy(lambda, m);
//                xout << "cw after updating mu constraint"<<std::endl<<cw->str(2)<<std::endl;
          cc = cw->dot(cw);
          cw->axpy(1 / std::sqrt(cc) - 1, *cw);
//                xout << "cw after renormalising"<<std::endl<<cw->str(2)<<std::endl;
        }
      }
    }
  }
};

Run::Run(std::string fcidump)
    : m_hamiltonian(Operator::construct(FCIdump(fcidump))) {
#ifdef HAVE_MPI_H
  MPI_Comm_rank(MPI_COMM_COMPUTE, &parallel_rank);
  MPI_Comm_size(MPI_COMM_COMPUTE, &parallel_size);
  xout << "Parallel run of " << parallel_size << " processes" << std::endl;
  int lendata = 0;
  if (parallel_rank == 0) {
    options = Options(FCIdump(fcidump).data());
    lendata = (int) options.data().size();
  }
  MPI_Bcast(&lendata, (int) 1, MPI_INT, 0, MPI_COMM_COMPUTE);
  auto *buf = (char *) malloc(static_cast<size_t>(lendata + 1));
  if (parallel_rank == 0) for (auto i = 0; i < lendata; i++) buf[i] = options.data()[i];
  MPI_Bcast(&buf[0], lendata, MPI_CHAR, 0, MPI_COMM_COMPUTE);
  buf[lendata] = (char) 0;
  options = Options(buf);
  free(buf);
#else
  parallel_rank=0; parallel_size=1;
  options = Options(FCIdump(fcidump).data());
#endif
//  xout << "gci::Run::options="<<options.data()<<std::endl;
//  xout << "IUHF "<< options.parameter("IUHF",std::vector<int>{0})[0]<<std::endl;
//  xout << "NELEC "<< options.parameter("NELEC",std::vector<int>{0})[0]<<std::endl;
//  xout << "FUNKY "<< options.parameter("FUNKY",std::vector<int>{999})[0]<<std::endl;
}

std::unique_ptr<Profiler> gci::profiler = nullptr;
std::vector<double> Run::run() {
  if (profiler == nullptr) profiler.reset(new Profiler("GCI"));
  _nextval_counter.reset(new sharedCounter());
  profiler->reset("GCI");
  xout << "PROGRAM * GCI (General Configuration Interaction)     Author: Peter Knowles, 2014" << std::endl;
  std::vector<double> energies;
  std::string method = options.parameter("METHOD", std::vector<std::string>(1, "")).at(0);
  if (method == "MBPT" || method == "MOLLER") method = "RSPT";
  xout << "METHOD=" << method << std::endl;

  options.addParameter("EXPLICIT1", "1"); // because Operator no longer supports embedding 1-electron in 2-electron
  parallel_stringset = options.parameter("PARALLEL_STRINGSET") != 0;
  size_t referenceLocation;
  Determinant referenceDeterminant;
  State prototype;
  { // so that w goes out of scope
    auto p = profiler->push("find reference");
    Wavefunction w(State(m_hamiltonian.m_orbitalSpaces[0],
                         options.parameter("NELEC"),
                         options.parameter("ISYM") - 1,
                         options.parameter("MS2")));
    w.diagonalOperator(m_hamiltonian);
    referenceLocation = w.minloc();
    referenceDeterminant = w.determinantAt(referenceLocation);
    xout.precision(8);
    xout << std::fixed;
    xout << "Lowest energy determinant " << referenceDeterminant << " with energy " << w.at(referenceLocation)
         << std::endl;
    prototype = State(m_hamiltonian.m_orbitalSpaces[0], w.nelec, w.symmetry, w.ms2);
  }
  if (options.parameter("EXPLICIT1") == 0 && method != "RSPT") throw std::runtime_error("EXPLICIT1 has been retired");
  //hh.constructBraKet(
//        referenceDeterminant.nelec+referenceDeterminant.ms2,
//        referenceDeterminant.nelec-referenceDeterminant.ms2
//        );

  if (method == "RSPT") {
    xout << "Rayleigh-Schroedinger perturbation theory with the Fock hamiltonian" << std::endl;
    double scs_opposite = options.parameter("SCS_OPPOSITE", std::vector<double>(1, (double) 1)).at(0);
    double scs_same = options.parameter("SCS_SAME", std::vector<double>(1, (double) 1)).at(0);
    xout << "First-order hamiltonian contains " << scs_opposite << " of opposite-spin and " << scs_same
         << " of same spin" << std::endl;
    xout << "Second-order hamiltonian contains " << 1 - scs_opposite << " of opposite-spin and " << 1 - scs_same
         << " of same spin" << std::endl;
    gci::Operator h0 = m_hamiltonian.fockOperator(referenceDeterminant);
//    xout <<"h0="<<h0<<std::endl;
    gci::Operator ssh = m_hamiltonian.sameSpinOperator(referenceDeterminant);
//    xout <<"ssh="<<ssh<<std::endl;
    gci::Operator osh = m_hamiltonian;
    osh -= ssh;
    osh -= h0; // spinUnrestricted not yet implemented
//    xout <<"osh="<<osh<<std::endl;
    gci::Operator h1 = osh * scs_opposite;
    h1 += ssh * scs_same;
//    xout <<"h1="<<h1<<std::endl;
    gci::Operator h2(m_hamiltonian); // spinUnrestricted not yet implemented
    h2 -= h1;
    h2 -= h0;
//    xout <<"h2="<<h2<<std::endl;
    std::vector<gci::Operator *> hams;
    hams.push_back(&h0);
    hams.push_back(&h1);
    if (scs_opposite != (double) 1 || scs_same != (double) 1) hams.push_back(&h2);
    std::vector<double> emp = RSPT(hams, prototype);
//    std::vector<double> emp = ISRSPT(hh, h0, prototype);
    xout << std::fixed << std::setprecision(8);
    xout << "MP energies";
    for (double i : emp) xout << " " << i;
    xout << std::endl;
    xout << "MP total energies";
    double totalEnergy = 0;
    for (double &i : emp) xout << " " << (i = totalEnergy += i);
    xout << std::endl;
    energies.resize(1);
    energies[0] = totalEnergy;
#ifdef MOLPRO
    itf::SetVariables( "ENERGY_MP", &(emp.at(1)), (unsigned int) emp.size()-1, (unsigned int) 0, "" );
#endif
  } else if (method == "IPT") {
    xout << "Ionisation perturbation theory" << std::endl;
//    Operator h0 = m_hamiltonian.fockOperator(referenceDeterminant);
    IPT(m_hamiltonian, prototype, referenceLocation);
//    xout <<std::fixed << std::setprecision(8);
//    xout <<"MP energies" ; for (int i=0; i<(int)emp.size(); i++) xout <<" "<<emp[i]; xout <<std::endl;
//    xout <<"MP total energies" ; double totalEnergy=0; for (int i=0; i<(int)emp.size(); i++) xout <<" "<<(emp[i]=totalEnergy+=emp[i]); xout <<std::endl;
//    energies.resize(1); energies[0] = totalEnergy;
  } else if (method == "ISRSPT") {
    xout << "Rayleigh-Schroedinger perturbation theory with the Fock hamiltonian" << std::endl;
    Operator h0 = m_hamiltonian.fockOperator(referenceDeterminant);
    std::vector<double> emp = ISRSPT(m_hamiltonian, h0, prototype);
    xout << std::fixed << std::setprecision(8);
    xout << "MP energies";
    for (double i : emp) xout << " " << i;
    xout << std::endl;
    xout << "MP total energies";
    double totalEnergy = 0;
    for (double &i : emp) xout << " " << (i = totalEnergy += i);
    xout << std::endl;
    energies.resize(1);
    energies[0] = totalEnergy;
#ifdef MOLPRO
    itf::SetVariables( "ENERGY_MP", &(emp.at(1)), (unsigned int) emp.size()-1, (unsigned int) 0, "" );
#endif
  } else if (method == "DAVIDSON") {
    energies = Davidson(m_hamiltonian, prototype);
  } else if (method == "CS") {
    energies = CSDavidson(m_hamiltonian, prototype);
#ifdef MOLPRO
    //    itf::SetVariables( "ENERGY_METHOD", &(emp.at(1)), (unsigned int) emp.size()-1, (unsigned int) 0, "" );
#endif
  } else if (method == "DIIS") {
    energies = DIIS(m_hamiltonian, prototype);
  } else if (method == "HAMILTONIAN")
    HamiltonianMatrixPrint(m_hamiltonian, prototype);
  else if (method == "PROFILETEST") {
    double a = 1.234;
    for (int i = 0; i < 100000000; i++) a = (a + 1 / std::sqrt(a));
    energies.resize(1);
    energies[0] = a;
  } else {
    xout << "Unknown method in GCI, " << method << std::endl;
  }

  {
    auto profile = options.parameter("PROFILER", std::vector<int>(1, -1)).at(0);
    if (profile > 1) xout << profiler->str(profile, false) << std::endl;
    xout << profiler->str(profile, true) << std::endl;
  }
  _nextval_counter.reset(nullptr);

  if (false) { // just for a test
    energies.clear();
    for (const auto &w: m_wavefunctions)
      energies.push_back(w->m_properties["ENERGY"]);
  }

  {
    auto reference_energies = options.parameter("ENERGY", std::vector<double>(0));
    double diff = 0;
    for (size_t i = 0; i < reference_energies.size() && i < energies.size(); i++)
      diff += std::fabs(energies[i] - reference_energies[i]);
    if (diff > .0000001) {
      xout << "Disagreement of calculated energies:\n";
      for (auto r : energies) xout << " " << r;
      xout << "\n with reference energies:\n";
      for (auto r : reference_energies) xout << " " << r;
      xout << std::endl;
      throw std::runtime_error("Disagreement of results with reference energies");
    }
  }

  if (options.parameter("DENSITY") > 0)
    for (size_t state = 0; state < m_wavefunctions.size(); state++) {
      m_densityMatrices.emplace_back(
          m_wavefunctions[state]->density(options.parameter("DENSITY")));
      if (options.parameter("DENSITY") == 2)
        xout << "Density . hamiltonian =" << (m_densityMatrices[state] & m_hamiltonian) << std::endl;
    }
  if (options.parameter("EXCITATIONLEVEL") > 0) {
    xout << "Excitation level analysis" << std::endl;
//      auto dm_hermitian=m_wavefunctions[0]->density(2,false,true);
//      xout << "dm_hermitian\n"<<dm_hermitian<<std::endl;
    auto dm = m_wavefunctions[0]->density(2, false, false);
    auto metric = dm.metric();
    auto metricInverse = dm.metric().inverse(1e-5);
//      xout << "metric\n" << metric<<std::endl;
//      xout << "metric.inverse(1e-5)\n" << metric.inverse(1e-5)<<std::endl;
//      xout << "metric.inverse(1e-5)\n" << metricInverse<<std::endl;
//      xout << "metric*metric.inverse(1e-5)\n" << metric*metric.inverse(1e-5)<<std::endl;

    for (const auto &w : m_wavefunctions) {
//          xout << "wavefunction:"<<w->values()<<std::endl;
//          auto td0=w->density(1,false,true);
//          xout << "td0\n"<<td0<<std::endl;
      auto td = w->density(1, false, false, m_wavefunctions[0].get());
//          xout << "td\n"<<td<<std::endl;
      SMat td1(dims_t{{td.O1().size()}, {1}});
      std::copy(td.O1().block(0).begin(), td.O1().block(0).end(), td1.block(0).begin());
//          xout << "td1\n"<<td1<<std::endl;
//          xout << "metricInverse*td1\n"<<(metricInverse*td1)<<std::endl;
      xout << "Norm of wavefunction projected to ground+1 space=" << (transpose(td1) & (metricInverse * td1))
           << std::endl;
    }
  }

  return energies;
}

Run::~Run() {
  profiler.release();
  _nextval_counter.release();
}

#ifdef MOLPRO
#include "gciMolpro.h"
using namespace itf;
#endif

std::vector<double> Run::DIIS(const Operator &ham, const State &prototype, double energyThreshold, int maxIterations) {
  std::unique_ptr<gci::Operator> residual_Q;
  profiler->start("DIIS");
  profiler->start("DIIS preamble");
//  xout << "on entry to Run::DIIS energyThreshold="<<energyThreshold<<std::endl;
  if (maxIterations < 0)
    maxIterations = options.parameter("MAXIT", std::vector<int>(1, 1000)).at(0);
  xout << "MAXIT=" << maxIterations << std::endl;
  if (energyThreshold <= (double) 0)
    energyThreshold = options.parameter("TOL", std::vector<double>(1, (double) 1e-12)).at(0);
  //  xout << "after options.parameter in Run::DIIS energyThreshold="<<energyThreshold<<std::endl;
  _residual_q = options.parameter("CHARGE", std::vector<double>{0}).at(0);
  if (_residual_q > 0) {
    xout << "q=" << _residual_q << std::endl;
    residual_Q.reset(ham.projector("Q", true));
//      xout << "Q operator" <<std::endl<<residual_Q<<std::endl;
  }
//  Operator P("P",hamiltonian,true);
//  xout << "P operator" <<std::endl<<P<<std::endl;
  Wavefunction d(prototype);
  d.diagonalOperator(ham);
  size_t reference = d.minloc();
  ParameterVectorSet gg;
  gg.push_back(std::make_shared<Wavefunction>(prototype));
  ParameterVectorSet ww;
  ww.push_back(std::make_shared<Wavefunction>(prototype));
  std::static_pointer_cast<Wavefunction>(ww.back())->set((double) 0);
  std::static_pointer_cast<Wavefunction>(ww.back())->set(reference, (double) 1);
//  double e0=d.at(reference);
  //  g -= (e0-(double)1e-10);
//      xout << "Diagonal H: " << d.str(2) << std::endl;
  updater precon(d, true);
  residual resid(ham, true, residual_Q.get());
  LinearAlgebra::DIIS<scalar> solver;
  solver.m_verbosity = options.parameter("SOLVER_VERBOSITY", std::vector<int>(1, 1)).at(0);
  solver.m_thresh = energyThreshold;
  solver.m_maxIterations = static_cast<unsigned int>(maxIterations);
  for (size_t iteration = 0; iteration < maxIterations; iteration++) {
    resid(ww, gg);
    solver.addVector(ww, gg);
    std::vector<double> shift;
    shift.push_back(0);
    precon(ww, gg, shift);
    if (solver.endIteration(ww, gg)) break;
  }
  //      xout << "Final w: "<<w.str(2)<<std::endl;
  //      xout << "Final g: "<<g.str(2)<<std::endl;
  if (_residual_q > 0) {
    return std::vector<double>{_lastEnergy, _lastEnergy + _mu};
  }
  return std::vector<double>{_lastEnergy};
}

std::vector<double> Run::Davidson(
    const Operator &ham,
    const State &prototype,
    double energyThreshold, int nState, int maxIterations) {
  auto p = profiler->push("Davidson");
//  profiler->start("Davidson preamble");
  //  xout << "on entry to Run::Davidson energyThreshold="<<energyThreshold<<std::endl;
  if (maxIterations < 0)
    maxIterations = options.parameter("MAXIT", std::vector<int>(1, 1000)).at(0);
  if (nState < 0)
    nState = options.parameter("NSTATE", std::vector<int>(1, 1)).at(0);
  if (energyThreshold <= (double) 0)
    energyThreshold = options.parameter("TOL", std::vector<double>(1, (double) 1e-12)).at(0);
  xout << "Davidson eigensolver, maximum iterations=" << maxIterations;
  if (nState > 1) xout << "; number of states=" << nState;
  xout << "; energy threshold=" << std::scientific << std::setprecision(1) << energyThreshold << std::endl;
  xout << std::fixed << std::setprecision(8);
  Wavefunction d(prototype);
  d.diagonalOperator(ham);
  updater update(d, false);
  residual resid(ham, false);
  LinearAlgebra::LinearEigensystem<scalar> solver;
  solver.m_thresh = energyThreshold;
  ParameterVectorSet gg;
  ParameterVectorSet ww;
  using Pvector = std::map<size_t, double>;
  for (int root = 0; root < nState; root++) {
    std::shared_ptr<Wavefunction> w = std::make_shared<Wavefunction>(prototype);
    ww.push_back(w);
    w->allocate_buffer();
    std::shared_ptr<Wavefunction> g = std::make_shared<Wavefunction>(prototype);
    g->allocate_buffer();
    g->settilesize(
        options.parameter("TILESIZE", std::vector<int>(1, -1)).at(0),
        options.parameter("ALPHATILESIZE", std::vector<int>(1, -1)).at(0),
        options.parameter("BETATILESIZE", std::vector<int>(1, -1)).at(0));
    gg.push_back(g);
  }
  auto initialNP =
      std::min(std::max(options.parameter("PSPACE_INITIAL", nState), nState), static_cast<int>(ww.front()->size()));
  auto NP =
      std::min(std::max(options.parameter("PSPACE", 200), nState), static_cast<int>(ww.front()->size()));
  std::vector<double> initialHPP(initialNP * initialNP, (double) 0);
  std::vector<Pvector> initialP;
  for (auto p = 0; p < initialNP; p++) {
    auto det1 = d.minloc(static_cast<size_t>(p + 1));
//    xout << "P: " << det1 << " : " << d.at(det1) << std::endl;
    initialP.emplace_back(Pvector{{det1, (double) 1}});
    Wavefunction wsparse(prototype);
    wsparse.m_sparse = true;
    Wavefunction gsparse(prototype);
    gsparse.m_sparse = true;
    wsparse.set(det1, (double) 1);
    gsparse.operatorOnWavefunction(ham, wsparse);
    for (int p1 = 0; p1 <= p; p1++) {
      auto jdet1 = initialP[p1].begin()->first;
      if (gsparse.buffer_sparse.count(jdet1))
        initialHPP[p1 + p * initialNP] = initialHPP[p + p1 * initialNP] = gsparse.buffer_sparse.at(jdet1);
      else
        initialHPP[p1 + p * initialNP] = initialHPP[p + p1 * initialNP] = 0;
    }
  }
  solver.m_verbosity = options.parameter("SOLVER_VERBOSITY", std::vector<int>(1, 1)).at(0);
  solver.m_thresh = energyThreshold;
  solver.m_maxIterations = static_cast<unsigned int>(maxIterations);
  solver.m_roots = static_cast<size_t>(nState);
  std::vector<std::vector<double> > Pcoeff;
  Presidual Presid(ham);
  solver.addP(initialP, initialHPP.data(), ww, gg, Pcoeff);
//  for (auto i = 0; i < ww.size(); i++)
//    xout << "g . w after addP " << gg[i]->dot(*ww[i]) << std::endl;
//  for (auto i = 0; i < ww.size(); i++)
//    xout << "square norm of w after addP " << ww[i]->dot(*ww[i]) << std::endl;
//  for (auto i = 0; i < gg.size(); i++)
//    xout << "square norm of g after addP " << gg[i]->dot(*gg[i]) << std::endl;
//  for (auto i=0; i<gg.size(); i++)
//    gg[i]->axpy(-solver.eigenvalues()[i],*ww[i]);
//  for (auto i=0; i<gg.size(); i++)
//    xout << "square norm of g after eigenvalue subtract "<<gg[i]->dot(*gg[i])<<std::endl;
//  for (auto i = 0; i < gg.size(); i++)
//    xout << "residual active " << gg.m_active[i] << std::endl;
//  xout << "after addP, Pcoeff="; for (const auto & p : Pcoeff) for (const auto & pp : p) xout <<" "<<pp; xout <<std::endl;
  for (const auto &pp : initialP)
    Presid.pvec.push_back(pp.begin()->first);
  //  profiler->stop("Davidson preamble");

  for (size_t iteration = 1; iteration <= maxIterations; iteration++) {
//    for (auto i = 0; i < gg.size(); i++)
//      xout << "square norm of g before Presid " << gg[i]->dot(*gg[i]) << std::endl;
    Presid(Pcoeff, gg); // augment residual with contributions from P space
//    for (auto i = 0; i < gg.size(); i++)
//      xout << "square norm of g after Presid " << gg[i]->dot(*gg[i]) << std::endl;
    ParameterVectorSet pp;
//    for (const auto &pc : Pcoeff)
    std::vector<double> shift;
//    for (auto root = 0; root < nState; root++)
//      xout << "eigenvalue " << solver.eigenvalues()[root] << std::endl;
    for (auto root = 0; root < nState; root++)
      shift.push_back(-solver.eigenvalues()[root] + 1e-10);
//    for (const auto &s : shift) xout << "shift " << s << std::endl;
//    for (auto root = 0; root < nState; root++) {
//      xout << "square norm of solution " << ww[root]->dot(*ww[root]) << std::endl;
//             xout << "ww before precon "<<ww[root]->str(2)<<std::endl;
//             xout << "gg before precon "<<gg[root]->str(2)<<std::endl;
//      xout << "g[1708407]: " << static_cast<Wavefunction*>(gg[root].get())->at(1708407) << std::endl;
//      xout << "g[1708429]: " << static_cast<Wavefunction*>(gg[root].get())->at(1708429) << std::endl;
//      xout << "w[1708407]: " << static_cast<Wavefunction*>(ww[root].get())->at(1708407) << std::endl;
//      xout << "w[1708429]: " << static_cast<Wavefunction*>(ww[root].get())->at(1708429) << std::endl;
//      xout << "update" << std::endl;
//    }
    update(ww, gg, shift, true);
//    for (auto root = 0; root < nState; root++) {
//      xout << "square norm of solution " << ww[root]->dot(*ww[root]) << std::endl;
//      xout << "g[1708407]: " << static_cast<Wavefunction*>(gg[root].get())->at(1708407) << std::endl;
//      xout << "g[1708429]: " << static_cast<Wavefunction*>(gg[root].get())->at(1708429) << std::endl;
//      xout << "w[1708407]: " << static_cast<Wavefunction*>(ww[root].get())->at(1708407) << std::endl;
//      xout << "w[1708429]: " << static_cast<Wavefunction*>(ww[root].get())->at(1708429) << std::endl;
//           xout << "ww after precon "<<ww[root]->str(2)<<std::endl;
//           xout << "gg after precon "<<gg[root]->str(2)<<std::endl;
//    }
    if (solver.endIteration(ww, gg)) break;
//    xout << "after endIteration " << std::endl;
//    for (const auto &e: solver.errors()) xout << e << std::endl;
    resid(ww, gg);
//   xout << "ww after resid "<<ww.front()->str(2)<<std::endl;
//   xout << "gg after resid "<<gg.front()->str(2)<<std::endl;
    solver.addVector(ww, gg, Pcoeff);
//    xout << "after addVector, Pcoeff="; for (const auto & p : Pcoeff) for (const auto & pp : p) xout <<" "<<pp; xout <<std::endl;
    if (iteration == 1 && NP > initialNP) { // find some more P space
      auto newP = solver.suggestP(ww, gg, NP - initialNP);
      const auto addNP = newP.size(), newNP = initialNP + addNP;
      std::vector<double> addHPP(newNP * addNP, (double) 0);
      std::vector<Pvector> addP;
      for (auto p0 = 0; p0 < addNP; p0++) {
//    xout << "P: " << det1 << " : " << d.at(det1) << std::endl;
        addP.emplace_back(Pvector{{newP[p0], (double) 1}});
        Wavefunction wsparse(prototype);
        wsparse.m_sparse = true;
        Wavefunction gsparse(prototype);
        gsparse.m_sparse = true;
        wsparse.set(newP[p0], (double) 1);
        gsparse.operatorOnWavefunction(ham, wsparse);
        for (int p1 = 0; p1 < initialNP; p1++) {
          auto jdet1 = initialP[p1].begin()->first;
          if (gsparse.buffer_sparse.count(jdet1))
            addHPP[p1 + p0 * newNP] = gsparse.buffer_sparse.at(jdet1);
        }
        for (int p1 = 0; p1 <= p0; p1++) {
          auto jdet1 = addP[p1].begin()->first;
          xout << "p0=" << p0 << ", p1=" << p1 << ", jdet1=" << jdet1 << std::endl;
          xout << addHPP.size() << " > " << p1 + initialNP + p0  * newNP << std::endl;
          xout << gsparse.buffer_sparse.count(jdet1) << std::endl;
          if (gsparse.buffer_sparse.count(jdet1))
            xout << gsparse.buffer_sparse.at(jdet1) << std::endl;
          if (gsparse.buffer_sparse.count(jdet1))
            addHPP[p1 + initialNP + p0 * newNP] = gsparse.buffer_sparse.at(jdet1);
        }
      }
      Pcoeff.resize(newNP);
      solver.addP(addP, addHPP.data(), ww, gg, Pcoeff);
      for (const auto &pp : addP)
        Presid.pvec.push_back(pp.begin()->first);
    }
  }
  for (auto root = 0; root < nState; root++) {
    m_wavefunctions.push_back(std::static_pointer_cast<Wavefunction>(ww[root]));
    m_wavefunctions.back()->m_properties["ENERGY"] = solver.eigenvalues()[root];
//      if (options.parameter("DENSITY",0)>0)
//        m_wavefunctions.back()->density = m_wavefunctions.back()->density(options.parameter("DENSITY",0));
  }
//  std::cout << "Final wavefunction\n"<<dynamic_cast<std::shared_ptr<Wavefunction> >(ww[0])->str(2)<<std::endl;
//  std::cout << "get density"<<std::endl;
//  auto dens1 = std::static_pointer_cast<Wavefunction>(ww[0])->Wavefunction::density(1);
//  xout << "density:\n"<<dens1<<std::endl;
//  dens1.FCIDump("density1.fcidump");
//  auto natorb = dynamic_cast<std::shared_ptr<Wavefunction> >(ww[0])->Wavefunction::naturalOrbitals();
//  xout << "natorb:\n"<<natorb<<std::endl;
  return solver.eigenvalues();
}

std::vector<double> Run::CSDavidson(const Operator &ham,
                                    const State &prototype,
                                    double energyThreshold, int nState, int maxIterations) {
  profiler->start("Davidson");
  profiler->start("Davidson preamble");
  // xout << "on entry to Run::Davidson energyThreshold="<<energyThreshold<<std::endl;
  if (nState < 0)
    nState = options.parameter("NSTATE", std::vector<int>(1, 1)).at(0);
  xout << "nState " << nState << std::endl;
  if (maxIterations < 0)
    maxIterations = options.parameter("MAXIT", std::vector<int>(1, 1000)).at(0);
//  xout << "MAXIT="<<maxIterations<<std::endl;
  if (energyThreshold <= (double) 0)
    energyThreshold = options.parameter("TOL", std::vector<double>(1, (double) 1e-8)).at(0);
  // xout << "after options.parameter in Run::Davidson energyThreshold="<<energyThreshold<<std::endl;
  double compressionK = options.parameter("COMPRESSIONK", std::vector<double>(1, 2)).at(0);
  int compressionL = options.parameter("COMPRESSIONL", std::vector<int>(1, 1)).at(0);
  bool compressive = compressionK != 2; // whether to use compressive sampling penalty
  if (compressive)
    xout << "Compressive sampling algorithm, k=" << compressionK << ", l=" << compressionL << ", epsilon="
         << energyThreshold << std::endl;
  Wavefunction w(prototype);
  Wavefunction g(w);
  g.diagonalOperator(ham);
  size_t reference = g.minloc();
  double e0 = g.at(reference);
  g -= (e0 - (double) 1e-10);
  std::vector<double> e;
  //  xout << "Denominators: " << g.str(2) << std::endl;
  gci::File h0file;
  h0file.name = "H0";
  g.putw(h0file);
  gci::File wfile;
  wfile.name = "Wavefunction vectors";
  gci::File gfile;
  gfile.name = "Action vectors";
  w.set((double) 0);
  w.set(reference, (double) 1);
  std::vector<double> reducedHamiltonian;
  std::vector<double> elast(static_cast<unsigned long>(nState), e0 + 1);
  profiler->stop("Davidson preamble");
  for (int n = 0; n < maxIterations; n++) {
    w.putw(wfile, n);
    g.set((double) 0);
    profiler->start("Davidson Hc");
    g.operatorOnWavefunction(ham, w);
    profiler->stop("Davidson Hc");
    g.putw(gfile, n);
    reducedHamiltonian.resize((size_t) (n + 1) * (n + 1));
    profiler->start("Davidson build rH");
    {
      for (int i = n - 1; i > -1; i--)
        for (int j = n - 1; j > -1; j--)
          reducedHamiltonian[j + i * (n + 1)] = reducedHamiltonian[j + i * n];
      for (int i = 0; i < n + 1; i++) {
        g.getw(gfile, i);
        reducedHamiltonian[i + n * (n + 1)] = g * w;
      }
      gsum(&reducedHamiltonian[n * (n + 1)], n + 1);
      for (int i = 0; i < n + 1; i++)
        reducedHamiltonian[n + i * (n + 1)] = reducedHamiltonian[i + n * (n + 1)];
    }
    profiler->stop("Davidson build rH");
    // { xout << "Reduced hamiltonian:"<<std::endl; for (int i=0; i < n+1; i++) { for (int j=0; j < n+1; j++) xout <<" "<<reducedHamiltonian[j+(n+1)*i]; xout << std::endl; } }
    std::vector<double> eigenvectors(reducedHamiltonian);
    std::vector<double> eigenvalues(n + 1);
    Diagonalize(&eigenvectors[0], &eigenvalues[0], (unsigned int) (n + 1), (unsigned int) (n + 1));
    e.resize((nState > n + 1 ? n + 1 : nState));
    e.assign(eigenvalues.begin(), eigenvalues.begin() + e.size());
    xout << "Iteration " << n << ", energies:";
    xout << std::fixed;
    xout.precision(8);
    for (int i = 0; i < (int) e.size(); i++) xout << " " << eigenvalues[i];
    xout << "; ";
    // xout << std::endl << "Eigenvectors:"<<std::endl; for (int i=0; i < nState; i++) { for (int j=0; j < n+1; j++) xout <<" "<<eigenvectors[j+(n+1)*i]; xout << std::endl; }
    int track = 0;
    double tracktest = 0;
    for (int i = 0; i < (int) e.size(); i++) {
      if (std::fabs(eigenvectors[n + 1 + i * (n + 1)]) > tracktest) {
        track = i;
        tracktest = std::fabs(eigenvectors[n + 1 + i * (n + 1)]);
      }
    }

    double energy = eigenvalues[track];
    std::vector<double> alpha;
    for (int i = 0; i <= n; i++)
      alpha.push_back(eigenvectors[i + track * (n + 1)]);
    std::vector<double> hamiltonianInverse;  // to contain (H-E)^{-1} in subspace
    for (int i = 0; i <= n; i++)
      for (int j = 0; j <= n; j++) {
        double hinvij = 0;
        for (int k = 0; k <= n; k++)
          if (k != track)
            hinvij += eigenvectors[i + k * (n + 1)] * eigenvectors[j + k * (n + 1)] / (eigenvalues[k] - energy);
//        xout << "hinv["<<i<<","<<j<<"]="<<hinvij<<std::endl;
        hamiltonianInverse.push_back(hinvij);
      }

    // determine penalty magnitude, first solving for d.alpha/d.mu
    w.set((double) 0);
    for (int i = 0; i <= n; i++) {
//      xout << "alpha "<<alpha[i]<<std::endl;
      g.getw(wfile, i);
      w.axpy(alpha[i], g);
    } // w contains current wavefunction
    double l2norm = w.norm((double) 2);
    double lknorm = w.norm(compressionK);
//    xout << "l2norm="<<l2norm<<" "<<w*w<<std::endl;
    // xout << "lknorm="<<lknorm<<std::endl;
    double factor = pow(lknorm, compressionL) * pow(l2norm, -compressionK * compressionL * (double) 0.5);
    if (compressionL * (2 - compressionK) < 0) factor = -factor;
    xout << "factor=" << factor << std::endl;
    // double Pkl=(factor -
    // ((compressionL*(2-compressionK) > 0) ? 1 : -1)) / (compressionK*compressionL);
    // xout << "Pkl from eigenvectors = " <<Pkl<<std::endl;
    // construct dP/dmu in g
    g.set((double) 0);
    g.addAbsPower(w, compressionK - 2, factor / lknorm);
    g.axpy(-factor / l2norm, w);
    // project dP/dmu onto subspace
    std::vector<double> dPdmu(n + 1);
    std::vector<double> dalphadmu(n + 1);
    for (size_t i = 0; i <= (size_t) n; i++) {
      w.getw(wfile, i);
      dPdmu[i] = g * w;
      //      xout << "dPdmu[] "<<dPdmu[i]<<std::endl;
    }
    auto d2Edmu2 = (double) 0;
    for (size_t i = 0; i <= (size_t) n; i++) {
      dalphadmu[i] = (double) 0;
      for (size_t j = 0; j <= (size_t) n; j++) {
        dalphadmu[i] -= hamiltonianInverse[j + i * (n + 1)] * dPdmu[j];
        //if (!i) xout << "dPdmu["<<j<<"]="<<dPdmu[j]<<std::endl;
      }
      d2Edmu2 -= dalphadmu[i] * dPdmu[i];
      //xout << "dalphadmum["<<i<<"]="<<dalphadmu[i]<<std::endl;
    }
    double mu = d2Edmu2 == (double) 0 ? (double) 0 : sqrt(2 * energyThreshold / d2Edmu2)
        * options.parameter("PENALTY_SCALE", std::vector<double>(1, 1)).at(0);
    // xout << "d2Edmu2="<< d2Edmu2<<", mu="<<mu<<std::endl;

    // penalised equation solver here

    profiler->start("Davidson residual");
    if (compressive) {
      g.set((double) 0);
      for (int i = 0; i <= n; i++) {
        w.getw(wfile, i);
        g.axpy(alpha[i], w);
      } // g contains the current wavefunction
      w.set((double) 0);
      w.addAbsPower(g, compressionK - 2, mu * factor / (2 * lknorm));
      w.axpy(-mu * factor / (2 * l2norm), g);
    } else // !compressive
      w.set((double) 0);
    for (int i = 0; i <= n; i++) {
      g.getw(wfile, i);
//      w += energy*alpha[i] * g;
      w.axpy(energy * alpha[i], g);
      g.getw(gfile, i);
//      w -= alpha[i] * g;
      w.axpy(-alpha[i], g);
    }
    // at this point we have the residual
    g.getw(h0file);
    //g -= (energy-e0); // Davidson
    // form update
    w.putw(h0file, 1); // save a copy
    double etruncate = options.parameter("ETRUNCATE", std::vector<double>(1, -1)).at(0);
//    xout << "energyThreshold="<<energyThreshold<<std::endl;
    if (etruncate < 0) etruncate = energyThreshold;
//    xout << "etruncate="<<etruncate<<std::endl;
    double discarded;
    double ePredicted = w.update(g, discarded, etruncate);
    // xout << "discarded="<<discarded<<std::endl;
    for (double etrunc = etruncate * .3; etrunc > 1e-50 && discarded > etruncate; etrunc *= 0.3) {
      w.getw(h0file, 1); // retrieve original
      ePredicted = w.update(g, discarded, etrunc);
      // xout << "etrunc="<<etrunc<<", discarded="<<discarded<<std::endl;
    }
    // orthogonalize to previous expansion vectors
    for (int i = 0; i <= n; i++) {
      g.getw(wfile, i);
      double factor = -(g * w) / (g * g);
      gsum(&factor, 1);
//      w += factor*g;
      w.axpy(factor, g);
    }
    profiler->stop("Davidson residual");
    double norm2 = w * w;
    gsum(&norm2, 1);

    double econv = 0;
    for (int i = 0; i < (int) e.size(); i++) {
      //if (i != track)
      econv += std::fabs(e[i] - elast[i]);
      //else {
      //econv += std::fabs(ePredicted);
      //e[i] += ePredicted;
      xout << "econv=" << econv << " (ePredicted=" << ePredicted << ", e-elast=" << e[i] - elast[i] << std::endl;
      //}
    }

    elast = e;
    // normalise
    w *= ((double) 1 / std::sqrt(norm2));
    if (norm2 > (double) 1e-30 && econv > energyThreshold) continue;

    {
      profiler->start("Histogram");
      w.set((double) 0);
      for (int i = 0; i <= n; i++) {
        g.getw(wfile, i);
        w.axpy(alpha[i], g);
      }
      w.replicate();
      double histmin = 1e-14, histmax = 1.1;
      size_t nhist = 25;
      double ratio = std::pow(histmin / histmax, 1 / ((double) nhist));
      std::vector<double> edges(nhist);
      edges[0] = histmax * ratio;
      for (size_t i = 1; i < nhist; i++)
        edges[i] = edges[i - 1] * ratio;
      std::vector<std::size_t> cumulative = w.histogram(edges);
      while (!cumulative.empty() && *cumulative.end() == *(cumulative.end() - 1)) cumulative.pop_back();
      std::vector<double> fcumulative(cumulative.size());
      for (size_t i = 0; i < nhist; i++) {
        fcumulative[i] = ((double) cumulative[i]) / (double) w.size();
        if (fcumulative[i] < 1e-8) continue;
        xout << "Histogram: " << fcumulative[i] * 100 << "% > " << edges[i] << std::endl;
        if (fcumulative[i] > 1 - 1e-8 || (edges[i] < 1e-4 && i > 0 && fcumulative[i] == fcumulative[i - 1])) break;
      }
#ifdef MOLPRO
      // put the histogram to Molpro variable space
      SetVariables ("HISTOGRAM_X",&edges[0],(unsigned int)nhist, 0, "");
      SetVariables ("HISTOGRAM_Y",&fcumulative[0],(unsigned int)nhist, 0, "");
#endif
      profiler->stop("Histogram");
    }
    break;
  }
  profiler->stop("Davidson");
  return e;
}

std::vector<double> Run::RSPT(const std::vector<Operator *> &hams,
                              const State &prototype,
                              int maxOrder,
                              double energyThreshold,
                              int maxIterations) {
  if (maxOrder < 0)
    maxOrder = options.parameter("MAXORDER", std::vector<int>(1, 1000)).at(0);
  if (maxIterations < 0)
    maxIterations = options.parameter("MAXIT", std::vector<int>(1, 1000)).at(0);
  if (energyThreshold < (double) 0)
    energyThreshold = options.parameter("TOL", std::vector<double>(1, (double) 1e-8)).at(0);
  std::vector<double> e(maxOrder + 1, (double) 0);
//  for (int k=0; k<(int)hamiltonians.size(); k++)
//    HamiltonianMatrixPrint(*hamiltonians[k],prototype);
//  return e;
  if (hams.empty()) throw std::logic_error("not enough hamiltonians");
//  for (int k=0; k<(int)hamiltonians.size(); k++) xout << "H("<<k<<"): " << *hamiltonians[k] << std::endl;
  Wavefunction w(prototype);
  xout << "RSPT wavefunction size=" << w.size() << std::endl;
  Wavefunction g(w);
  g.diagonalOperator(*hams[0]);
  size_t reference = g.minloc();
  e[0] = g.at(reference);
  g -= e[0];
  g.set(reference, (double) 1);
//  xout << "Møller-Plesset denominators: " << g.str(2) << std::endl;
  gci::File h0file;
  h0file.name = "H0";
  g.putw(h0file);
  w.set((double) 0);
  w.set(reference, (double) 1);
  gci::File wfile;
  wfile.name = "Wavefunction vectors";
  w.putw(wfile, 0);
  gci::File gfile;
  gfile.name = "Action vectors";
  for (int k = 0; k < (int) hams.size(); k++) {
    g.set((double) 0);
//    xout << "hamiltonian about to be applied to reference: "<< *hams[k] <<std::endl;
    g.operatorOnWavefunction(*hams[k], w);
//    xout << "hamiltonian on reference: " << g.str(2) << std::endl;
    g.putw(gfile, k);
  }
  int nmax = maxOrder < maxIterations ? maxOrder : maxIterations + 1;
  e.resize(nmax + 1);
  for (int n = 1; n < maxOrder && n <= maxIterations; n++) {
    // construct  |n> = -(H0-E0)^{-1} ( -sum_k^{n-1} E_{n-k} |k> + sum_{k=n-h}^{n-1} H|k>) where h is the maximum order of hamiltonian
//    xout <<std::endl<<std::endl<<"MAIN ITERATION n="<<n<<std::endl;
    g.set((double) 0);
//            xout <<std::endl<< "g after set 0: " << g.str(2) <<std::endl;
    for (int k = n; k > 0; k--) {
      w.getAll(wfile, n - k);
      if (k < (int) hams.size()) {
//            xout <<"k="<<k<< " g before H.w: " << g.str(2) <<std::endl;
        g.operatorOnWavefunction(*hams[k], w);
//            xout << "g after H.w: " << g.str(2) <<std::endl;
        if (n == k) e[n] += g.at(reference);
//        if (n == k) xout << "k, E:"<<k<<" "<<e[k]<<std::endl;
      }
//        xout << "k, E:"<<k<<" "<<e[k]<<", g before -E.w: " << g.str(2) <<std::endl;
//        xout <<"w="<<w.str(2)<<std::endl;
//      g -= w * e[k];
      g.axpy(-e[k], w);
//        xout << "k, E:"<<k<<" "<<e[k]<<", g after -E.w: " << g.str(2) <<std::endl;
    }
    {
      w = -g;
      g.getw(h0file);
//    xout <<std::endl<< "Perturbed wavefunction before precondition: " << w.str(2) <<std::endl;
      w.set(reference, (double) 0);
      w /= g;
//     xout <<std::endl<< "Perturbed wavefunction, order="<<n<<": " << w.str(2) <<std::endl;
      w.putw(wfile, n);
      for (int k = 1; k < (int) hams.size(); k++) {
        if (n + k > maxOrder) break;
        g.getw(gfile, k);
//      xout <<"gfile "<<g.str(2)<<std::endl;
//      xout <<"contribution from n="<<n<<", k="<<k<<" to E("<<n+k<<")="<<g*w<<std::endl;
        e[n + k] += g * w;
      }
      gsum(&e[n + 1], (size_t) (hams.size() - 1));
    }
    xout << "n=" << n << ", E(n+1)=" << e[n + 1] << std::endl;
    if ((e[n + 1] < 0 ? -e[n + 1] : e[n + 1]) < energyThreshold && e[n + 1] != (double) 0) {
      e.resize(n + 2);
      break;
    }
  }
  return e;
}

void Run::IPT(const gci::Operator &ham, const State &prototype, const size_t referenceLocation) {
  int maxOrder = options.parameter("MAXORDER", std::vector<int>(1, 3)).at(0);
  std::vector<double> energies;
  Wavefunction d(prototype);
  xout << "IPT wavefunction size=" << d.size() << std::endl;
  _IPT_Q = std::unique_ptr<gci::Operator>(ham.projector("Q", true));
  int continuumOrbitalSymmetry;
  size_t continuumOrbitalOffset;
  for (continuumOrbitalSymmetry = 0; continuumOrbitalSymmetry < 8; continuumOrbitalSymmetry++)
    for (continuumOrbitalOffset = 0; continuumOrbitalOffset < _IPT_Q->O1().dimension(continuumOrbitalSymmetry);
         continuumOrbitalOffset++)
      if (_IPT_Q->element(continuumOrbitalOffset,
                          continuumOrbitalSymmetry,
                          continuumOrbitalOffset,
                          continuumOrbitalSymmetry,
                          false) != 0)
        goto continuumFound;
  throw std::runtime_error("Continuum orbital cannot be found");
  continuumFound:
  xout << "IPT Q operator" << *_IPT_Q << std::endl;
  _IPT_c.clear();
  _IPT_c.emplace_back(Wavefunction(prototype));
  _IPT_c[0].set((double) 0);
  _IPT_c[0].set(referenceLocation, (double) 1);
  _IPT_Fock.clear();
  xout << "gamma00 " << _IPT_c[0].density(1, true) << std::endl;;
  _IPT_Fock.emplace_back(ham.fock(_IPT_c[0].density(1, true), true, "F00"));
  _IPT_Epsilon.clear();
  auto referenceDeterminant = _IPT_c[0].determinantAt(referenceLocation);
  d.diagonalOperator(_IPT_Fock[0]);
  _IPT_Epsilon.push_back(d.at(referenceLocation));
  _IPT_Epsilon.push_back(0);
  _IPT_c.emplace_back(Wavefunction(prototype)); // c[1]
  gci::Operator excK(_IPT_Fock[0]);
  excK.O1() *= 0;
  excK.O1(false) *= 0;
  excK.m_description = "excK";
  double io = options.parameter("IO", std::vector<double>(1, 1.1)).at(0);
  int ioo = io - 1;
  int ios = 10 * (io - ioo - 1) - 1;
  xout << "Ionise orbital " << ioo + 1 << "." << ios + 1 << std::endl;
  xout << "Continuum orbital " << continuumOrbitalOffset + 1 << "." << continuumOrbitalSymmetry + 1 << std::endl;
  _IPT_eta.clear();
  _IPT_eta.push_back(-_IPT_Fock[0].element(ioo, ios, ioo, ios)); // eta[0]
//  _IPT_eta.push_back(0); // eta[1]
  excK.element(ioo, ios, continuumOrbitalOffset, continuumOrbitalSymmetry, false) = 1;
  excK.m_description = "Excitor";
  xout << excK << std::endl;
  _IPT_c[1].set(0);
  _IPT_c[1].operatorOnWavefunction(excK, _IPT_c[0], parallel_stringset);
  _IPT_Fock.emplace_back(ham.fock(_IPT_c[0].density(1, true, true, &_IPT_c[1]), false, "F01")); // F01
  xout << "c[0]: " << _IPT_c[0].values() << std::endl;
  xout << "c[1]: " << _IPT_c[1].values() << std::endl;
  xout << "gamma01 " << _IPT_c[0].density(1, true, true, &_IPT_c[1]) << std::endl;
  xout << _IPT_Fock[0] << std::endl;
  xout << _IPT_Fock[1] << std::endl;
  {
    auto g = Wavefunction(prototype);
    g.set(0);
    g.operatorOnWavefunction(_IPT_Fock[0], _IPT_c[0]);
  }
  auto dnull = _IPT_c[0].density(1);
  dnull.O1() *= 0;
  auto h = ham.fock(dnull, true, "One-electron hamiltonian");
  xout << h << std::endl;
  energies.push_back(ham.m_O0);
  {
    Wavefunction g0(prototype);
    g0.set(0);
    g0.operatorOnWavefunction(_IPT_Fock[0], _IPT_c[0]);
    g0.operatorOnWavefunction(h, _IPT_c[0]);
    energies[0] += 0.5 * g0.dot(_IPT_c[0]);
  }
  energies.push_back(0);
//  xout << "diagonal d"<<d.str(3)<<std::endl;
  for (int m = 2; m <= maxOrder; m++) {
    xout << "Start orbital relaxation order m=" << m << std::endl;
    // construct F0m*
    _IPT_Fock.emplace_back(gci::Operator(_IPT_Fock[0]));
    _IPT_Fock.back().zero();
    _IPT_Fock.back().m_description = "F0" + std::to_string(m) + "*";
    for (int j = 1; j < m; j++) {
//                xout << "density"<<j<<m-j<<_IPT_c[j].density(1, true , true, &_IPT_c[m-j], "gamma "+std::to_string(j)+std::to_string(m-j), parallel_stringset) <<std::endl;
      _IPT_Fock.back() += ham.fock(
          _IPT_c[j].density(1, true, true, &_IPT_c[m - j], "", parallel_stringset),
          false);
    }
    xout << _IPT_Fock.back() << std::endl;
    // construct b0m
    _IPT_b0m.reset(new Wavefunction(prototype));
    _IPT_b0m->set(0);
    _IPT_b0m->operatorOnWavefunction(_IPT_Fock.back(), _IPT_c[0]);
    xout << "b0m after F0m*: " << _IPT_b0m->values() << std::endl;
    for (int k = 1; k < m; k++) {
      _IPT_b0m->operatorOnWavefunction(_IPT_Fock[k], _IPT_c[m - k]);
      xout << "b0m after F0k: " << _IPT_b0m->values() << std::endl;
      _IPT_b0m->axpy(-_IPT_Epsilon[k], _IPT_c[m - k]);
      xout << "b0m after epsilon: " << _IPT_b0m->values() << std::endl;
    }
    for (int k = 1; k < m - 1; k++) {
      gci::Wavefunction Qc(prototype);
      Qc.set(0);
      Qc.operatorOnWavefunction(*_IPT_Q, _IPT_c[m - k]);
      xout << "eta " << k << " " << _IPT_eta[k] << std::endl;
      _IPT_b0m->axpy(-_IPT_eta[k], Qc);
      xout << "b0m after eta: " << _IPT_b0m->values() << std::endl;
    }
    for (int k = 0; k < m - 1; k++)
      _IPT_b0m->axpy(_IPT_eta[k], _IPT_c[m - k - 2]);
    xout << "b0m after eta: " << _IPT_b0m->values() << std::endl;

    for (int k = 0; k < 2; k++) {
      auto b0mc0 = _IPT_c[k] * *_IPT_b0m;
      _IPT_b0m->axpy(-b0mc0, _IPT_c[k]);
    }
    xout << "b0m after project: " << _IPT_b0m->values() << std::endl;

    xout << "b0m: " << _IPT_b0m->values() << std::endl;
    xout << "solve for c0" + std::to_string(m) << std::endl;
    // solve for c0m
    ParameterVectorSet gg;
    gg.push_back(std::make_shared<Wavefunction>(prototype));
    ParameterVectorSet ww;
    ww.push_back(std::make_shared<Wavefunction>(prototype));
    updater update(d, false);
    meanfield_residual resid(ham, false);
    if (false) { // print A
      _IPT_b0m->set(0);
      for (size_t i = 0; i < ww.back()->size(); i++) {
        std::static_pointer_cast<Wavefunction>(ww.back())->set((double) 0);
        std::static_pointer_cast<Wavefunction>(ww.back())->set(i, (double) 1);
        std::static_pointer_cast<Wavefunction>(gg.back())->set((double) 0);
        std::static_pointer_cast<Wavefunction>(gg.back())->operatorOnWavefunction(_IPT_Fock[0],
                                                                                  *std::static_pointer_cast<Wavefunction>(
                                                                                      ww.back()));
        xout << "F0[" << i << ":] = " <<
             std::static_pointer_cast<Wavefunction>(gg[0])->values() << std::endl;
        resid(ww, gg);
        xout << "A[" << i << ":] = " <<
             std::static_pointer_cast<Wavefunction>(gg[0])->values() << std::endl;
      }
      exit(0);
    }
    std::static_pointer_cast<Wavefunction>(ww.back())->set((double) 0);
    std::static_pointer_cast<Wavefunction>(ww.back())->set(referenceLocation + m % 2, (double) 1);
    LinearAlgebra::DIIS<scalar> solver;
    solver.m_verbosity = options.parameter("SOLVER_VERBOSITY", std::vector<int>(1, 1)).at(0);
    solver.m_thresh = options.parameter("TOL", std::vector<double>(1, (double) 1e-8)).at(0);
    solver.m_maxIterations = options.parameter("MAXIT", std::vector<int>(1, 1000)).at(0);
    solver.m_linear = true;
//      solver.solve(gg,ww);
    for (size_t iteration = 0; iteration < solver.m_maxIterations; iteration++) {
      resid(ww, gg);
      solver.addVector(ww, gg);
      std::vector<double> shift;
      shift.push_back(0);
      update(ww, gg, shift);
      if (solver.endIteration(ww, gg)) break;
    }
    xout << "Final g: " << std::static_pointer_cast<Wavefunction>(gg[0])->values() << std::endl;
//      xout << "Final w: "<<ww[0]->str(2)<<std::endl;
    _IPT_c.push_back(*std::static_pointer_cast<Wavefunction>(ww[0]));
    xout << "c0" + std::to_string(m) << _IPT_c.back().values() << std::endl;
    // set reference component of c0m
    double refc0m = 0;
    for (int k = 1; k < m; k++)
      refc0m -= 0.5 * _IPT_c[k].dot(_IPT_c[m - k]);
    _IPT_c.back().set(referenceLocation, refc0m);
    xout << "c0m after setting reference component: " << _IPT_c[m].values() << std::endl;
    // set Koopmans component of c0m
    {
      double refc1m = 0;
      for (int k = 1; k < m + 1; k++) {
        gci::Wavefunction Qc(_IPT_c[k]);
        Qc.set(0);
        Qc.operatorOnWavefunction(*_IPT_Q, _IPT_c[m + 1 - k]);
        refc1m -= 0.5 * Qc.dot(_IPT_c[k]);
      }
      _IPT_c[m].axpy(refc1m, _IPT_c[1]);
//      _IPT_c.back += _IPT_c[1] * (refc1m - _IPT_c[m].dot(_IPT_c[1]));
      xout << "c0m after setting Koopmans component: " << _IPT_c[m].values() << std::endl;
    }
//      xout << "density"<<0<<m-0<<_IPT_c[0].density(1, true , true, &_IPT_c[m-0], "gamma "+std::to_string(0)+std::to_string(m-0), parallel_stringset) <<std::endl;


    // evaluate F0m
    _IPT_Fock.back() += ham.fock(
        _IPT_c[0].density(1, true, true, &_IPT_c[m], "", parallel_stringset),
        false) * 2;
    _IPT_Fock.back().m_description = "F0" + std::to_string(m);
    xout << _IPT_Fock.back() << std::endl;
    // evaluate E0m
    energies.push_back(0);
    for (int k = 0; k <= m; k++) {
      for (int j = 0; j <= m - k; j++) {
//              if (j%2 && (m-j)%2) {
        Wavefunction g(prototype);
        g.set(0);
        g.operatorOnWavefunction(_IPT_Fock[k], _IPT_c[m - j - k], parallel_stringset);
        if (k == 0)
          g.operatorOnWavefunction(h, _IPT_c[m - j - k], parallel_stringset);
        energies.back() += 0.5 * g.dot(_IPT_c[j]);
      }
//            }
    }
    // evaluate Epsilon0m
    _IPT_Epsilon.push_back(0);
    for (int k = 1; k <= m; k++) {
      Wavefunction g(prototype);
      g.set(0);
      g.operatorOnWavefunction(_IPT_Fock[k], _IPT_c[m - k], parallel_stringset);
      _IPT_Epsilon.back() += g.dot(_IPT_c[0]);
    }
//      xout << "Epsiloon after Fock "<<_IPT_Epsilon.back()<<std::endl;
    for (int k = 1; k < m; k++)
      _IPT_Epsilon.back() -= _IPT_c[m - k].dot(_IPT_c[0]) * _IPT_Epsilon[k];
//      xout << "Epsiloon after Epsilon "<<_IPT_Epsilon.back()<<std::endl;
    for (int k = 0; k < m - 1; k++)
      _IPT_Epsilon.back() += _IPT_c[m - k - 2].dot(_IPT_c[0]) * _IPT_eta[k];
//      xout << "eta[0] "<<_IPT_eta[0]<<std::endl;
//      xout << "Epsiloon after eta "<<_IPT_Epsilon.back()<<std::endl;
    // evaluate eta0m
    _IPT_eta.push_back(0);
    for (int k = 0; k <= m; k++) {
      Wavefunction g(prototype);
      g.set(0);
      g.operatorOnWavefunction(_IPT_Fock[k], _IPT_c[m - k], parallel_stringset);
      _IPT_eta.back() += g.dot(_IPT_c[1]);
    }
    for (int k = 0; k < m; k++)
      _IPT_eta.back() -= _IPT_c[m - k].dot(_IPT_c[1]) * _IPT_Epsilon[k];
    for (int k = 0; k < m - 1; k++)
      _IPT_eta.back() -= _IPT_c[m - k].dot(_IPT_c[1]) * _IPT_eta[k];
    for (int k = 0; k < m - 1; k++)
      _IPT_eta.back() += _IPT_c[m - k - 2].dot(_IPT_c[1]) * _IPT_eta[k];
    xout << "Energies:";
    for (auto e : energies) xout << " " << e;
    xout << std::endl;
    xout << "Energies:";
    for (auto e = energies.begin(); e != energies.end(); e++)
      xout << " " << std::accumulate(energies.begin(),
                                     e + 1,
                                     (double) 0);
    xout << std::endl;
    xout << "Energies:";
    for (auto e = energies.begin(); e != energies.end(); e++)
      xout << " " << std::accumulate(energies.begin() + 1,
                                     e + 1,
                                     (double) 0);
    xout << std::endl;
    xout << "Epsilon:";
    for (auto e : _IPT_Epsilon) xout << " " << e;
    xout << std::endl;
    xout << "eta:";
    for (auto e : _IPT_eta) xout << " " << e;
    xout << std::endl;

    if (true) { // check residual
      Wavefunction gcheck(_IPT_c[0]);
      gcheck.zero();
      for (int k = 0; k <= m; k++) {
        gcheck.operatorOnWavefunction(_IPT_Fock[k], _IPT_c[m - k]);
        gcheck.axpy(-_IPT_Epsilon[k], _IPT_c[m - k]);
        if (k < m) {
          gci::Wavefunction Qc(_IPT_c[0]);
          Qc.set(0);
          Qc.operatorOnWavefunction(*_IPT_Q, _IPT_c[m - k]);
          gcheck.axpy(-_IPT_eta[k], Qc);
        }
        if (k < m - 1)
          gcheck.axpy(_IPT_eta[k], _IPT_c[m - k - 2]);
      }
      xout << "g0" + std::to_string(m) << gcheck.values() << std::endl;
    }

    if (true) { // check idempotency
      gci::Operator idem = _IPT_c[0].density(1, true, true, &_IPT_c[0], "", parallel_stringset);
      idem.zero();
      for (int k = 0; k < m; k++) {
        idem -=
            _IPT_c[k].density(1, true, true, &_IPT_c[m - k], "", parallel_stringset);
        for (int l = 0; l < m - k; l++) {
          auto square = idem;
          square.zero();
          auto d1 = _IPT_c[k].density(1, true, true, &_IPT_c[k], "", parallel_stringset);
          auto d2 = _IPT_c[m - l - k].density(1, true, true, &_IPT_c[m], "", parallel_stringset);
          square.O1(true) = d1.O1(true) * d2.O1(true);
          square.O1(false) = d1.O1(false) * d2.O1(false);
          idem += square;
        }
      }
      xout << "idempotency " << idem << std::endl;
    }
  }
  Wavefunction neutralstate(_IPT_c[0]);
  Wavefunction ionstate(_IPT_c[1]);
  double eion = 0;
  double eneutral = 0;
  for (int m = 0; m <= maxOrder; m++) {
    if (false) {
      {
        double norm = 0;
        for (int k = 0; k <= m; k++) norm += _IPT_c[k] * _IPT_c[m - k];
        xout << "check normalisation m=" << m << " " << norm << std::endl;
      }

      {
        double norm = 0;
        for (int k = 0; k <= m; k++) {
          gci::Wavefunction Qc(_IPT_c[k]);
          Qc.set(0);
          Qc.operatorOnWavefunction(*_IPT_Q, _IPT_c[m - k]);
          norm += _IPT_c[k] * Qc;
        }
        xout << "check charge normalisation m=" << m << " " << norm << std::endl;
      }
    }
    xout << "c0" + std::to_string(m) << _IPT_c[m].values() << std::endl;
    if (m % 2 && m > 1) {
      ionstate += _IPT_c[m];
      xout << "ionstate" << ionstate.values() << std::endl;
    } else if (m > 1) {
      neutralstate += _IPT_c[m];
      xout << "neutralstate" << neutralstate.values() << std::endl;
    }
    double e = 0;
    double nn = 0;
    for (auto k = 1; k <= m; k += 2) {
      Wavefunction gggg(_IPT_c[0]);
      gggg.zero();
      gggg.operatorOnWavefunction(ham, _IPT_c[k]);
      e += gggg * _IPT_c[m - k];
      nn += _IPT_c[k] * _IPT_c[m - k];
    }
    xout << "state-summed ion energy contribution: " << e << " " << nn << std::endl;
    eion += e;
    xout << "state-summed ion energy : " << eion << std::endl;

    e = 0;
    nn = 0;
    for (auto k = 0; k <= m; k += 2) {
      Wavefunction gggg(_IPT_c[0]);
      gggg.zero();
      gggg.operatorOnWavefunction(ham, _IPT_c[k]);
      e += gggg * _IPT_c[m - k];
      nn += _IPT_c[k] * _IPT_c[m - k];
    }
    xout << "state-summed neutral energy contribution: " << e << " " << nn << std::endl;
    eneutral += e;
    xout << "state-summed neutral energy : " << eneutral << std::endl;
  }
  xout << "ionstate " << ionstate * ionstate << std::endl;
}

std::vector<double> Run::ISRSPT(
    const gci::Operator &ham,
    const gci::Operator &ham0,
    const State &prototype,
    int maxOrder,
    double energyThreshold,
    int maxIterations) {
  if (maxOrder < 0)
    maxOrder = options.parameter("MAXORDER", std::vector<int>(1, 1000)).at(0);
  if (maxIterations < 0)
    maxIterations = options.parameter("MAXIT", std::vector<int>(1, 1000)).at(0);
  if (energyThreshold < (double) 0)
    energyThreshold = options.parameter("TOL", std::vector<double>(1, (double) 1e-8)).at(0);
  std::vector<double> e(maxOrder + 1, (double) 0);
  Wavefunction d(prototype);
  xout << "RSPT wavefunction size=" << d.size() << std::endl;
  d.diagonalOperator(ham0);
  size_t reference = d.minloc();
  ParameterVectorSet gg;
  gg.push_back(std::make_shared<Wavefunction>(prototype));
  ParameterVectorSet ww;
  ww.push_back(std::make_shared<Wavefunction>(prototype));
  std::static_pointer_cast<Wavefunction>(ww.back())->set((double) 0);
  std::static_pointer_cast<Wavefunction>(ww.back())->set(reference, (double) 1);
  updater update(d, false);
  residual resid(ham, false);
  LinearAlgebra::LinearEigensystem<double> solver; // TODO
  solver.m_verbosity = options.parameter("SOLVER_VERBOSITY", std::vector<int>(1, 1)).at(0);
  solver.m_thresh = energyThreshold;
  solver.m_maxIterations = maxIterations;
//  solver.solve(gg,ww);
  for (int iteration = 0; iteration < maxIterations; iteration++) {
    resid(ww, gg);
    solver.addVector(ww, gg);
    std::vector<double> shift;
    shift.push_back(0);
    update(ww, gg, shift);
    if (solver.endIteration(ww, gg)) break;
  }
  //      xout << "Final w: "<<w.str(2)<<std::endl;
  //      xout << "Final g: "<<g.str(2)<<std::endl;
//  return solver.incremental_energies(); // TODO
  return solver.eigenvalues();
}

void Run::HamiltonianMatrixPrint(Operator &hamiltonian, const State &prototype, int verbosity) {
  Wavefunction w(&hamiltonian.m_orbitalSpaces[0], prototype.nelec, prototype.symmetry, prototype.ms2);
  Wavefunction g(w);
  xout << std::endl << "Full Hamiltonian matrix" << std::endl;
  if (verbosity >= 0) {
    for (size_t i = 0; i < w.size(); i++) {
      w.set((double) 0);
      w.set(i, (double) 1);
      g.set((double) 0);
      g.operatorOnWavefunction(hamiltonian, w);
      for (size_t j = 0; j < w.size(); j++)
        xout << (std::abs(g.at(j)) > 1e-7 ? g.at(j) : 0) << " ";
      xout << std::endl;
    }
  }
}
