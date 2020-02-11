#include "gciRun.h"
//#ifndef MOLPRO
#include "gciMolpro.h"
//#endif
#include <iomanip>
#include <algorithm>
#include <IterativeSolver.h>
#include <Operator.h>
#include "gciMixedOperator.h"
#include "gciMixedWavefunction.h"
#include "gciDavidson.h"

using namespace gci;


int gci::parallel_size = 1;
int gci::parallel_rank = 0;
bool gci::molpro_plugin = false;
MPI_Comm gci::molpro_plugin_intercomm = MPI_COMM_NULL;
std::unique_ptr<sharedCounter> gci::_nextval_counter = nullptr;
int64_t gci::__my_first_task = 0;
int64_t gci::__task = 0;
int64_t gci::__task_granularity = 1;

using ParameterVectorSet = std::vector<Wavefunction>;
using scalar = double;

static double _lastEnergy;
static double _mu;
static double _residual_q;
static bool parallel_stringset;
struct residual {
 protected:
  const SymmetryMatrix::Operator &m_hamiltonian;
  const bool m_subtract_Energy;
  const SymmetryMatrix::Operator *m_Q;
 public:
  residual(const SymmetryMatrix::Operator &hamiltonian, bool subtract_Energy, SymmetryMatrix::Operator *Q = nullptr) : m_hamiltonian(
      hamiltonian), m_subtract_Energy(subtract_Energy), m_Q(Q) {}
  void operator()(const ParameterVectorSet &psx, ParameterVectorSet &outputs, std::vector<bool> active, bool append = false) const {
    for (size_t k = 0; k < psx.size(); k++) {
      const Wavefunction& x = psx[k];
      Wavefunction& g = outputs[k];
//        profiler->start("density");
//        SymmetryMatrix::SMat natorb=x->naturalOrbitals();
      //    activeHamiltonian->rotate(&natorb);
//        profiler->stop("density");
      if (not append)
        g.zero();
//        xout << "x in residual "<<x->str(2)<<std::endl;
//        xout << "g "<<g->str(2)<<std::endl;
//        xout <<"g->buffer"<<g->data()<<std::endl;
//        xout << "activeHamiltonian "<<activeHamiltonian->str(2)<<std::endl;
      // HERE!!
//          g->operatorOnWavefunction(*activeHamiltonian, *x, parallel_stringset);
//        xout << "g "<<g->str(2)<<std::endl;
//        g->zero();
      if (active[k]) {
        auto prof = profiler->push("Hc");
        g.operatorOnWavefunction(m_hamiltonian, x, parallel_stringset);
//        xout << "g "<<g->str(2)<<std::endl;
//        xout << "g=Hc "<<g->str(2)<<std::endl;
      }
      if (m_subtract_Energy) {
        double cc = x.dot(x);
        double cg = x.dot(g);
        _lastEnergy = cg / cc;
        double epsilon = cg / cc;
        if (m_Q != nullptr) {
//                xout << "@ m_Q in _residual"<<std::endl<<*m_Q<<std::endl;
          Wavefunction m(g);
          m.zero();
          m.operatorOnWavefunction(*m_Q, x);
//                xout << "m "<<m.str(2)<<std::endl;
          double cm = x.dot(m);
          double gm = g.dot(m);
          _mu = cm == 0 ? 0 : (cg * cm - cc * gm) / (cm * cm - cm * cc);
          epsilon = (cg - cm * _mu + cc * _mu * _residual_q) / (cc);
          g.axpy(-_mu, m);
//                xout << "cm="<<cm<<std::endl;
//                xout << "gm="<<gm<<std::endl;
//                xout << "mu="<<_mu<<std::endl;
//                xout << "epsilon="<<epsilon<<", cg/cc="<<cg/cc<<std::endl;
//                xout << "residual after subtracting m "<<g->str(2)<<std::endl;
          // FIXME idempotency constraint to follow
          _lastEnergy = epsilon - _mu * _residual_q;
        }
//            xout << "_lastEnergy "<<_lastEnergy<<std::endl;
        g.axpy(-_lastEnergy, x);
      }
//        xout << "final residual "<<g->str(2)<<std::endl;
    }
  }
};

using Pvector = std::map<size_t, double>;
struct Presidual {
 private:
  const SymmetryMatrix::Operator &m_hamiltonian;
  const std::vector<Pvector> &m_P;
 public:
  Presidual(const SymmetryMatrix::Operator &hamiltonian, const std::vector<Pvector> &P) : m_hamiltonian(hamiltonian), m_P(P) {}
  void operator()(const std::vector<std::vector<double> > &Pcoeff, ParameterVectorSet &outputs) const {
    for (size_t k = 0; k < Pcoeff.size(); k++) {
      assert(m_P.size() == Pcoeff[k].size());
      Wavefunction& g = outputs[k];
      Wavefunction w(g);
      w.m_sparse = true;
      for (size_t i = 0; i < m_P.size(); i++)
        w.buffer_sparse.insert({m_P[i].begin()->first, Pcoeff[k][i]});
      auto prof = profiler->push("HcP");
      g.operatorOnWavefunction(m_hamiltonian, w);
    }
  }
};

static std::vector<SymmetryMatrix::Operator> _IPT_Fock;
static std::vector<double> _IPT_Epsilon;
static std::vector<double> _IPT_eta;
static std::vector<Wavefunction> _IPT_c;
static std::unique_ptr<Wavefunction> _IPT_b0m;
static std::unique_ptr<SymmetryMatrix::Operator> _IPT_Q;
struct meanfield_residual : residual {
 public:
  meanfield_residual(const SymmetryMatrix::Operator &hamiltonian, bool subtract_Energy, SymmetryMatrix::Operator *Q = nullptr)
      : residual(hamiltonian, subtract_Energy, Q) {}
  void operator()(const ParameterVectorSet &psx,
                  ParameterVectorSet &outputs,
                  const std::vector<double> &shift = std::vector<double>(),
                  bool append = false) const {
    for (size_t k = 0; k < psx.size(); k++) {
      const Wavefunction& x = psx[k];
      Wavefunction& g = outputs[k];
      auto p = profiler->push("Mean field residual");
//        xout << "_meanfield_residual: append"<<append<<std::endl;
//        xout << "_meanfield_residual: b0m"<<_IPT_b0m->values()<<std::endl;
//        xout <<&g->buffer[0]<<" "<<&_IPT_b0m->buffer[0]<<std::endl;
//        xout << "_meanfield_residual: b0m"<<_IPT_b0m->values()<<std::endl;
//        xout << "_meanfield_residual: x"<<x->values()<<std::endl;
      if (not append)
        g.set(0);
      g -= *_IPT_b0m;
//        xout << "_meanfield_residual: g after b0m"<<g->values()<<std::endl;
      g.operatorOnWavefunction(_IPT_Fock[0], x, parallel_stringset);
//        xout << "_meanfield_residual: g after fock"<<g->values()<<std::endl;
      g.axpy(-_IPT_Epsilon[0], x);
      gci::Wavefunction Qc(g);
      Qc.set(0);
      Qc.operatorOnWavefunction(*_IPT_Q, x);
//        xout << "_meanfield_residual: Qc "<<Qc.values()<<std::endl;
      g.axpy(-_IPT_eta[0], Qc);
      g.operatorOnWavefunction(m_hamiltonian.fock(x.density(1, false, true, &_IPT_c[0]) * 2, false),
                                _IPT_c[0],
                                parallel_stringset);
      g.operatorOnWavefunction(m_hamiltonian.fock(
          _IPT_c[0].density(1, false, true, &_IPT_c[0]) * (2 * (-(x * _IPT_c[0]))),
          false),
                                _IPT_c[0],
                                parallel_stringset);
      auto dd = _IPT_c[0] * g;
      g.axpy(-dd, _IPT_c[0]);
//        xout << "final residual "<<g->str(2)<<std::endl;
//        xout << "_meanfield_residual: c00"<<_IPT_c[0].values()<<std::endl;

      if (false) {
        // shift away ref and Koopmans to avoid singularities
        for (int k = 0; k < 2; k++)
          g.axpy(_IPT_c[k] * x, _IPT_c[k]);
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
      Wavefunction& cw = psc[state];
      const Wavefunction& gw = psg[state];
      if (shift[state] == 0) {
        cw.times(&gw, &m_diagonals);
      } else {
        shifts[state] += 2 * std::numeric_limits<scalar>::epsilon()
            * std::fmax(1, std::fabs(m_diagonals.at(m_diagonals.minloc(state + 1)))); // to guard against zero
//                xout << "initial gw  in preconditioner"<<gw->str(2)<<std::endl;
//                xout << "initial cw  in preconditioner"<<cw->str(2)<<std::endl;
//                xout << "diag  in preconditioner"<<diag->str(2)<<std::endl;
//                xout << "append "<<append<<std::endl;
        cw.divide(&gw, &m_diagonals, shifts[state], append, true);
//                xout << "cw after divide in preconditioner"<<cw->str(2)<<std::endl;

//            if (m_Q != nullptr) {
        if (false) { //FIXME needs reimplementing
          //FIXME this is fragile to the case that cw does not have any component in Q
          // but this has to be dealt with by providing an appropriate trial function
          // however we have the diagonals right here so we can do it.
          Wavefunction m(cw);
          m.zero();
//                xout << "cw  in preconditioner"<<cw->str(2)<<std::endl;

//                m.operatorOnWavefunction(*m_Q,*cw); // FIXME reinstate this!

//                xout << "m  in preconditioner"<<m.str(2)<<std::endl;
          double cm = cw.dot(m);
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
            cw.axpy(lambda, m);
            m.axpy(lambda - 1, m);
//                    xout << "cw after initial generation"<<std::endl<<cw->str(2)<<std::endl;
//                    xout << "m after initial generation"<<std::endl<<m.str(2)<<std::endl;
            cm = cw.dot(m);
          }
          double cc = cw.dot(cw);
          double lambda = -1 + std::sqrt(_residual_q * (cc - cm) / ((1 - _residual_q) * cm));
//                xout << "cc="<<cc<<std::endl;
//                xout << "cm="<<cm<<std::endl;
          cw.axpy(lambda, m);
//                xout << "cw after updating mu constraint"<<std::endl<<cw->str(2)<<std::endl;
          cc = cw.dot(cw);
          cw.axpy(1 / std::sqrt(cc) - 1, cw);
//                xout << "cw after renormalising"<<std::endl<<cw->str(2)<<std::endl;
        }
      }
    }
  }
};

Run::Run(std::string fcidump)
    : m_hamiltonian(constructOperator(FCIdump(fcidump))) {
#ifdef HAVE_MPI_H
    MPI_Comm_rank(MPI_COMM_COMPUTE, &parallel_rank);
    MPI_Comm_size(MPI_COMM_COMPUTE, &parallel_size);
    xout << "Parallel run of " << parallel_size << " processes" << std::endl;
    int lendata = 0;
    if (parallel_rank == 0) {
        options = Options(FCIdump(fcidump).data());
        options.addParameter("FCIDUMP", fcidump);
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
  int activeLvl = options.parameter("PROFACTIVE",-1);
  if (activeLvl >= 0)  profiler->active(activeLvl);
  xout << "PROGRAM * GCI (General Configuration Interaction)     Author: Peter Knowles, 2014" << std::endl;
  std::vector<double> energies;
  std::string method = options.parameter("METHOD", std::vector<std::string>(1, "DAVIDSON")).at(0);
  if (method == "MBPT" || method == "MOLLER") method = "RSPT";
  xout << "METHOD=" << method << std::endl;

  options.addParameter("EXPLICIT1", "1"); // because Operator no longer supports embedding 1-electron in 2-electron
  parallel_stringset = options.parameter("PARALLEL_STRINGSET") != 0;
  size_t referenceLocation;
  Determinant referenceDeterminant;
  State prototype;
  { // so that w goes out of scope
    auto p = profiler->push("find reference");
    Wavefunction w(State(m_hamiltonian,
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
    prototype = State(m_hamiltonian, w.nelec, w.symmetry, w.ms2);
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
    SymmetryMatrix::Operator h0 = fockOperator(m_hamiltonian,referenceDeterminant);
//    xout <<"h0="<<h0<<std::endl;
    SymmetryMatrix::Operator ssh = sameSpinOperator(m_hamiltonian,referenceDeterminant);
//    xout <<"ssh="<<ssh<<std::endl;
    SymmetryMatrix::Operator osh = m_hamiltonian;
    osh -= ssh;
    osh -= h0; // spinUnrestricted not yet implemented
//    xout <<"osh="<<osh<<std::endl;
    SymmetryMatrix::Operator h1 = osh * scs_opposite;
    h1 += ssh * scs_same;
//    xout <<"h1="<<h1<<std::endl;
    SymmetryMatrix::Operator h2(m_hamiltonian); // spinUnrestricted not yet implemented
    h2 -= h1;
    h2 -= h0;
//    xout <<"h2="<<h2<<std::endl;
    std::vector<SymmetryMatrix::Operator *> hams;
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
    SymmetryMatrix::Operator h0 = fockOperator(m_hamiltonian,referenceDeterminant);
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
        if (options.parameter("VIBRONIC", 0)) {
            if (options.parameter("SECOND_QUANT", 0)) {
                auto ham = MixedOperatorSecondQuant(FCIdump(options.parameter("FCIDUMP", "fcidump")));
                auto wfn = MixedWavefunction(options);
                run::Davidson<MixedWavefunction, MixedOperatorSecondQuant> solver(std::move(wfn), std::move(ham),
                                                                                  options);
                solver.run();
            } else {
                auto ham = MixedOperator(FCIdump(options.parameter("FCIDUMP", "fcidump")));
                auto wfn = MixedWavefunction(options);
                run::Davidson<MixedWavefunction, MixedOperator> solver(std::move(wfn), std::move(ham), options);
                solver.run();
//        auto ham2 = m_hamiltonian;
//        auto wfn2 = Wavefunction(prototype);
//        auto solver2 = run::Davidson<Wavefunction, SymmetryMatrix::Operator>(std::move(wfn2), std::move(ham2), options);
//        solver2.run();
            }
        } else {
            auto ham = m_hamiltonian;
            auto wfn = Wavefunction(prototype);
            auto solver = run::Davidson<Wavefunction, SymmetryMatrix::Operator>(std::move(wfn), std::move(ham),
                                                                                options);
            solver.run();
//        energies = Davidson(m_hamiltonian, prototype);
        }
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
      SymmetryMatrix::SMat td1(SymmetryMatrix::dims_t{{td.O1().size()}, {1}});
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
  std::cout << *profiler.get() << std::endl;
  profiler.release();
  _nextval_counter.release();
}

namespace gci {
/*!
 * \brief Construct an operator templated on source, but with a special specification
 * \param special
 * - "Q" projector onto space containing satellite orbital
 * - "P" 1-Q
 * \param forceSpinUnrestricted whether to force conversion to a UHF object
 */
static
SymmetryMatrix::Operator* projector(const SymmetryMatrix::Operator& source,
                                    std::string special,
                                    bool forceSpinUnrestricted);
}

#ifdef MOLPRO
#include "gciMolpro.h"
using namespace itf;
#endif

std::vector<double> Run::DIIS(const SymmetryMatrix::Operator &ham, const State &prototype, double energyThreshold, int maxIterations) {
  std::unique_ptr<SymmetryMatrix::Operator> residual_Q;
  profiler->start("DIIS");
  profiler->start("DIIS preamble");
//  xout << "on entry to Run::DIIS energyThreshold="<<energyThreshold<<std::endl;
  if (maxIterations < 0)
    maxIterations = options.parameter("MAXIT", 1000);
  xout << "MAXIT=" << maxIterations << std::endl;
  if (energyThreshold <= (double) 0)
    energyThreshold = options.parameter("TOL", 1e-8);
  //  xout << "after options.parameter in Run::DIIS energyThreshold="<<energyThreshold<<std::endl;
  _residual_q = options.parameter("CHARGE", 0.0);
  if (_residual_q > 0) {
    xout << "q=" << _residual_q << std::endl;
    residual_Q.reset(projector(ham,"Q", true));
//      xout << "Q operator" <<std::endl<<residual_Q<<std::endl;
  }
//  Operator P("P",hamiltonian,true);
//  xout << "P operator" <<std::endl<<P<<std::endl;
  Wavefunction d(prototype);
  d.diagonalOperator(ham);
  size_t reference = d.minloc();
  ParameterVectorSet gg;
  gg.emplace_back(prototype);
  (gg.back()).allocate_buffer();
  ParameterVectorSet ww;
  ww.emplace_back(prototype);
  (ww.back()).set((double) 0);
  (ww.back()).set(reference, (double) 1);
//  double e0=d.at(reference);
  //  g -= (e0-(double)1e-10);
//      xout << "Diagonal H: " << d.str(2) << std::endl;
  updater precon(d, true);
  residual resid(ham, true, residual_Q.get());
  IterativeSolver::DIIS<Wavefunction> solver;
  solver.m_verbosity = options.parameter("SOLVER_VERBOSITY", 1);
  solver.m_thresh = energyThreshold;
  solver.m_maxIterations = static_cast<unsigned int>(maxIterations);
  for (size_t iteration = 0; iteration < static_cast<size_t>(maxIterations); iteration++) {
    resid(ww, gg, solver.active());
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
    const SymmetryMatrix::Operator &ham,
    const State &prototype,
    double energyThreshold, int nState, int maxIterations) {
  auto p = profiler->push("Davidson");
//  profiler->start("Davidson preamble");
  //  xout << "on entry to Run::Davidson energyThreshold="<<energyThreshold<<std::endl;
  if (maxIterations < 0)
    maxIterations = options.parameter("MAXIT", 1000);
  if (nState < 0)
    nState = options.parameter("NSTATE", 1);
  if (energyThreshold <= (double) 0)
    energyThreshold = options.parameter("TOL", 1e-8);
  xout << "Davidson eigensolver, maximum iterations=" << maxIterations;
  if (nState > 1) xout << "; number of states=" << nState;
  xout << "; energy threshold=" << std::scientific << std::setprecision(1) << energyThreshold << std::endl;
  xout << std::fixed << std::setprecision(8);
  Wavefunction d(prototype);
  d.diagonalOperator(ham);
  updater update(d, false);
  residual resid(ham, false);
  IterativeSolver::LinearEigensystem<Wavefunction> solver;
  solver.m_thresh = energyThreshold;
  ParameterVectorSet gg;
  ParameterVectorSet ww;
  for (int root = 0; root < nState; root++) {
    ww.emplace_back(prototype);
    ww.back().allocate_buffer();
    gg.emplace_back(prototype);
    gg.back().allocate_buffer();
    gg.back().settilesize(
        options.parameter("TILESIZE", std::vector<int>(1, 128)).at(0),
        options.parameter("ALPHATILESIZE", std::vector<int>(1, -1)).at(0),
        options.parameter("BETATILESIZE", std::vector<int>(1, -1)).at(0));
  }
  solver.m_verbosity = options.parameter("SOLVER_VERBOSITY", 1);
  solver.m_thresh = energyThreshold;
  solver.m_maxIterations = static_cast<unsigned int>(maxIterations);
  solver.m_roots = static_cast<size_t>(nState);

  std::vector<std::vector<double> > Pcoeff(nState);
  std::vector<Pvector> P;
  Presidual Presid(ham, P);

  auto initialNP =
      std::min(std::max(options.parameter("PSPACE_INITIAL", nState), nState), static_cast<int>(ww.front().size()));
  auto maxNP =
      std::min(std::max(options.parameter("PSPACE", 100), nState), static_cast<int>(ww.front().size()));

  for (size_t iteration = 1; iteration <= static_cast<size_t>(maxIterations); iteration++) {
    if (options.parameter("PSPACE_REBUILD", 0)) { // clear out the P space and rebuild it
      solver.clearP();
      P.clear();
      Pcoeff.clear();
      Pcoeff.resize(nState);
    }
    if (static_cast<size_t>(maxNP) > P.size()) { // find some more P space
      size_t NP = P.size();
      if (iteration == 1) { // initial
        for (auto p = 0; p < initialNP; p++) {
          auto det1 = d.minloc(static_cast<size_t>(p + 1));
          P.emplace_back(Pvector{{det1, 1}});
        }
      } else { // subsequent
        Presid(Pcoeff, gg); // augment residual with contributions from P space
        auto newP = solver.suggestP(ww, gg, (maxNP - NP));
        for (const auto &pp : P)
          newP.erase(std::remove(newP.begin(), newP.end(), pp.begin()->first),
                     newP.end()); // remove anything already in P
        for (const auto &det1 : newP)
          P.emplace_back(Pvector{{det1, 1}});
      }
      const auto newNP = P.size();
      if (solver.m_verbosity > 1 && newNP > NP)
        xout << "Adding " << newNP - NP << " P-space configurations (total " << newNP << ")" << std::endl;
      std::vector<double> addHPP(newNP * (newNP - NP), (double) 0);
      for (size_t p0 = NP; p0 < newNP; p0++) {
        Wavefunction wsparse(prototype);
        wsparse.m_sparse = true;
        Wavefunction gsparse(prototype);
        gsparse.m_sparse = true;
        wsparse.set(P[p0].begin()->first, (double) 1);
        gsparse.operatorOnWavefunction(ham, wsparse);
        for (size_t p1 = 0; p1 < newNP; p1++) {
          auto jdet1 = P[p1].begin()->first;
          if (gsparse.buffer_sparse.count(jdet1))
            addHPP[p1 + (p0 - NP) * newNP] = gsparse.buffer_sparse.at(jdet1);
        }
      }
      solver.addP(std::vector<Pvector>(P.begin() + NP, P.end()), addHPP.data(), ww, gg, Pcoeff);
    }
    Presid(Pcoeff, gg); // augment residual with contributions from P space
    std::vector<double> shift;
    for (auto root = 0; root < nState; root++)
      shift.push_back(-solver.eigenvalues()[root] + 1e-10);
    update(ww, gg, shift, true);
    if (solver.endIteration(ww, gg)) break;
    resid(ww, gg, solver.active());
    solver.addVector(ww, gg, Pcoeff);
  }
  if (solver.m_verbosity > 0)
    xout << "Number of actions of matrix on vector = " << solver.actions() << std::endl;
  for (auto root = 0; root < nState; root++) {
    m_wavefunctions.push_back(std::make_shared<Wavefunction>(ww[root]));
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

std::vector<double> Run::CSDavidson(const SymmetryMatrix::Operator &ham,
                                    const State &prototype,
                                    double energyThreshold, int nState, int maxIterations) {
  profiler->start("Davidson");
  profiler->start("Davidson preamble");
  // xout << "on entry to Run::Davidson energyThreshold="<<energyThreshold<<std::endl;
  if (nState < 0)
    nState = options.parameter("NSTATE", 1);
  xout << "nState " << nState << std::endl;
  if (maxIterations < 0)
    maxIterations = options.parameter("MAXIT", 1000);
//  xout << "MAXIT="<<maxIterations<<std::endl;
  if (energyThreshold <= (double) 0)
    energyThreshold = options.parameter("TOL", 1e-8);
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

std::vector<double> Run::RSPT(const std::vector<SymmetryMatrix::Operator *> &hams,
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
      g.getAll(h0file);
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

void Run::IPT(const SymmetryMatrix::Operator &ham, const State &prototype, const size_t referenceLocation) {
  int maxOrder = options.parameter("MAXORDER", std::vector<int>(1, 3)).at(0);
  std::vector<double> energies;
  Wavefunction d(prototype);
  xout << "IPT wavefunction size=" << d.size() << std::endl;
  _IPT_Q = std::unique_ptr<SymmetryMatrix::Operator>(projector(ham,"Q", true));
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
  SymmetryMatrix::Operator excK(_IPT_Fock[0]);
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
    _IPT_Fock.emplace_back(SymmetryMatrix::Operator(_IPT_Fock[0]));
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
    gg.emplace_back(prototype);
    ParameterVectorSet ww;
    ww.emplace_back(prototype);
    updater update(d, false);
    meanfield_residual resid(ham, false);
    if (false) { // print A not sure what this is all about!
      _IPT_b0m->set(0);
      for (size_t i = 0; i < ww.back().size(); i++) {
        ww.back().set((double) 0);
        ww.back().set(i, (double) 1);
        gg.back().set((double) 0);
        gg.back().operatorOnWavefunction(_IPT_Fock[0], ww.back());
        xout << "F0[" << i << ":] = " <<
             gg[0].values() << std::endl;
        resid(ww, gg);
        xout << "A[" << i << ":] = " <<
             gg[0].values() << std::endl;
      }
      exit(0);
    }
    ww.back().set((double) 0);
    ww.back().set(referenceLocation + m % 2, (double) 1);
    IterativeSolver::DIIS<Wavefunction> solver;
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
    xout << "Final g: " << gg[0].values() << std::endl;
//      xout << "Final w: "<<ww[0]->str(2)<<std::endl;
    _IPT_c.push_back(ww[0]);
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
      SymmetryMatrix::Operator idem = _IPT_c[0].density(1, true, true, &_IPT_c[0], "", parallel_stringset);
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
    const SymmetryMatrix::Operator &ham,
    const SymmetryMatrix::Operator &ham0,
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
  gg.emplace_back(prototype);
  ParameterVectorSet ww;
  ww.emplace_back(prototype);
  ww.back().set((double) 0);
  ww.back().set(reference, (double) 1);
  updater update(d, false);
  residual resid(ham, false);
  IterativeSolver::LinearEigensystem<Wavefunction> solver; // TODO
  solver.m_verbosity = options.parameter("SOLVER_VERBOSITY", std::vector<int>(1, 1)).at(0);
  solver.m_thresh = energyThreshold;
  solver.m_maxIterations = maxIterations;
//  solver.solve(gg,ww);
  for (int iteration = 0; iteration < maxIterations; iteration++) {
    resid(ww, gg, solver.active());
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

void Run::HamiltonianMatrixPrint(SymmetryMatrix::Operator &hamiltonian, const State &prototype, int verbosity) {
  Wavefunction w(hamiltonian, prototype.nelec, prototype.symmetry, prototype.ms2);
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

SymmetryMatrix::Operator* gci::projector(const SymmetryMatrix::Operator& source, const std::string special, const bool forceSpinUnrestricted) {
  SymmetryMatrix::dim_t dims;
  for (auto s=0; s<8; s++)
    dims.push_back(source.dimension(s,0,0));
  auto result = new SymmetryMatrix::Operator(dims,
                                             1,
                                             source.m_uhf > 0 || forceSpinUnrestricted,
                                             0,
                                             true,
                                             false,
                                             special + " projector");
  result->m_O0 = 0;

  if (special == "P" || special == "Q") {
    if (special == "P")
      result->O1(true).setIdentity();
    else
      result->O1(true).assign(0);
    if (source.m_uhf)
      result->O1(false) = result->O1(true);
    // determine the uncoupled orbital
    size_t uncoupled_orbital = 0;
    unsigned int uncoupled_orbital_symmetry = 0;
    double min_rowsum = 1e50;
    for (unsigned int sym = 0; sym < 8; sym++) {
      for (size_t k = 0; k < source.dimension(sym); k++) {
        double rowsum = 0;
        for (size_t l = 0; l < source.dimension(sym); l++)
          rowsum += source.O1(true).block(sym)[k > l ? k * (k + 1) / 2 + l : l * (l + 1) / 2 + k];
        if (std::fabs(rowsum) < min_rowsum) {
          min_rowsum = std::fabs(rowsum);
          uncoupled_orbital = k;
          uncoupled_orbital_symmetry = sym;
        }
      }
    }
    xout << "non-interacting orbital is " << uncoupled_orbital + 1 << "." << uncoupled_orbital_symmetry + 1
         << std::endl;
    result->O1(false).block(uncoupled_orbital_symmetry)[(uncoupled_orbital + 2) * (uncoupled_orbital + 1) / 2 - 1] =
        (special == "P" ? 0 : 1);
  }
  return result;
}

Eigen::VectorXd gci::int1(const SymmetryMatrix::Operator& hamiltonian, int spin) {
  size_t basisSize = 0;
  for (auto si = 0; si < 8; si++)
    basisSize += hamiltonian.O1(spin>0).dimension(si);
  Eigen::VectorXd result(basisSize);
  size_t off=0;
  for (auto si=0; si<8; si++)
    for (size_t oi=0; oi<hamiltonian.dimension(si); oi++)
      result[off++] = hamiltonian.O1(spin>0).block(si)[(oi+1)*(oi+2)/2-1];
  return result;
}

Eigen::MatrixXd gci::intJ(const SymmetryMatrix::Operator &hamiltonian, int spini, int spinj) {
    if (spinj > spini) return intJ(hamiltonian, spinj, spini).transpose();
    size_t basisSize = 0;
    for (auto si = 0; si < 8; si++) {
        basisSize += hamiltonian.O1().dimension(si);
    }
    auto oddPacked = hamiltonian.m_hermiticity[0] == -1 && hamiltonian.m_hermiticity[1] == -1;
    if (oddPacked) return Eigen::MatrixXd::Zero(basisSize, basisSize);
    Eigen::MatrixXd result(basisSize, basisSize);
    auto hamO2 = hamiltonian.O2(spini > 0, spinj > 0);
    size_t i = 0;
    for (uint si = 0; si < 8; si++) {
        for (size_t oi = 0; oi < hamiltonian.dimension(si, 0, spini > 0); oi++, i++) {
            size_t j = 0;
            for (uint sj = 0; sj < 8; sj++) {
                for (size_t oj = 0; oj < hamiltonian.dimension(sj, 0, spinj > 0); oj++, j++) {
                    result(j, i) = hamO2.smat(0, si, oi, oi)->block(sj)[(oj + 2) * (oj + 1) / 2 - 1];
                }
            }
        }
    }
    return result;
}

Eigen::MatrixXd gci::intK(const SymmetryMatrix::Operator &hamiltonian, int spin) {
    size_t basisSize = 0;
    for (auto si = 0; si < 8; si++) {
        basisSize += hamiltonian.O1().dimension(si);
    }
    auto oddPacked = hamiltonian.m_hermiticity[0] == -1 && hamiltonian.m_hermiticity[1] == -1;
    Eigen::MatrixXd result(basisSize, basisSize);
    auto hamO2 = hamiltonian.O2(spin > 0, spin > 0);
    size_t i = 0;
    for (uint si = 0; si < 8; si++) {
        for (size_t oi = 0; oi < hamiltonian.dimension(si, 0, spin > 0); oi++, i++) {
            size_t j = 0;
            for (uint sj = 0; sj < 8; sj++) {
                for (size_t oj = 0; oj < hamiltonian.dimension(sj, 0, spin > 0); oj++, j++) {
                    if (si < sj) {
                        result(j, i) = hamO2.smat(si ^ sj, si, oi, oj)->blockMap(si)(oi, oj);
                    } else {
                        if (si > sj) {
                            result(j, i) = hamO2.smat(si ^ sj, sj, oj, oi)->blockMap(sj)(oj, oi);
                        } else {
                            if (i > j) {
                                if (oddPacked) {
                                    if (oi == oj) {
                                        result(j, i) = 0;
                                    } else {
                                        result(j, i) = hamO2.smat(0, si, oi, oj)->block(si)[(oi * (oi - 1)) / 2 + oj];
                                    }
                                } else {
                                    result(j, i) = hamO2.smat(0, si, oi, oj)->block(si)[(oi * (oi + 1)) / 2 + oj];
                                }

                            } else {
                                if (oddPacked) {
                                    if (oi == oj) {
                                        result(j, i) = 0;
                                    } else {
                                        result(j, i) =
                                                hamO2.smat(0, sj, oj, oi)->block(sj)[(oj * (oj - 1)) / 2 + oi];
                                    }
                                } else {

                                    result(j, i) =
                                            hamO2.smat(0, sj, oj, oi)->block(sj)[(oj * (oj + 1)) / 2 + oi];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

SymmetryMatrix::Operator gci::constructOperator(const FCIdump &dump) {
  std::vector<char> portableByteStream;
  int lPortableByteStream;
  int rank = 0;
#ifdef HAVE_MPI_H
  MPI_Comm_rank(MPI_COMM_COMPUTE, &rank);
#endif
  if (rank == 0) {
    int verbosity = 0;
    std::vector<int> orbital_symmetries = dump.parameter("ORBSYM");
    SymmetryMatrix::dim_t dim(8);
    for (const auto &s : orbital_symmetries)
      dim.at(s - 1)++;
    SymmetryMatrix::Operator result(dim, 2, dump.parameter("IUHF")[0] > 0, 0, true, false, "Hamiltonian");
    result.zero();

    dump.rewind();
    double value;
    FCIdump::integralType type;
    auto &integrals_a = result.O1(true);
    integrals_a.assign(0);
    auto &integrals_b = result.O1(false);
    integrals_b.assign(0);
    auto &integrals_aa = result.O2(true, true);
    auto &integrals_ab = result.O2(true, false);
    auto &integrals_bb = result.O2(false, false);
    if (verbosity > 0) {
      xout << "integral addresses " << &integrals_a << " " << &integrals_b << std::endl;
      xout << "integral addresses " << &integrals_a.block(0)[0] << " " << &integrals_b.block(0)[0] << std::endl;
      xout << "integral addresses " << &integrals_aa << " " << &integrals_ab << " " << &integrals_bb << std::endl;
      xout << "integral sizes " << integrals_aa.size() << " " << integrals_ab.size() << " " << integrals_bb.size()
           << std::endl;
    }
    unsigned int si,sj,sk,sl;
    size_t oi,oj,ok,ol;
    while ((type = dump.nextIntegral(si, oi, sj, oj, sk, ok, sl, ol, value)) != FCIdump::endOfFile) {
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;
      if (si < sj || (si == sj && oi < oj)) {
        std::swap(oi, oj);
        std::swap(si, sj);
      }
      if (sk < sl || (sk == sl && ok < ol)) {
        std::swap(ok, ol);
        std::swap(sk, sl);
      }
      unsigned int sij = si ^sj;
//      xout << "\nvalue: "<<value<<std::endl;
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;

      if (type == FCIdump::I2aa) {
        (sij ? integrals_aa.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) : integrals_aa.smat(sij, si, oi, oj)->block(sk)[
            ok * (ok + 1) / 2 + ol]) = value;
        (sij ? integrals_aa.smat(sij, sk, ok, ol)->blockMap(si)(oi, oj) : integrals_aa.smat(sij, sk, ok, ol)->block(si)[
            oi * (oi + 1) / 2 + oj]) = value;
      } else if (type == FCIdump::I2ab) {
        (sij ? integrals_ab.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) : integrals_ab.smat(sij, si, oi, oj)->block(sk)[
            ok * (ok + 1) / 2 + ol]) = value;
      } else if (type == FCIdump::I2bb) {
        (sij ? integrals_bb.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) : integrals_bb.smat(sij, si, oi, oj)->block(sk)[
            ok * (ok + 1) / 2 + ol]) = value;
        (sij ? integrals_bb.smat(sij, sk, ok, ol)->blockMap(si)(oi, oj) : integrals_bb.smat(sij, sk, ok, ol)->block(si)[
            oi * (oi + 1) / 2 + oj]) = value;
      } else if (type == FCIdump::I1a) {
        integrals_a.block(si).at(oi * (oi + 1) / 2 + oj) = value;
      } else if (type == FCIdump::I1b) {
        integrals_b.block(si).at(oi * (oi + 1) / 2 + oj) = value;
      } else if (type == FCIdump::I0)
        result.m_O0 = value;
    }
    if (verbosity > 0) xout << result << std::endl;
    if (verbosity > 1) xout << "int1:\n" << int1(result,1) << std::endl;
    if (verbosity > 1) xout << "intJ:\n" << intJ(result,1, 1) << std::endl;
    if (verbosity > 1) xout << "intK:\n" << intK(result,1) << std::endl;
    portableByteStream = result.bytestream().data();
    lPortableByteStream = portableByteStream.size();
  }
#ifdef HAVE_MPI_H
  MPI_Bcast(&lPortableByteStream, 1, MPI_INT, 0, MPI_COMM_COMPUTE);
#endif
  char *buf = (rank == 0) ? portableByteStream.data() : (char *) malloc(lPortableByteStream);
#ifdef HAVE_MPI_H
  MPI_Bcast(buf, lPortableByteStream, MPI_CHAR, 0, MPI_COMM_COMPUTE);
#endif
  class memory::bytestream bs(buf);
  auto result = SymmetryMatrix::Operator::construct(bs);
  if (rank != 0) free(buf);
  return result;
}

void gci::FCIDump(const SymmetryMatrix::Operator& op, const std::string filename, std::vector<int> orbital_symmetries) {
  FCIdump dump;
  int verbosity = 0;
  if (orbital_symmetries.empty())
    for (auto sym = 0; sym < 8; sym++)
      for (size_t i = 0; i < op.dimension(sym); i++)
        orbital_symmetries.push_back(sym+1);
  size_t n = orbital_symmetries.size();
  dump.addParameter("IUHF", op.m_uhf ? 1 : 0);
  dump.addParameter("ORBSYM", orbital_symmetries);
  dump.addParameter("NORB", int(orbital_symmetries.size()));

  dump.write(filename);
  dump.rewind();
  size_t i, j, k, l;
  const auto &integrals_a = op.O1(true);
  const auto &integrals_b = op.O1(false);
  const auto &integrals_aa = op.O2(true, true);
  const auto &integrals_ab = op.O2(true, false);
  const auto &integrals_bb = op.O2(false, false);
  if (verbosity > 0) {
    xout << "integral addresses " << &integrals_a << " " << &integrals_b << std::endl;
    xout << "integral addresses " << &integrals_a.block(0)[0] << " " << &integrals_b.block(0)[0] << std::endl;
    xout << "integral addresses " << &integrals_aa << " " << &integrals_ab << " " << &integrals_bb << std::endl;
//      xout << "integral sizes "<< integrals_aa.size()<<" "<<integrals_ab.size()<<" "<<integrals_bb.size()<<std::endl;
    xout << "n=" << n << std::endl;
    xout << "m_rank=" << op.m_rank << std::endl;
  }
  if (op.m_uhf) throw std::logic_error("UHF not supported");
  if (op.m_rank > 1) { // xout << "2 electron"<<std::endl;
    for (i = 1; i <= n; i++)
      for (j = 1; j <= i; j++) { // xout << "j="<<j<<std::endl;
        for (k = 1; k <= i; k++) { // xout << "k="<<k<<std::endl;
          for (l = 1; l <= (i == k ? j : k); l++) {
            auto oi = dump.orbital_offset(i);
            auto oj = dump.orbital_offset(j);
            auto ok = dump.orbital_offset(k);
            auto ol = dump.orbital_offset(l);
            auto si = dump.orbital_symmetry(i);
            auto sj = dump.orbital_symmetry(j);
            auto sk = dump.orbital_symmetry(k);
            auto sl = dump.orbital_symmetry(l);
//      xout << "ijkl "<<i<<j<<k<<l<<std::endl;
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;
            if (si < sj || (si == sj && oi < oj)) {
              std::swap(oi, oj);
              std::swap(si, sj);
            }
            if (sk < sl || (sk == sl && ok < ol)) {
              std::swap(ok, ol);
              std::swap(sk, sl);
            }
            unsigned int sij = si ^sj;
            unsigned int skl = sk ^sl;
            if ((sij ^ skl) != op.m_symmetry) continue;
//      xout << "\nvalue: "<<value<<std::endl;
//      xout << "s: ijkl "<<si<<sj<<sk<<sl<<std::endl;
//      xout << "o: ijkl "<<oi<<oj<<ok<<ol<<std::endl;
//      xout << (sij ? integrals_ab.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol) : integrals_ab.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol]) <<" "<< i << " "<<j<<" "<<k<<" "<<l <<std::endl;;
//      xout << (sij ? &integrals_ab.smat(sij,si,oi,oj)->blockMap(sk)(ok,ol) : &integrals_ab.smat(sij,si,oi,oj)->block(sk)[ok*(ok+1)/2+ol]) <<" "<< i << " "<<j<<" "<<k<<" "<<l <<std::endl;;
            auto value = (sij ? integrals_ab.smat(sij, si, oi, oj)->blockMap(sk)(ok, ol) : integrals_ab.smat(sij,
                                                                                                             si,
                                                                                                             oi,
                                                                                                             oj)->block(
                sk)[ok * (ok + 1) / 2 + ol]);
            if (std::abs(value) > 1e-15)
              dump.writeIntegral(i, j, k, l, value);
          }
        }
      }
  }
//  xout << "n="<<n<<std::endl;
  if (op.m_rank > 0)
    for (size_t ii = 1; ii <= n; ii++)
      for (size_t jj = 1; jj <= ii; jj++) {
//        xout << "ii="<<ii<<", jj="<<jj<<std::endl;
        auto oi = dump.orbital_offset(ii);
        auto oj = dump.orbital_offset(jj);
        auto si = dump.orbital_symmetry(ii);
        auto sj = dump.orbital_symmetry(jj);
        if ((si ^ sj) == op.m_symmetry) {
          auto value = integrals_a.block(si).at(oi * (oi + 1) / 2 + oj);
          if (std::abs(value) > 1e-15)
            dump.writeIntegral(ii, jj, 0, 0, value);
        }
      }
  dump.writeIntegral(0, 0, 0, 0, op.m_O0);

}

SymmetryMatrix::Operator gci::fockOperator(const SymmetryMatrix::Operator& hamiltonian, const Determinant &reference, const std::string description) {
  SymmetryMatrix::dim_t dims;
  for (auto s=0; s<8; s++)
    dims.push_back(hamiltonian.dimension(s,0,0));
  SymmetryMatrix::Operator f(dims,
             1,
             hamiltonian.m_uhf,
             hamiltonian.m_symmetry,
             true,
             false,
             description);
  // xout << "gci::Operator::fockOperator Reference alpha: "<<reference.stringAlpha<<std::endl;
  // xout << "gci::Operator::fockOperator Reference beta: "<<reference.stringBeta<<std::endl;
  //  bool closed = reference.stringAlpha==reference.stringBeta;
  //  xout << "gci::Operator::fockOperator Reference alpha=beta: "<<closed<<std::endl;
  auto refAlphaOrbitals = reference.stringAlpha.orbitals();
  auto refBetaOrbitals = reference.stringBeta.orbitals();
  f.m_O0 = hamiltonian.m_O0;
  f.O1(true) = hamiltonian.O1(true);
  if (hamiltonian.m_uhf) f.O1(false) = hamiltonian.O1(false);
  size_t basisSize = 0;
  for (auto si=0; si<8; si++)
    basisSize += hamiltonian.O1().dimension(si);
  // xout <<"reference.stringAlpha.orbitals ";for (size_t i=0; i < reference.stringAlpha.orbitals().size(); i++) xout <<reference.stringAlpha.orbitals()[i]<<" ";xout <<std::endl;
  for (const auto &o : refAlphaOrbitals) {
//       xout << "gci::Operator::fockOperator Reference alpha: "<<reference.stringAlpha<<std::endl;
//       xout<< "f alpha, alpha occ: " <<*o << std::endl;
    unsigned int os = reference.orbitalSpace->orbital_symmetries[o - 1];
    unsigned int oo = reference.orbitalSpace->orbitalIndex(o);
//    unsigned int oo = offset(o);
    for (unsigned int i = 1; i <= basisSize; i++)
      for (unsigned int j = 1; j <= i; j++) {
        if (reference.orbitalSpace->orbital_symmetries[i - 1] != reference.orbitalSpace->orbital_symmetries[j - 1]) continue;
        f.element(reference.orbitalSpace->orbitalIndex(i),
                  reference.orbitalSpace->orbital_symmetries[i - 1],
                  reference.orbitalSpace->orbitalIndex(j),
                  reference.orbitalSpace->orbital_symmetries[j - 1],
                  true) +=
            hamiltonian.element(reference.orbitalSpace->orbitalIndex(i),
                                reference.orbitalSpace->orbital_symmetries[i - 1],
                                reference.orbitalSpace->orbitalIndex(j),
                                reference.orbitalSpace->orbital_symmetries[j - 1],
                                oo,
                                os,
                                oo,
                                os,
                                true,
                                true) -
                hamiltonian.element(reference.orbitalSpace->orbitalIndex(i),
                                    reference.orbitalSpace->orbital_symmetries[i - 1],
                                    oo,
                                    os,
                                    oo,
                                    os,
                                    reference.orbitalSpace->orbitalIndex(j),
                                    reference.orbitalSpace->orbital_symmetries[j - 1],
                                    true,
                                    true);
      }
  }
  for (const auto &o : refBetaOrbitals) {
    // xout<< "f alpha, beta occ: " <<*o << std::endl;
    unsigned int os = reference.orbitalSpace->orbital_symmetries[o - 1];
    unsigned int oo = reference.orbitalSpace->orbitalIndex(o);
    for (unsigned int i = 1; i <= basisSize; i++)
      for (unsigned int j = 1; j <= i; j++) {
        if (reference.orbitalSpace->orbital_symmetries[i - 1] != reference.orbitalSpace->orbital_symmetries[j - 1]) continue;
        f.element(reference.orbitalSpace->orbitalIndex(i),
                  reference.orbitalSpace->orbital_symmetries[i - 1],
                  reference.orbitalSpace->orbitalIndex(j),
                  reference.orbitalSpace->orbital_symmetries[j - 1],
                  true) +=
            hamiltonian.element(reference.orbitalSpace->orbitalIndex(i),
                                reference.orbitalSpace->orbital_symmetries[i - 1],
                                reference.orbitalSpace->orbitalIndex(j),
                                reference.orbitalSpace->orbital_symmetries[j - 1],
                                oo,
                                os,
                                oo,
                                os,
                                true,
                                false);
      }
  }
  if (f.m_uhf) {
    for (const auto &o : refBetaOrbitals) {
      // xout<< "f beta, beta occ: " <<*o << std::endl;
      unsigned int os = reference.orbitalSpace->orbital_symmetries[o - 1];
      unsigned int oo = reference.orbitalSpace->orbitalIndex(o);
      for (unsigned int i = 1; i <= basisSize; i++)
        for (unsigned int j = 1; j <= i; j++) {
          if (reference.orbitalSpace->orbital_symmetries[i - 1] != reference.orbitalSpace->orbital_symmetries[j - 1]) continue;
          f.element(reference.orbitalSpace->orbitalIndex(i),
                    reference.orbitalSpace->orbital_symmetries[i - 1],
                    reference.orbitalSpace->orbitalIndex(j),
                    reference.orbitalSpace->orbital_symmetries[j - 1],
                    false) +=
              hamiltonian.element(reference.orbitalSpace->orbitalIndex(i),
                                  reference.orbitalSpace->orbital_symmetries[i - 1],
                                  reference.orbitalSpace->orbitalIndex(j),
                                  reference.orbitalSpace->orbital_symmetries[j - 1],
                                  oo,
                                  os,
                                  oo,
                                  os,
                                  false,
                                  false) -
                  hamiltonian.element(reference.orbitalSpace->orbitalIndex(i),
                                      reference.orbitalSpace->orbital_symmetries[i - 1],
                                      oo,
                                      os,
                                      oo,
                                      os,
                                      reference.orbitalSpace->orbitalIndex(j),
                                      reference.orbitalSpace->orbital_symmetries[j - 1],
                                      false,
                                      false);
        }
    }
    for (const auto &o : refAlphaOrbitals) {
      // xout<< "f beta, alpha occ: " <<*o << std::endl;
      unsigned int os = reference.orbitalSpace->orbital_symmetries[o - 1];
      unsigned int oo = reference.orbitalSpace->orbitalIndex(o);
      for (unsigned int i = 1; i <= basisSize; i++)
        for (unsigned int j = 1; j <= i; j++) {
          if (reference.orbitalSpace->orbital_symmetries[i - 1] != reference.orbitalSpace->orbital_symmetries[j - 1]) continue;
          f.element(reference.orbitalSpace->orbitalIndex(i),
                    reference.orbitalSpace->orbital_symmetries[i - 1],
                    reference.orbitalSpace->orbitalIndex(j),
                    reference.orbitalSpace->orbital_symmetries[j - 1],
                    false) +=
              hamiltonian.element(oo,
                                  os,
                                  oo,
                                  os,
                                  reference.orbitalSpace->orbitalIndex(i),
                                  reference.orbitalSpace->orbital_symmetries[i - 1],
                                  reference.orbitalSpace->orbitalIndex(j),
                                  reference.orbitalSpace->orbital_symmetries[j - 1],
                                  true,
                                  false);
        }
    }
  }
  return f;
}

SymmetryMatrix::Operator gci::sameSpinOperator(const SymmetryMatrix::Operator& hamiltonian, const Determinant &reference, const std::string description) {
  SymmetryMatrix::dim_t dims;
  for (auto s=0; s<8; s++)
    dims.push_back(hamiltonian.dimension(s,0,0));
  SymmetryMatrix::Operator result(dims,
                  hamiltonian.m_rank,
                  true,
                  hamiltonian.m_symmetry,
                  true,
                  false,
                  description);
  {
    Determinant ra = reference;
    ra.stringBeta.nullify();
    auto f = fockOperator(hamiltonian,ra);
    for (size_t i = 0; i < hamiltonian.O1(true).size(); i++)
      (*result.O1(true).data())[i] = (*hamiltonian.O1(true).data())[i] - (*f.O1(true).data())[i];
  }
  {
    Determinant ra = reference;
    ra.stringAlpha.nullify();
    auto f = fockOperator(hamiltonian,ra);
    for (size_t i = 0; i < hamiltonian.O1(false).size(); i++)
      (*result.O1(false).data())[i] = (*hamiltonian.O1(false).data())[i] - (*f.O1(false).data())[i];
  }

  *result.O2(true, true).data() = *hamiltonian.O2(true, true).data();
  *result.O2(false, false).data() = *hamiltonian.O2(false, false).data();
  for (auto &s : *result.O2(true, false).data()) s = 0;
  return result;
}

void gci::gsum(SymmetryMatrix::Operator& op) {
  std::vector<double> O0;
  O0.push_back(op.m_O0);
  ::gci::gsum(O0);
  op.m_O0 = O0[0];
  if (op.m_rank > 0) {
    ::gci::gsum(*op.O1(true).data());
    if (op.m_uhf)
      ::gci::gsum(*op.O1(false).data());
  }
  if (op.m_rank > 1) {
    ::gci::gsum(*op.O2(true, true).data());
    if (op.m_uhf) {
      ::gci::gsum(*op.O2(true, false).data());
      ::gci::gsum(*op.O2(false, false).data());
    }
  }
}
#ifdef MOLPRO
#ifdef __cplusplus
extern "C" {
#endif
void gcirun(double* energies, int nenergies, char* fcidump) {
  Run run(fcidump);
  try {
    std::vector<double>e = run.run();
    for (int i=0; i < (nenergies > (int)e.size() ? (int)e.size() : nenergies); i++)
      energies[i]=e[i];
  }
  catch(char const* c) { xout << "caught error: " <<c<<std::endl;}
  return;
}
#ifdef __cplusplus
}
#endif
#endif
