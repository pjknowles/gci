#include "gciMixedWavefunction.h"

MixedWavefunction::MixedWavefunction(const State &state, int nMode, int nModal, int modeCoupling)
        : nMode_(nMode), nModal_(nModal), modeCoupling_(modeCoupling) {
    assert((modeCoupling == 1 || modeCoupling == 2) && "Current Hamiltonian assumes 1 or 2 mode coupling");
    // Zeroth order
    wfn_.emplace_back(state);
    // First order coupling
    for (int iMode = 0; iMode < nMode_; ++iMode) {
        for (int iModal = 0; iModal < nModal_; ++iModal) {
            wfn_.emplace_back(state);
        }
    }
    // Second order coupling
    if (modeCoupling_ > 1) {
        throw std::logic_error("modeCoupling > 1 is not yet implemented");
    }
}

bool MixedWavefunction::compatible(const MixedWavefunction &w2) const {
    bool b_size = (wfn_.size() == w2.wfn_.size());
    if (!b_size) return b_size;
    bool b_electronic_wfn = true;
    for (int i = 0; i < wfn_.size(); ++i) {
        b_electronic_wfn = b_electronic_wfn && wfn_[i].compatible(w2.wfn_[i]);
    }
    bool b_vib_basis = (nMode_ == w2.nMode_)&&(nModal_ == w2.nModal_);
    return b_size && b_electronic_wfn && b_vib_basis;
}

void MixedWavefunction::axpy(double a, const LinearAlgebra::vector<double> &other) {
    const auto &x = dynamic_cast <const MixedWavefunction &> (other);
    for (int i = 0; i < wfn_.size(); ++i) {
        wfn_[i].axpy(a, x.wfn_[i]);
    }
}


void MixedWavefunction::scal(double a) {
    for (auto &el: wfn_) el.scal(a);
}

double MixedWavefunction::operator*(const MixedWavefunction &w1, const MixedWavefunction &w2) {
    if (!w1.compatible(w2))
        throw std::domain_error("attempt to form scalar product between incompatible Wavefunction objects");
    for (int i = 0; i < wfn_.size(); ++i) {
        wfn_[i].axpy(a, x.wfn_[i]);
    }
}

double MixedWavefunction::dot(const LinearAlgebra::vector<double> &other) const {
    return (*this) * ((dynamic_cast<const MixedWavefunction &>(other)));
}

void MixedWavefunction::times(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b) {

}

void
MixedWavefunction::divide(const LinearAlgebra::vector<double> *a, const LinearAlgebra::vector<double> *b, double shift,
                          bool append, bool negative) {

}

void MixedWavefunction::zero() {

}

