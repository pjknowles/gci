#include "gciExcitationSet.h"
#include <iostream>
#include <sstream>

//using ExcitationSet = gci::ExcitationSet;
namespace gci {

ExcitationSet::ExcitationSet(const String &from, const StringSet &to, int annihilations, int creations,
                             SymmetryMatrix::parity_t parity)
        : From(from), To(to) {
    constexpr String::orbital_type one = 1;
    int symexc = -1;
    //    xout <<"ExcitationSet taking From="<<From.str()<<" annihilations="<<annihilations<<" creations="<<creations<<std::endl;
    if (to.symmetry >= 0) symexc = from.symmetry ^ static_cast<unsigned int>(to.symmetry); // use symmetry if we can
//  if (false && annihilations == 0 && creations == 1) {
//    int ii = 0;
//    for (int i = 0; i < (int) from.orbitalSpace->orbital_symmetries.size(); i++) {
//      if (from.orbitalSpace->orbital_symmetries[i] == (unsigned int) symexc || symexc == -1) {
//        String tt = from;
//        int phase = tt.create(i + 1);
//        if (phase) {
//          tt.gci::String::keygen(to.PartialWeightArray);
//          size_t ti = to.addressMap.find(tt.key())->second;
//          buffer.emplace_back(ti, phase, ii);
//        }
//        ii++;
//      }
//    }
//  } else if (false && annihilations == 1 && creations == 0) {
//    int ii = 0;
//    for (int i = 0; i < (int) from.orbitalSpace->orbital_symmetries.size(); i++) {
//      if (from.orbitalSpace->orbital_symmetries[i] == (unsigned int) symexc || symexc == -1) {
//        String tt = from;
//        int phase = tt.destroy(i + 1);
//        if (phase) {
//          tt.gci::String::keygen(to.PartialWeightArray);
//          size_t ti = to.addressMap.find(tt.key())->second;
//          buffer.emplace_back(ti, phase, ii);
//        }
//        ii++;
//      }
//    }
//  } else
    if (annihilations + creations == 1) {
        int ii = 0;
        buffer.reserve(
                from.orbitalSpace->orbital_symmetries.size());
        for (String::orbital_type i = 0; i < (int) from.orbitalSpace->orbital_symmetries.size(); i++) {
            if (from.orbitalSpace->orbital_symmetries[i] == (unsigned int) symexc || symexc == -1) {
                String tt = from;
                int phase = (annihilations > 0) ? tt.destroy(i + one) : tt.create(i + one);
                if (phase) {
                    tt.String::keygen(to.PartialWeightArray);
                    size_t ti = to.addressMap.find(tt.key())->second;
                    buffer.emplace_back(ti, phase, ii);
                }
                ii++;
            }
        }
    } else if (annihilations + creations == 2) {
//    int parity = annihilations==creations ? 1 : -1 ; // could be changed
//            xout << "two-orbital excitation, parity="<<parity<<std::endl;
//            xout <<"from="<<from<<std::endl;
        for (String::orbital_type j = 0; j < (int) from.orbitalSpace->orbital_symmetries.size(); j++) {
            String a = from;
            int phasea = (annihilations > 0) ? a.destroy(j + one) : a.create(j + one);
            if (phasea) {
                for (String::orbital_type i = 0;
                     i < (int) (annihilations == creations ? from.orbitalSpace->orbital_symmetries.size() : j);
                     i++) {
                    if (from.orbitalSpace->orbital_symmetries[i]
                        == ((unsigned int) symexc ^ from.orbitalSpace->orbital_symmetries[j]) || symexc == -1) {
                        String tt = a;
                        int phase = (annihilations > 1) ? phasea * tt.destroy(i + one) : phasea * tt.create(i + one);
                        if (phase) {
                            tt.keygen(to.PartialWeightArray);
                            size_t ti = to.addressMap.find(tt.key())->second;
                            if (ti >= to.size()) {
                                xout << "i=" << i + one << " phase=" << phase << " ti=" << ti << " tt=" << tt.str()
                                     << std::endl;
                                throw std::range_error("index error in ExcitationSet");
                            }
                            if (from.orbitalSpace->orbital_symmetries[i] >
                                from.orbitalSpace->orbital_symmetries[j])
                                phase = -phase;
                            buffer.emplace_back(ti, phase, from.orbitalSpace->pairIndex(i + one, j + one, parity));
                        }
                    }
                }
            }
        }
    }
//  buffer.shrink_to_fit();
}

std::string ExcitationSet::str(int verbosity) const {
    if (To.empty() || verbosity < 0) return "";
    std::stringstream ss;
    ss << "Excitations for String " << From.str() << " into symmetry " << To[0].symmetry + 1 << ":";
    for (const auto &e : buffer) {
        ss << std::endl << " String index=" << e.stringIndex
           << "(" << To[e.stringIndex].str() << ")"
           << " phase=" << e.phase << " orbitalAddress=" << e.orbitalAddress;
    }
    return ss.str();
}

std::ostream &operator<<(std::ostream &os, ExcitationSet const &obj) {
    return os << obj.str();
}
}// namespace gci
