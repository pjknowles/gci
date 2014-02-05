#include "gciExcitationSet.h"
#include <iostream>
#include <sstream>
Excitation::Excitation(long StringIndex, int Phase, long OrbitalAddress)
{
    stringIndex = StringIndex;
    phase = Phase;
    orbitalAddress = OrbitalAddress;
}

ExcitationSet::ExcitationSet(String &from, StringSet &to, int annihilations, int creations)
{
    int symexc=-1;
    From = from;
    xout <<"ExcitationSet taking From="<<From.printable()<<" annihilations="<<annihilations<<" creations="<<creations<<std::endl;
    To = &to;
    if (to.symmetry>=0) symexc = from.symmetry ^ to.symmetry; // use symmetry if we can
    if (annihilations==1 && creations==0) {
        for (int i=0; i<(int)from.orbitals_.size(); i++) {
            if (from.hamiltonian->orbital_symmetries[i]==(unsigned int)symexc || symexc==-1) {
                String tt = from;
                int phase = tt.destroy(i+1);
                if (phase) {
                    xout << "found annihilation from "<<from.printable()<<" to "<<tt.printable()<<" phase="<<phase<<" orbital="<<i+1<<std::endl;
                    long ti=to.offset(tt);
                    push_back(Excitation(ti,phase,i));
                }
            }
        }
    }
    else if (annihilations==0 && creations==1) {
        for (int i=0; i<(int)from.orbitals_.size(); i++) {
            if (from.hamiltonian->orbital_symmetries[i]==(unsigned int)symexc || symexc==-1) {
                String tt = from;
                int phase = tt.create(i+1);
                if (phase) {
                    xout << "found creation from "<<from.printable()<<" to "<<tt.printable()<<" phase="<<phase<<" orbital="<<i+1<<std::endl;
                    long ti=to.offset(tt);
                    push_back(Excitation(ti,phase,i));
                }
            }
        }
    }
    else if (annihilations==1 && creations==1) {
//        String proto(from); proto.destroy(1);
//        StringSet ann1(proto,true,-1); // set of all annihilations
        for (int j=0; j<(int)from.orbitals_.size(); j++) {
            String a = from;
            int phasea = a.destroy(j+1);
//            xout <<"j, phase, a "<<j+1<<phasea<<a.printable()<<std::endl;
            if (phasea) {
                for (int i=0; i<(int)from.orbitals_.size(); i++) {
                    if (from.hamiltonian->orbital_symmetries[i]==((unsigned int)symexc^from.hamiltonian->orbital_symmetries[j]) || symexc==-1) {
                        String tt = a;
                        int phase = phasea*tt.create(i+1);
//                        xout <<"i, phase, tt "<<i+1<<phase<<std::endl;
                            if (phase) {
                                long ti=to.offset(tt);
                                if (ti < 0 || (size_type) ti>= to.size()) {
                                    xout <<"i="<<i+1<<" phase="<<phase<<" ti="<<ti<<" tt="<<tt.printable()<<std::endl;
                                    throw "index error in ExcitationSet";
                                }
                                xout << "found Excitation i="<<i+1<<", j="<<j+1<<",pairIndex="<<from.hamiltonian->pairIndex(i+1,j+1)<<std::endl;
                                push_back(Excitation(ti,phase,from.hamiltonian->pairIndex(i+1,j+1)));
                            }
                    }
                }
            }
        }
    }
}

std::string ExcitationSet::printable() {
    if ((*To).size()==0) return "";
    std::stringstream ss;
    ss<<"Excitations for String "<<From.printable() << " into symmetry " << (*To)[0].symmetry+1 <<":";
    for (ExcitationSet::iterator e=this->begin(); e!=this->end(); e++) {
        ss<<std::endl<< " String index="<< e->stringIndex
         <<"("<<(*To)[e->stringIndex].printable()<<")"
        <<" phase="<<e->phase<<" orbitalAddress="<<e->orbitalAddress;
    }
    return ss.str();
}
