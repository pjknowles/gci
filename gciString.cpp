#include "gciString.h"

#include <iostream>
#include <sstream>

String::String(State* State, int Spin)
{
    if (State == NULL) {
        hamiltonian=NULL;
    } else {
        hamiltonian=State->hamiltonian;
    }
    nullify();
    spin=Spin;
}

int String::create(unsigned int orbital) {
        xout  << "create orbital "<<orbital <<" " <<orbitals_.size()<<std::endl;
//        xout << "hamiltonian "<<(hamiltonian!=NULL)<<std::endl;
        if (hamiltonian==NULL)
            throw "String::create missing hamiltonian";
//        xout << "basisSize "<<hamiltonian->basisSize<<std::endl;
    if (hamiltonian==NULL || orbital==(unsigned int)0 || orbital > (unsigned int) hamiltonian->basisSize) throw "invalid orbital";
//    xout <<"make iterator "<<std::endl;
    std::vector<unsigned int>::iterator ilast=orbitals_.begin();
//    xout <<"iterator OK"<<std::endl;

    int phase=((orbitals_.size()/2)*2 == orbitals_.size()) ? 1 : -1;
    xout <<"phase="<<phase<<std::endl;
    xout <<"spin="<<spin<<std::endl;
    ms2+=spin;
    nelec++;
    if (orbitals_.size() == 0 || *ilast > orbital) {
        xout <<"first "<<orbital<<std::endl;
        orbitals_.insert(orbitals_.begin(),orbital);
        xout <<"orbitals_[0]="<<orbitals_[0]<<std::endl;
        return phase;
    }
    for (std::vector<unsigned int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
        if (*i==orbital) return 0; // exclusion principle
        if (*ilast < orbital && *i > orbital){
            orbitals_.insert(++ilast,orbital);
            return phase;
        }
        phase=-phase;
        ilast=i;
    }
    orbitals_.insert(orbitals_.end(),orbital);
    return phase;
}

int String::destroy(unsigned int orbital) {
    if (hamiltonian==NULL || orbital==(unsigned int)0 || orbital > (unsigned int) hamiltonian->basisSize ) throw "invalid orbital";
    if (orbitals_.size() <= 0) throw "too few electrons in String";
    int phase=1;
    ms2-=spin;
    nelec--;
    for (std::vector<unsigned int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
        if (*i==orbital)  {
            orbitals_.erase(i);
            return phase;
        }
        phase=-phase;
    }
    return (int)0; // exclusion principle
}

void String::nullify()
{
    orbitals_.clear();
    ms2=0;
    nelec=0;
    symmetry=1;
}

std::vector<unsigned int> String::orbitals() {
    return orbitals_;
}

std::string String::printable(int verbosity) {
    std::string result;
//    xout <<"String::printable orbitals_[0]" <<orbitals_[0]<<std::endl;
    if (verbosity >=0) {
        for (std::vector<unsigned int>::iterator i = orbitals_.begin(); i!=orbitals_.end(); ++i) {
            if (i!=orbitals_.begin()) result.append(",");
            std::stringstream ss;
            int ispin=(int)(*i)*(int)spin;
            ss << ispin;
            std::string rr;
            ss >> rr;
            result.append(rr);
        }
    }
    return result;
}

bool String::next() {
    unsigned int limit=hamiltonian->basisSize;
    std::vector<unsigned int>::iterator k;
    unsigned int floor;
    for (std::vector<unsigned int>::reverse_iterator i = orbitals_.rbegin(); i!=orbitals_.rend(); ++i) {
        floor=++(*i);
        k=i.base();
        if (*i <= limit) break;
        limit--;
    }
    if (limit <= hamiltonian->basisSize-orbitals_.size()) return false; // we ran out of boxes to put the objects into
    while (k!=orbitals_.rbegin().base())
        *(k++)=++floor;
    return true;
}

void String::first(int n) {
    xout << "String::first " << n <<" nelec="<<nelec<< std::endl;
    if (n <=0 ) n=nelec;
    if (n <=0 ) n=orbitals_.size();
    nullify();
    xout << n <<std::endl;
    for (unsigned int i=1;i<=(unsigned int)n;i++)
        create(i);
    nelec=n;
    xout <<"OK"<<nelec<<std::endl;
    xout << "String::first result " << printable(1) << std::endl;
}

String String::exhausted;

void String::buildStrings(State prototype, std::vector<String> *strings)
{
//    xout <<"buildStrings prototype"<<prototype.printable(1)<<std::endl;
    String string(&prototype);
    string.first(prototype.nelec);
    strings->erase(strings->begin(),strings->end());
    do {
//        xout << "in buildStrings about to push_back " << string.printable(1) <<std::endl;
        strings->push_back(string);
    } while (string.next());
//    xout << "in buildStrings final list: " <<std::endl ;
//    for (std::vector<String>::iterator s=strings->begin(); s!=strings->end(); s++) xout << s->printable()<<std::endl;

    vector<vector<int> > inter=PartialWeightArray(prototype.nelec,prototype.hamiltonian->basisSize);
}

vector< vector<int> > String::PartialWeightArray (int nitem,int nbox)
// sets up partial weight array for addressing binomial distributions
{
    vector< vector<int> > inter(nitem, vector<int>(nbox));
    xout << "in String::PartialWeightArray "<<nitem << " " << nbox <<std::endl;
    for (int k=0;k<nitem;k++) {
        for (int l=0;l<nbox;l++) {
            inter[k][l]=0;
        }
        for (int l=k;l<nbox-nitem+k;l++)
            inter[k][l+1] = binomial_coefficient(nbox-l-1,nitem-k-1)+inter[k][l];
    }
    xout << "Inter aka PartialWeightArray:"<<std::endl; for (int k=0;k<nitem;k++) { for (int l=0;l<nbox;l++) xout << " "<<inter[k][l]; xout <<std::endl; }
    for (int k=0;k<nitem-1;k++) {
        for (int l=k;l<nbox-nitem+k+1;l++) {
            inter[k][l] = inter[k][l] - inter[k+1][l+1];
        }
    }
    xout << "Inter aka PartialWeightArray:"<<std::endl; for (int k=0;k<nitem;k++) { for (int l=0;l<nbox;l++) xout << " "<<inter[k][l]; xout <<std::endl; }
    for (int l=nitem-1;l<nbox;l++) {
      inter[nitem-1][l] = l-nitem+1;
    }

    xout << "Inter aka PartialWeightArray:"<<std::endl; for (int k=0;k<nitem;k++) { for (int l=0;l<nbox;l++) xout << " "<<inter[k][l]; xout <<std::endl; }

    return inter;
}

long gci::binomial_coefficient(unsigned long n, unsigned long k) {
    unsigned long i;
    long b;
    if (0 == k || n == k) {
        return 1;
    }
    if (k > n) {
        return 0;
    }
    if (k > (n - k)) {
        k = n - k;
    }
    if (1 == k) {
        return n;
    }
    b = 1;
    for (i = 1; i <= k; ++i) {
        b *= (n - (k - i));
        if (b < 0) return -1; /* Overflow */
        b /= i;
    }
    return b;
}
