#include "gciStringSet.h"
#include <iostream>
#include <sstream>

StringSet::StringSet() : std::vector<String>()
{
//    xout <<"StringSet default constructor"<<std::endl;
}

StringSet::StringSet(String prototype, bool all, int sym) : std::vector<String>()
{
//    xout <<"StringSet prototype constructor "<<all<<std::endl;
    // copy prototype
    proto = prototype;
    { // set up partial weight array for addressing binomial distributions
        int nitem = proto.nelec;
        int nbox = proto.orbitalSpace->total();
        std::vector< std::vector<int> > inter(nitem, std::vector<int>(nbox));
        //    xout << "in StringSet::PartialWeightArray "<<nitem << " " << nbox <<std::endl;
        for (int k=0;k<nitem;k++) {
            for (int l=0;l<nbox;l++) {
                inter[k][l]=0;
            }
            for (int l=k;l<nbox-nitem+k;l++)
                inter[k][l+1] = binomial_coefficient(nbox-l-1,nitem-k-1)+inter[k][l];
        }
        for (int k=0;k<nitem-1;k++) {
            for (int l=k;l<nbox-nitem+k+1;l++) {
                inter[k][l] = inter[k][l] - inter[k+1][l+1];
            }
        }
        for (int l=nitem-1;l<nbox;l++) {
            inter[nitem-1][l] = l-nitem+1;
        }

        PartialWeightArray=inter;
            xout << "PartialWeightArray:"<<std::endl; for (int k=0;k<nitem;k++) { for (int l=0;l<nbox;l++) xout << " "<<PartialWeightArray[k][l]; xout <<std::endl; }
    }
    if (all) complete(sym);
    symmetry = sym;
//    calculateAddressMap();
    // compute address map
}


long StringSet::binomial_coefficient(unsigned long n, unsigned long k) {
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

void StringSet::complete(int sym)
{
//    xout <<"StringSet::complete prototype"<<proto.printable(1)<<std::endl;
    String string(&proto);
    this->erase(this->begin(),this->end());
    if (string.first(proto.nelec,sym)) {
        //    xout <<"StringSet::complete symmetry="<<sym<<" first String: "<<string.printable()<<std::endl;
        do {
            //        xout << "in StringSet::complete about to push_back " << string.printable(1) <<std::endl;
            this->push_back(string);
        } while (string.next(sym));
    }
    //    xout << "in StringSet::complete final list: " <<std::endl ;
//    for (iterator s=this->begin(); s!=this->end(); s++) xout << s->printable()<<std::endl;

}

void StringSet::push_back(String& s)
{
    s.key=0;
    for (int k=0; k<(int)s.orbitals_.size(); k++)
        s.key+= PartialWeightArray[k][s.orbitals_[k]-1];
    addressMap[s.key]=size();
    xout <<"StringSet::push_back " <<s <<" size()=" <<size()<<std::endl;
    std::vector<String>::push_back(s);
}

std::string StringSet::str(int verbosity) const
{
    std::ostringstream s;
    if (verbosity >= -1) {
        s << "StringSet " << size();
    }
    return s.str();
}

std::vector<ExcitationSet> StringSet::allExcitations(StringSet &to, int annihilations, int creations)
{
    std::vector<ExcitationSet> set;
    for (iterator f=begin(); f!=end(); f++) {
        String ff = *f;
        set.push_back(ExcitationSet(ff,to,annihilations,creations));
    }
    return set;
}

std::vector<double> StringSet::occupationNumbers()
{
    std::vector<double> result;
    if (this->size()) {
        String firstString=this->at(0);
        result.resize(this->size()* firstString.orbitalSpace->total(), (double) 0);
        int stringoffset=0;
        for (StringSet::iterator s=this->begin(); s!=this->end(); s++)
        {
            std::vector<unsigned int> orbitals = s->orbitals();
            for (std::vector<unsigned int>::iterator i=orbitals.begin(); i !=orbitals.end(); i++) {
//                xout << "StringSet::occupationNumbers stringoffset="<<stringoffset<<" *i="<<*i<<std::endl;
                result[stringoffset+(*i-1)*this->size()]=(double)1;
            }
            stringoffset++;
        }
    }
    return result;
}
