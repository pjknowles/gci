#include "gciStringSet.h"
#include <iostream>

StringSet::StringSet() : vector<String>()
{
}

StringSet::StringSet(String prototype, bool all) : vector<String>()
{
    // copy prototype
    proto = prototype;
    { // set up partial weight array for addressing binomial distributions
        int nitem = proto.nelec;
        int nbox = proto.hamiltonian->basisSize;
        vector< vector<int> > inter(nitem, vector<int>(nbox));
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
        //    xout << "PartialWeightArray:"<<std::endl; for (int k=0;k<nitem;k++) { for (int l=0;l<nbox;l++) xout << " "<<PartialWeightArray[k][l]; xout <<std::endl; }
    }
    if (all) complete();

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

void StringSet::complete()
{
    xout <<"StringSet::complete prototype"<<proto.printable(1)<<std::endl;
    String string(&proto);
    string.first(proto.nelec);
    this->erase(this->begin(),this->end());
    do {
//        xout << "in StringSet::complete about to push_back " << string.printable(1) <<std::endl;
        this->push_back(string);
    } while (string.next());
    xout << "in StringSet::complete final list: " <<std::endl ;
    for (iterator s=this->begin(); s!=this->end(); s++) xout << s->printable()<<std::endl;

}

