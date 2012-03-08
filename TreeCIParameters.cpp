#include "TreeCIParameters.h"
using namespace std;
#include <fstream>
#include <istream>
#include <sstream>

TreeCIParameters::TreeCIParameters(string filename)
{
    if (filename!="") load(filename);
}

TreeCIParameters::~TreeCIParameters()
{
}


void TreeCIParameters::load(string filename) { // dirty sucking in from FCIDUMP namelist
//    cout << "loadParameters " <<file << endl;
    ifstream s;
    s.open(filename.c_str());
    if ( (s.rdstate() & ifstream::failbit ) != 0 ) {
        cerr << "Error opening " << filename <<endl;
        throw "file missing";
    }
    string ss, sslast="";
    while (s >> ss && ss != "&END") {
        string first = ss.substr(0,1);
//        cout << "ss="<<ss<<" sslast="<<sslast<<" first="<<first <<endl;
        const char * d1=first.c_str();
        if (isdigit(*d1)) {
//            cout <<"first is digit"<<endl;
            int i;
            std::istringstream b(ss);
            b >>i;
//            cout <<"i="<<i<<endl;
            char ccc[100]; sprintf(ccc,"%d",i);
            string sslast2(sslast);sslast2.append(ccc);
            sslast=ss.substr(ss.find(',')+1);
            ss=sslast2;
        } else {
            sslast=ss;
        }
//        cout <<"process ss="<<ss<<endl;
        int i, ipos;
        ipos=ss.find('=');
        istringstream bb(ss.substr(ipos+1));
        bb >> i;
        string key(ss.substr(0,ipos));
//        cout <<"found integer "<<i <<" for "<<key<<endl;
        if (key == "NORB") basisSize=i;
        if (key == "NELEC") nelec=i;
        if (key == "MS2") ms2=i;
    }
    cout <<"nelec="<<nelec<<endl;
    cout <<"ms2="<<ms2<<endl;
    cout <<"basisSize="<<basisSize<<endl;
    s.close();
}
