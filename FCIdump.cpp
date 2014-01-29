#include "FCIdump.h"
#include <iostream>

FCIdump::FCIdump(std::string filename)
{
    fileName = filename;
}

#include "gciState.h"
using namespace std;
#include <fstream>
#include <istream>
#include <sstream>



std::vector<int> FCIdump::parameter(std::string key, std::vector<int> def) { // dirty sucking in from FCIDUMP namelist
    std::vector<int> answer;
    ifstream s;
    s.open(fileName.c_str());
    if ( (s.rdstate() & ifstream::failbit ) != 0 ) {
        cerr << "Error opening " << fileName <<endl;
        throw "FCIDUMP::parameter file missing";
    }
    string ss;
    string keyin("");// the current keyword
    while (s >> ss && ss != "&END" && ss != "/") {
        while (ss.size()) {
            // gobble up leading spaces
            while (ss.size() && ss[0]==' ') ss.erase(0,1);
            // gobble the next thing up to equals, comma, end of line
            size_t iend = ss.find_first_of("=,"); if (iend==std::string::npos) iend=ss.size();
            string thisword = ss.substr(0,iend);
            ss=ss.substr(iend >= ss.size()-1 ? ss.size() : iend+1 );
//            xout << "thisword=" <<thisword <<"; new ss=" <<ss <<std::endl;
            if (isdigit(thisword[0])) {
                int i;
                std::istringstream b(thisword);
                b >>i;
//                xout <<"first is digit =" <<i <<endl;
                if (key==keyin) answer.push_back(i);
//                if (key==keyin) xout << "taking this for key="<<key<<endl;
            } else {
                keyin=thisword;
//                xout <<"keyin="<<keyin<<endl;
            }
        }
    }
    s.close();

    xout << "parameter "<<key<<"="; for (std::vector<int>::iterator s=answer.begin(); s < answer.end(); s++) xout <<*s ; xout <<endl;
    if (answer.size()==0) return def;
    return answer;
}
