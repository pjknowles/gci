#include "gciSymmetry.h"
#include <iostream>
#include <sstream>

SymmetryOffset::SymmetryOffset(std::string title, int maxrank)
{
    this->maxrank = maxrank;
    Title=title;
    resize(8,0);
}

std::string SymmetryOffset::toString(int verbosity) const
{
    if (verbosity < 0) return std::string("");
    std::ostringstream s;
    s << "SymmetryOffset object " << Title << ", dimensions:";
    for (const_iterator x=begin(); x!=end(); x++) {
        s << " " << *x;
    }
    if (verbosity > 0 && ! offsets.empty()) {
        s << std::endl << "Vector offsets:";
        for (int i=0; i<9; i++) s << " " << offsets[i] ;
        s << std::endl << "General matrix offsets:";
        for (int i=0; i<9; i++) {
            s <<std::endl;
            for (int j=0; j<8; j++) s << " " << this->offset(j,i,0);
        }
        s << std::endl << "Symmetric matrix offsets:";
        for (int i=0; i<9; i++) {
            s <<std::endl;
            for (int j=0; j<8; j++) s << " " << this->offset(j,i,1);
        }
        if (verbosity > 1) {
            s <<std::endl;
            s << "Raw offsets:";
            for (int i=0; i<offsets.size(); i++) s << " " << offsets[i] ;
        }
    }
    return s.str();
}

void SymmetryOffset::calculateOffsets()
{
    size_t dimension=9;
    if (maxrank > 1) dimension += 8*9*2;
    offsets.resize(dimension,0);
    if (maxrank >= 1) {
        for (int s=0; s<8; s++)
            offsets[s+1] = offsets[s] + this->at(s);
    }
    xout << "calculateOffsets sizes="; for (int i=0; i<8; i++) xout << this->at(i) << " "; xout <<std::endl;
    if (maxrank >= 2) {
        for (int sym=0; sym<8; sym++) {
                size_t ntqg=0; size_t ntdg=0;
            for (int isym=0; isym<8; isym++) {
                int jsym = sym ^ isym;
                offsets[9 + sym*9 + isym] = ntqg;
                if (isym >= jsym)  {
                    offsets[9 + 8*9 + sym*9 + isym] = ntdg;
                    offsets[9 + 8*9 + sym*9 + jsym] = ntdg;
                xout << "sym, isym, jsym: " << sym << isym << jsym<<"; ntdg=" <<ntdg << std::endl;
                }
                ntqg += this->at(isym) * this->at(jsym);
                if (isym > jsym) {
                    ntdg += this->at(isym) * this->at(jsym);
                } else if (isym == jsym) {
                    ntdg += (this->at(isym) * (this->at(isym)+1))/2;
                }
            }
            offsets[9 + sym*9 + 8] = ntqg;
            offsets[9 + 8*9 + sym*9 + 8] = ntdg;
        }
    }
}

size_t SymmetryOffset::offset(int sym1) const {
    if (offsets.empty()) throw "SymmetryOffset::calculateOffsets has not been called";
    return offsets[sym1];
}

size_t SymmetryOffset::offset(int sym1, int sym2, int parity) const
{
    if (offsets.empty()) throw "SymmetryOffset::calculateOffsets has not been called";
    return offsets[9 + 8*9 * (parity ? 1 : 0) + sym1*9 + sym2];
}

size_t SymmetryOffset::total() const {
    return offset(8);
}

size_t SymmetryOffset::total(int sym1, int parity) const
{
    return offset(sym1,8,parity);
}
