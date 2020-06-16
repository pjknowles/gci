#include "gciSymmetrySpace.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <utility>

// using SymmetrySpace = gci::SymmetrySpace;
namespace gci {

SymmetrySpace::SymmetrySpace(std::string title, int maxrank_)
    : std::vector<size_t>(8, 0), maxrank(maxrank_), Title(std::move(title)) {}

std::string SymmetrySpace::str(int verbosity, unsigned int columns) const {
  if (verbosity < 0)
    return std::string("");
  std::ostringstream s;
  s << "SymmetrySpace object " << Title << ", dimensions:";
  for (unsigned long x : *this) {
    s << " " << x;
  }
  if (verbosity > 0 && !offsets.empty()) {
    s << std::endl << "Vector offsets:";
    for (int i = 0; i < 9; i++)
      s << " " << offsets[i];
    s << std::endl << "General matrix offsets:";
    for (int i = 0; i < 9; i++) {
      s << std::endl;
      for (int j = 0; j < 8; j++)
        s << " " << this->offset(j, i, 0);
    }
    s << std::endl << "Symmetric matrix offsets:";
    for (int i = 0; i < 9; i++) {
      s << std::endl;
      for (int j = 0; j < 8; j++)
        s << " " << this->offset(j, i, 1);
    }
    if (verbosity > 1) {
      s << std::endl;
      s << "Raw offsets:";
      for (unsigned long offset : offsets)
        s << " " << offset;
    }
  }
  return s.str();
}

void SymmetrySpace::calculateOffsets() {
  size_t dimension = 9;
  if (maxrank > 1)
    dimension += 8 * 9 * 3;
  offsets.resize(dimension, 0);
  //    xout <<"dimension="<<dimension<<std::endl;
  if (maxrank >= 1) {
    for (int s = 0; s < 8; s++)
      offsets[s + 1] = offsets[s] + this->at(static_cast<unsigned long>(s));
  }
  //    xout << "calculateOffsets sizes="; for (int i=0; i<8; i++) xout << this->at(i) << " "; xout <<std::endl;
  if (maxrank >= 2) {
    for (int sym = 0; sym < 8; sym++) {
      size_t ntqg = 0;
      size_t ntdg = 0;
      size_t ntag = 0;
      for (int isym = 0; isym < 8; isym++) {
        int jsym = sym ^ isym;
        offsets[9 + 8 * 9 + sym * 9 + isym] = ntqg;
        if (isym >= jsym) {
          offsets[9 + 2 * 8 * 9 + sym * 9 + isym] = ntdg;
          offsets[9 + 2 * 8 * 9 + sym * 9 + jsym] = ntdg;
          offsets[9 + 0 * 8 * 9 + sym * 9 + isym] = ntag;
          offsets[9 + 0 * 8 * 9 + sym * 9 + jsym] = ntag;
          //                xout << "sym, isym, jsym: " << sym << isym << jsym<<"; ntdg=" <<ntdg << std::endl;
        }
        ntqg += this->at(isym) * this->at(jsym);
        if (isym > jsym) {
          ntdg += this->at(isym) * this->at(jsym);
          ntag += this->at(isym) * this->at(jsym);
        } else if (isym == jsym) {
          ntdg += (this->at(isym) * (this->at(isym) + 1)) / 2;
          ntag += (this->at(isym) * (this->at(isym) - 1)) / 2;
        }
      }
      offsets[9 + sym * 9 + 8] = ntag;
      offsets[9 + 8 * 9 + sym * 9 + 8] = ntqg;
      offsets[9 + 2 * 8 * 9 + sym * 9 + 8] = ntdg;
    }
  }
}

size_t SymmetrySpace::offset(int sym1) const {
  if (offsets.empty())
    throw std::logic_error("SymmetrySpace::calculateOffsets has not been called");
  return offsets[sym1];
}

size_t SymmetrySpace::offset(int sym1, int sym2, int parity) const {
  if (offsets.empty())
    throw std::logic_error("SymmetrySpace::calculateOffsets has not been called");
  return offsets[9 + 8 * 9 * (std::max(parity, -1) + 1) + sym1 * 9 + sym2];
}

size_t SymmetrySpace::total() const { return offset(8); }

size_t SymmetrySpace::total(int sym1, int parity) const { return offset(sym1, 8, parity); }
} // namespace  gci
