#ifndef GCISYMMETRY_H
#define GCISYMMETRY_H
#include <string>
#include <vector>

class SymmetryOffset : std::vector<size_t> {
public:
    SymmetryOffset();
    std::string printable();
};

#endif // GCISYMMETRY_H
