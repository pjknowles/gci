#include "SMat.h"
void default_get_orbital_space(char c, size_t nt[]) {}
static void (*get_orbital_space)(char c, size_t nt[]) = &default_get_orbital_space;
void molpro::SymmetryMatrix::register_get_orbital_space(void (*func)(char c, size_t nt[])) { get_orbital_space = func; }
molpro::dims_t molpro::SymmetryMatrix::spaces(std::string space) {
  dims_t result;
  dim_t nt1(8);
  for (const auto& s : space) {
    get_orbital_space(s, &(nt1[0]));
    result.push_back(nt1);
  }
  return result;
}
