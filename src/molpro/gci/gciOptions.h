#ifndef GCIOPTIONS_H
#define GCIOPTIONS_H
#include <molpro/Options.h>
namespace molpro::gci {
class Options : public molpro::Options {
public:
  explicit Options(std::string input = "") : molpro::Options("GCI", std::move(input)) {}
};
} // namespace molpro::gci

#endif // GCIOPTIONS_H
