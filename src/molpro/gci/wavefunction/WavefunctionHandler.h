#ifndef GCI_SRC_MOLPRO_GCI_WAVEFUNCTION_WAVEFUNCTIONHANDLER_H
#define GCI_SRC_MOLPRO_GCI_WAVEFUNCTION_WAVEFUNCTIONHANDLER_H
#include <molpro/linalg/array/ArrayHandler.h>

#include <map>

namespace molpro {
namespace gci {
namespace wavefunction {
template <class AL, class AR>
class WavefunctionHandler : public linalg::array::ArrayHandler<AL, AR> {
public:
  using typename linalg::array::ArrayHandler<AL, AR>::value_type_L;
  using typename linalg::array::ArrayHandler<AL, AR>::value_type_R;
  using typename linalg::array::ArrayHandler<AL, AR>::value_type;
  using typename linalg::array::ArrayHandler<AL, AR>::value_type_abs;
  using typename linalg::array::ArrayHandler<AL, AR>::ProxyHandle;

  ProxyHandle lazy_handle() override { return this->lazy_handle(*this); };

  using linalg::array::ArrayHandler<AL, AR>::lazy_handle;

  AL copy(const AR &source) override { return AL{source}; };

  void scal(value_type alpha, AL &x) override { x.scal(alpha); }

  void fill(value_type alpha, AL &x) override { x.fill(alpha); }

  void axpy(value_type alpha, const AR &x, AL &y) override { y.axpy(alpha, x); }

  value_type dot(const AL &x, const AR &y) override { return x.dot(y); }

  std::map<size_t, value_type_abs> select_max_dot(size_t n, const AL &x, const AR &y) override {
    return x.distr_buffer.select_max_dot(n, y.distr_buffer);
  }

protected:
  using linalg::array::ArrayHandler<AL, AR>::error;
};

} // namespace wavefunction
} // namespace gci
} // namespace molpro

#endif // GCI_SRC_MOLPRO_GCI_WAVEFUNCTION_WAVEFUNCTIONHANDLER_H
