#ifndef GCI_GCIVIBOPERATOR_H
#define GCI_GCIVIBOPERATOR_H

#include <molpro/Operator.h>

#include "molpro/gci/gciHProduct.h"
#include "molpro/gci/gciVibExcitation.h"

namespace molpro {
namespace gci {
namespace ns_VibOperator {
enum class parity_t : int {
  none = 0, //!< no symmetry
  even = 1, //!< symmetric
  odd = -1  //!< antisymmetric
};

/*!
 * @brief Unique hash value for this excitation string for a given parity
 * @note Exchange of excitation strings for different modes is alway symmetric.
 * @param parity Symmetry for exchange of modals in excitation string
 * @return
 */
size_t hash(const VibExcitation &exc, int nMode, int nModal, parity_t hermiticity, parity_t exchange);

/*!
 * @brief Hash value for 1 mode operator with conjugation symmetry
 * Order:
 *      1	0	0	0
 *      2	3	0	0
 *      4	5	6 	0
 *      7	8	9	10
 */
size_t hash_mc1_sym(const VibExcitation &exc, int nModal);

/*!
 * @brief Hash value for 1 mode operator with no conjugation symmetry
 * @note this ordering has no value anymore. I thought I could get around needing nModal, but couldn't.
 * Order:
 *      1	4	9	16
 *      2	3	8	15
 *      5	6	7	14
 *      10	11	12	13
 */
size_t hash_mc1_nosym(const VibExcitation &exc, int nModal);

/*!
 * @brief Hash value for 1 mode operator with no conjugation symmetry
 * @note this ordering has no value anymore. I thought I could get around needing nModal, but couldn't.
 * Order:
 *      1	4	9	16
 *      2	3	8	15
 *      5	6	7	14
 *      10	11	12	13
 */
size_t hash_mc1_nosym_old(const VibExcitation &exc, int nModal);
} //  namespace ns_VibOperator

/*!
 * @brief An element of the tensor corresponding to a particular vibrational excitation string
 * @note Assumes 1MC
 * @tparam Container
 */
template <class Container> class VibTensEl {
public:
  using parity_t = ns_VibOperator::parity_t;
  Container oper;    //!< tensor element operator
  VibExcitation exc; //!< excitation

  VibTensEl(const Container &op, VibExcitation vibExc, parity_t hermiticity, parity_t exchange)
      : oper(op), exc(std::move(vibExc)), m_hermiticity(hermiticity), m_exchange(exchange) {}

  auto hermiticity() const { return m_hermiticity; }

  /*!
   * @brief Returns the conjugate transpose of this element including change in phase
   * @todo treat cases with no symmetry
   */
  VibTensEl conjugate() const {
    auto conjOp = Container(oper);
    if (m_hermiticity == parity_t::odd)
      conjOp *= -1;
    auto conjExc = VibExcitation(exc);
    conjExc.conjugate();
    return {conjOp, conjExc, m_hermiticity, m_exchange};
  }

protected:
  parity_t m_hermiticity; //!< conjugation symmetry E_{ij} -> E_{ji}
  parity_t m_exchange;    //!< symmetry under echange of mode indices E_{ij}^A E_{kl}^B -> E_{kl}^B E_{ij}^A
};

/*!
 * @brief Integrals over electronic and vibrational basis, stored with corresponding vibrational excitation operators.
 */
template <class Container> class VibOperator {
public:
  //    using tensor_el_t = std::pair<Container, VibExcitation>;
  using tensor_el_t = VibTensEl<Container>;
  using tensor_t = std::map<size_t, tensor_el_t>;
  using parity_t = ns_VibOperator::parity_t;

  tensor_t tensor;  //!< Mixed tensor
  std::string name; //!< Name of the tensor

  explicit VibOperator(int nMode, int nModal, parity_t hermiticity = parity_t::even, parity_t exchange = parity_t::even,
                       std::string name_ = "")
      : name(std::move(name_)), m_nMode(nMode), m_nModal(nModal), m_hermiticity(hermiticity), m_exchange(exchange) {}

  auto hermiticity() const { return m_hermiticity; }

  /*!
   * @brief Appends a new term to the tensor.
   * @note If a parity equivalent tensor element exists, than nothing is done.
   * @param op Electronic operator
   * @param vibExc Vibrational excitation string
   */
  void append(const Container &op, const VibExcitation &vibExc) {
    auto h = hash(vibExc);
    auto el = tensor_el_t(op, vibExc, m_hermiticity, m_exchange);
    //        auto &&el = VibTensEl<Container>(op, vibExc, m_hermiticity, m_exchange);
    tensor.insert({h, el});
  }

  /*!
   * @brief Access tensor element corresponding to a vibrational excitation string
   * @note Throws std::out_of_range error if element does not exist
   * @param exc Vibrational excitation string
   * @return tensor element
   */
  Container &at(const VibExcitation &exc) {
    auto h = hash(exc);
    return tensor.at(h).oper;
  }

  class VibOperator &operator+=(const VibOperator<Container> &other) {
    for (const auto &otherEl : other.tensor) {
      bool found = false;
      for (auto &thisEl : tensor) {
        if (thisEl.first != otherEl.first)
          continue;
        found = true;
        thisEl.second.oper += otherEl.second.oper;
      }
      if (!found)
        append(otherEl.second.oper, otherEl.second.exc);
    }
    return *this;
  }

  class VibOperator operator+(const VibOperator<Container> &other) const {
    auto newOp = *this;
    newOp += other;
    return newOp;
  }

protected:
  int m_nMode;            //!< number of modes
  int m_nModal;           //!< number of modals per mode (assumed to be the same for each mode)
  parity_t m_hermiticity; //!< conjugation symmetry E_{ij} -> E_{ji}
  parity_t m_exchange;    //!< symmetry under echange of mode indices E_{ij}^A E_{kl}^B -> E_{kl}^B E_{ij}^A

  size_t hash(const VibExcitation &exc) {
    return ns_VibOperator::hash(exc, m_nMode, m_nModal, m_hermiticity, m_exchange);
  }
};

} // namespace gci
} // namespace molpro

#endif // GCI_GCIVIBOPERATOR_H
