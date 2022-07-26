#ifndef GCI_GCIPERSISTENTOPERATOR_H
#define GCI_GCIPERSISTENTOPERATOR_H

#ifdef HAVE_HDF5
#include <hdf5.h>
#else // HAVE_HDF5
#define hid_t int
#endif // HAVE_HDF5
#include <memory>
#include <molpro/symmetry_matrix/Operator.h>

namespace molpro {
namespace gci {

/*!
 * @brief stores SymmetryMatrix::Operator either in memory for small operators or on disk in an hdf5 file.
 * @note the hdf5 file should have already been created and should be open whenever ``PersistentOperator``
 * is used
 */
class PersistentOperator {
protected:
  std::shared_ptr<molpro::Operator> m_op; //!< small operators are stored in memory
  hid_t m_file_id;           //!< hdf5 file must already be open and remain open while access to operator is needed
  std::string m_description; //!< description of the operator, used to idenitfy it in hdf5 file
  bool m_on_disk;            //!< whether the operator was saved to disk
public:
  PersistentOperator();
  /*!
   * @brief takes ownership of the operator, writing it to hdf5 if it's large or keeping it in memory
   *
   * @note only root process writes to disk, other process participate in updating structure of hdf5 file
   *
   * @param op pointer to the operator that will be stored
   * @param _description unique description of the operator
   * @param root rank of the root process
   * @param id id of hdf5 file
   * @param in_memory flags operator to be stored in memory, irrespective of its size
   */
  PersistentOperator(const std::shared_ptr<molpro::Operator> &op, std::string _description, int root, hid_t id,
                     bool in_memory = false);
  //! operator is stored on a previously created hdf5 file
  PersistentOperator(hid_t id, std::string _description);
  //! gets a copy of the stored operator
  std::shared_ptr<molpro::Operator> get() const;
  //! unique description of the operator
  std::string description();

protected:
  //! converts operator to bytestream and writes it to hdf5 file
  void store_bytestream(const std::shared_ptr<molpro::Operator> &op, int root);
};

} // namespace gci
} // namespace molpro
#endif // GCI_GCIPERSISTENTOPERATOR_H
