#include "gciPersistentOperator.h"
#include "gci.h"
#include "gciUtils.h"

namespace gci {
PersistentOperator::PersistentOperator() : m_file_id(-1), m_on_disk(false) {}

PersistentOperator::PersistentOperator(const std::shared_ptr<molpro::Operator> &op, std::string _description, int root,
                                       hid_t id, bool in_memory)
    : m_file_id(id), m_description(std::move(_description)), m_on_disk(false) {
  auto rank = op == nullptr ? 2 : op->m_rank;
  if (rank < 2 || in_memory) {
    m_op = op;
  } else {
    if (!utils::hdf5_file_open(id))
      throw std::runtime_error("PersistentOperator::PersistentOperator(): hdf5 file is not open");
    store_bytestream(op, root);
    m_on_disk = true;
  }
}

std::string PersistentOperator::description() { return m_description; }

PersistentOperator::PersistentOperator(hid_t id, std::string _description)
    : m_file_id(id), m_description(std::move(_description)), m_on_disk(true) {}

void PersistentOperator::store_bytestream(const std::shared_ptr<molpro::Operator> &op, int root) {
  hid_t dataset;
  if (gci::parallel_rank == root) {
    if (op == nullptr)
      throw std::runtime_error("PersistentOperator::store_bytestream() *operator is nullptr");
    auto bs = op->bytestream();
    unsigned long size = bs.size();
    MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, root, gci::mpi_comm_compute);
    dataset = utils::open_or_create_hdf5_dataset(m_file_id, m_description, H5T_NATIVE_CHAR, bs.size());
    auto plist = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist, H5FD_MPIO_INDEPENDENT);
    auto status = H5Dwrite(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, plist, &bs.data().front());
    if (status < 0)
      throw std::runtime_error("PersistentOperator::store_bytestream(): write failed");
    H5Pclose(plist);
  } else {
    unsigned long size;
    MPI_Bcast(&size, 1, MPI_UNSIGNED_LONG, root, gci::mpi_comm_compute);
    dataset = utils::open_or_create_hdf5_dataset(m_file_id, m_description, H5T_NATIVE_CHAR, size);
  }
  H5Dclose(dataset);
}

std::shared_ptr<molpro::Operator> PersistentOperator::get() const {
  if (!m_on_disk)
    return m_op;
  auto buffer = std::vector<char>();
  // get size of the dataset and allocate a buffer
  auto dataset = H5Dopen(m_file_id, m_description.c_str(), H5P_DEFAULT);
  auto dspace = H5Dget_space(dataset);
  auto ndims = H5Sget_simple_extent_ndims(dspace);
  if (ndims != 1)
    throw std::runtime_error("PersistentOperator::get(): ndims != 1");
  hsize_t dims[ndims];
  H5Sget_simple_extent_dims(dspace, dims, nullptr);
  buffer.resize(dims[0]);
  auto plist_id = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
  auto status = H5Dread(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, plist_id, buffer.data());
  if (status < 0)
    throw std::runtime_error("PersistentOperator::store_bytestream(): write failed");
  H5Sclose(dspace);
  H5Dclose(dataset);
  H5Pclose(plist_id);
  auto op = std::make_shared<molpro::Operator>(molpro::Operator::construct(buffer.data()));
  return op;
}

} // namespace gci
