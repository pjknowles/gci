#include "gciPersistentOperator.h"
#include "gciUtils.h"
#include "gci.h"

namespace gci {

PersistentOperator::PersistentOperator(std::shared_ptr<SymmetryMatrix::Operator> &op, std::string fname_hdf5)
        : m_fname_hdf5(std::move(fname_hdf5)), m_description(op->m_description), m_on_disk(false) {
    if (m_fname_hdf5.empty()) throw std::runtime_error("PersistentOperator::PersistentOperator(): empty file name");
    if (op->m_rank < 2)
        m_op = op;
    else {
        store_bytestream(*op);
        m_on_disk = true;
    }
}

std::string PersistentOperator::description() {return m_description;}

PersistentOperator::PersistentOperator(std::string fname_hdf5, std::string _description)
        : m_fname_hdf5(std::move(fname_hdf5)), m_description(std::move(_description)), m_on_disk(true) { }

void PersistentOperator::store_bytestream(const SymmetryMatrix::Operator &op_) {
    bool create = !utils::file_exists(m_fname_hdf5);
    auto id_f = gci::utils::open_hdf5_file(m_fname_hdf5, gci::mpi_comm_compute, create);
    auto bs = op_.bytestream();
    auto dataset = gci::utils::open_or_create_hdf5_dataset(id_f, m_description, H5T_NATIVE_CHAR, bs.size());
    auto plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    auto status = H5Dwrite(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, plist_id, &bs.data().front());
    if (status < 0) throw std::runtime_error("PersistentOperator::store_bytestream(): write failed");
    H5Fclose(plist_id);
    H5Fclose(dataset);
    H5Fclose(id_f);
}

std::shared_ptr<SymmetryMatrix::Operator> PersistentOperator::get() {
    if (!m_on_disk) return m_op;
    bool create = !utils::file_exists(m_fname_hdf5);
    auto id_f = gci::utils::open_hdf5_file(m_fname_hdf5, gci::mpi_comm_compute, create);
    auto buffer = std::vector<char>();
    // get size of the dataset and allocate a buffer
    auto dataset = H5Dopen(id_f, m_description.c_str(), H5P_DEFAULT);
    auto dspace = H5Dget_space(dataset);
    auto ndims = H5Sget_simple_extent_ndims(dspace);
    if (ndims != 1) throw std::runtime_error("PersistentOperator::get(): ndims != 1");
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, nullptr);
    buffer.resize(dims[0]);
    auto plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    auto status = H5Dread(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, plist_id, buffer.data());
    if (status < 0) throw std::runtime_error("PersistentOperator::store_bytestream(): write failed");
    H5Fclose(plist_id);
    H5Fclose(dspace);
    H5Fclose(dataset);
    H5Fclose(id_f);
    auto op = std::make_shared<SymmetryMatrix::Operator>(SymmetryMatrix::Operator::construct(buffer.data()));
    return op;
}

} // namespace gci
