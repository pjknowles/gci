#include "gciPersistentOperator.h"
#include "gciUtils.h"
#include "gci.h"

namespace gci {

PersistentOperator::PersistentOperator(std::shared_ptr<SymmetryMatrix::Operator> &op, hid_t id)
        : m_file_id(id), m_description(op->m_description), m_on_disk(false) {
    if (op->m_rank < 2)
        m_op = op;
    else {
        if (!utils::hdf5_file_open(id))
            throw std::runtime_error("PersistentOperator::PersistentOperator(): hdf5 file is not open");
        store_bytestream(*op);
        m_on_disk = true;
    }
}

PersistentOperator::PersistentOperator(std::string _description, size_t size, hid_t id)
        : m_file_id(id), m_description(std::move(_description)), m_on_disk(false) {
    if (!utils::hdf5_file_open(id))
        throw std::runtime_error("PersistentOperator::PersistentOperator(): hdf5 file is not open");
    store_bytestream(_description, size);
    m_on_disk = true;
}

std::string PersistentOperator::description() {return m_description;}

PersistentOperator::PersistentOperator(hid_t id, std::string _description)
        : m_file_id(id), m_description(std::move(_description)), m_on_disk(true) { }

void PersistentOperator::store_bytestream(const SymmetryMatrix::Operator &op_) {
    auto bs = op_.bytestream();
    auto dataset = utils::open_or_create_hdf5_dataset(m_file_id, m_description, H5T_NATIVE_CHAR, bs.size());
    auto plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
    auto status = H5Dwrite(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, plist_id, &bs.data().front());
    if (status < 0) throw std::runtime_error("PersistentOperator::store_bytestream(): write failed");
    H5Dclose(dataset);
    H5Pclose(plist_id);
}

void PersistentOperator::store_bytestream(const std::string &_description, size_t size) {
    auto dataset = utils::open_or_create_hdf5_dataset(m_file_id, m_description, H5T_NATIVE_CHAR, size);
    H5Dclose(dataset);
}

std::shared_ptr<SymmetryMatrix::Operator> PersistentOperator::get() {
    if (!m_on_disk) return m_op;
    auto buffer = std::vector<char>();
    // get size of the dataset and allocate a buffer
    auto dataset = H5Dopen(m_file_id, m_description.c_str(), H5P_DEFAULT);
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
    H5Sclose(dspace);
    H5Dclose(dataset);
    H5Pclose(plist_id);
    auto op = std::make_shared<SymmetryMatrix::Operator>(SymmetryMatrix::Operator::construct(buffer.data()));
    return op;
}

} // namespace gci
