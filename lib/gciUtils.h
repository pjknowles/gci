#ifndef GCI_GCIUTILS_H
#define GCI_GCIUTILS_H

#include <string>
#include <hdf5.h>
#include <mpi.h>

namespace gci {
namespace utils {
hid_t open_hdf5_file(const std::string &fname, MPI_Comm communicator, bool create);

hid_t open_or_create_hdf5_dataset(const hid_t &location, const std::string &dataset_name, const hid_t &dtype_id,
                                  const size_t &length);
} // namespace utils
} // namespace gci
#endif //GCI_GCIUTILS_H
