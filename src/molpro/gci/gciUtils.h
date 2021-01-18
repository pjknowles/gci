#ifndef GCI_GCIUTILS_H
#define GCI_GCIUTILS_H

#ifdef HAVE_HDF5
#include <hdf5.h>
#else // HAVE_HDF5
#define hid_t int
#endif // HAVE_HDF5
#include <mpi.h>
#include <string>

namespace molpro {
namespace gci {
namespace utils {

//! checks if a file exists. Optionally, writes a message to std::out
bool file_exists(const std::string &fname, const std::string &message = "");

//! opens or creates an hdf5 file
hid_t open_hdf5_file(const std::string &fname, MPI_Comm communicator, bool create);

//! checks if hdf5 file is open using its id
bool hdf5_file_open(hid_t file_id);

bool hdf5_dataset_exists(hid_t location, const std::string &dataset_name);

//! opens or creates a dataset in hdf5 file
hid_t open_or_create_hdf5_dataset(const hid_t &location, const std::string &dataset_name, const hid_t &dtype_id,
                                  const size_t &length);
} // namespace utils
} // namespace gci
} // namespace molpro
#endif // GCI_GCIUTILS_H
