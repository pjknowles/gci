#include "gciFile.h"
#include <unistd.h>
#include <molpro/linalg/array/util/temp_file.h>

// using File = gci::File;
namespace molpro {
namespace gci {

File::File() : f(molpro::linalg::array::util::temp_file_name("tmp", "tmp"),
                 std::fstream::binary | std::fstream::in | std::fstream::out) {
}

void File::read(std::vector<double> &buf, size_t address) { read(&buf[0], buf.size(), address); }

void File::read(double *buf, size_t length, size_t address) {
  f.seekg(static_cast<long long int>(address * 8), std::ios_base::beg);
  f.read((char *) &buf[0], length * 8);
}

void File::write(std::vector<double> &buf, size_t address) { write(&buf[0], buf.size(), address); }

void File::write(double *buf, size_t length, size_t address) {
  f.seekp(static_cast<long long int>(address * 8), std::ios_base::beg);
  f.write((char *) &buf[0], length * 8);
}
} // namespace gci
} // namespace molpro
