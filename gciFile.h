#ifndef GCIFILE_H
#define GCIFILE_H
#ifdef MOLPRO
#define GCIMOLPROFILE
#endif
#include <vector>
#ifdef GCIMOLPROFILE
#include <cic/ItfFortranInt.h>
#else
#include <fstream>
#endif

namespace gci {

/*!
   * \brief Manage binary I/O on scratch files
 */
class File
{
public:
  /*!
   * \brief Create a new file
   */
  File();
  ~File();

  /*!
   * \brief write data to the file
   * \param buf the buffer to be written
   * \param address offset on the file
   */
  void write(std::vector<double>& buf, size_t address=0);
  /*!
   * \brief read data from the file
   * \param buf the buffer to be read
   * \param address offset on the file
   */
  void read(std::vector<double>& buf, size_t address=0);
private:
#ifdef GCIMOLPROFILE
  static int baseRecord;
  itf::FMolproStorageBlock* f;
#else
  std::fstream f;
#endif
};

}

using namespace gci;

#endif // GCIFILE_H
