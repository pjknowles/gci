#ifndef GCIFILE_H
#define GCIFILE_H
#include <vector>
#include <fstream>

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
  void write(std::vector<double>& buf, size_t address=0);  /*!
   * \brief write data to the file
   * \param buf the buffer to be written
   * \param length size of buffer
   * \param address offset on the file
   */
  void write(double* buf, size_t length, size_t address=0);
  /*!
   * \brief read data from the file
   * \param buf the buffer to be read
   * \param address offset on the file
   */
  void read(std::vector<double>& buf, size_t address=0);
  /*!
   * \brief read data from the file
   * \param buf the buffer to be read
   * \param length size of buffer
   * \param address offset on the file
   */
  void read(double* buf, size_t length, size_t address=0);
  /*!
   * \brief name An optional identifier for the file
   */
  std::string name;
private:
  std::fstream f;
};

}

using namespace gci;

#endif // GCIFILE_H
