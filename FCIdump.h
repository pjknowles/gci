#ifndef FCIDUMP_H
#define FCIDUMP_H
#include <string>
#include <vector>

/*!
 * \brief Class that provides access to FCIdump files
 */
class FCIdump
{
public:
  /*!
     * \brief Construct FCIdump object
     * \param filename The file containing the FCIDUMP data
     */
  FCIdump(std::string filename="FCIDUMP");

  /*!
     * \brief parameter Obtain an integer namelist parameter from the FCIDUMP data.
     * \param key The name of the parameter
     * \param def Default value if the parameter is not found.
     * \return  The result as a vector of integers.
     */
  std::vector<int> parameter(std::string key, std::vector<int> def=std::vector<int>(1,0));

  /*!
     * \brief parameter Obtain a real namelist parameter from the FCIDUMP data.
     * \param key The name of the parameter
     * \param def Default value if the parameter is not found.
     * \return  The result as a vector of integers.
     */
  std::vector<double> parameter(std::string key, std::vector<double> def);

  /*!
     * \brief parameter Obtain a string namelist parameter from the FCIDUMP data.
     * \param key The name of the parameter
     * \param def Default value if the parameter is not found.
     * \return  The result as a vector of integers.
     */
  std::vector<std::string> parameter(std::string key, std::vector<std::string> def);
  /*!
   * \brief addParameter add a parameter
   * \param key key
   * \param values values
   */
  void addParameter(const std::string& key, const std::vector<std::string>& values);

   /*!
     * \brief fileName The file containing the FCIDUMP data
     */
  std::string fileName;
private:
  std::string namelistData;
};

#endif // FCIDUMP_H
