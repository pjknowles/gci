#ifndef PROFILER_H
#define PROFILER_H
#include <string>
#include <map>

/*!
 * \brief The Profiler class: framework for timing code sections
 */
class Profiler
{
public:
  Profiler();
  /*!
   * \brief Profiler construct a named instance
   * \param name the title of this object
   */
  Profiler(std::string name);
  /*!
   * \brief start begin timing a code segment
   * \param name name of the code segment
   */
  void start(const std::string name);
  /*!
   * \brief stop finish timing a code segment
   * \param name if given, must match the argument of the previous start call
   */
  void stop(const std::string name="");
  /*!
   * \brief Generate a printable representation of the object
   * \param verbosity how much to print
   * \return
   */
  std::string str(const int verbosity=0);

  struct times {double cpu; double wall;
                struct Profiler::times& operator+=(const struct Profiler::times &other);
                            struct Profiler::times& operator-=(const struct Profiler::times &other);
               };
private:
  std::string Name;
  std::string current;
  struct times startTimes;
  struct times getTimes();
  std::map<std::string,struct times> results;
};
  std::ostream& operator<<(std::ostream& os, Profiler const& obj);

#endif // PROFILER_H
