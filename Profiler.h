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
   * \brief reset the object
   * \param name the title of this object
   */
  void reset(const std::string name);
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
  std::string str(const int verbosity=0) const;

  struct times {double cpu; double wall;
                struct Profiler::times& operator+=(const struct Profiler::times &other);
                            struct Profiler::times& operator-=(const struct Profiler::times &other);
                            bool operator<(const struct times & b);
               };
  typedef std::map<std::string,struct times> resultMap;
private:
  std::string Name;
  std::string current;
  struct times startTimes;
  struct times getTimes();
  resultMap results;
};
  std::ostream& operator<<(std::ostream& os, Profiler const& obj);
  bool operator<(const Profiler::resultMap::value_type& a, const Profiler::resultMap::value_type& b);

#endif // PROFILER_H
