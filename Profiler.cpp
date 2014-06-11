#include "Profiler.h"
#include <sstream>
#include <iostream>

Profiler::Profiler()
{
}

Profiler::Profiler(std::string name)
{
  reset(name);
}

void Profiler::reset(std::string name)
{
  Name=name;
}

struct Profiler::times& Profiler::times::operator+=( const struct Profiler::times &w2)
{
  cpu += w2.cpu;
  wall += w2.wall;
  return *this;
}
struct Profiler::times& Profiler::times::operator-=( const struct Profiler::times &w2)
{
  cpu -= w2.cpu;
  wall -= w2.wall;
  return *this;
}

struct Profiler::times operator+(const struct Profiler::times &w1, const struct Profiler::times &w2)
{
  struct Profiler::times result=w1;
  result += w2;
  return result;
}
struct Profiler::times operator-(const struct Profiler::times &w1, const struct Profiler::times &w2)
{
  struct Profiler::times result=w1;
  result -= w2;
  return result;
}


void Profiler::start(std::string name)
{
  current=name;
  startTimes=getTimes();
}

#include <assert.h>
void Profiler::stop(std::string name)
{
  assert(name=="" || name == current);
  this->results[current] += getTimes()-startTimes;
  startTimes=getTimes();
}

std::string Profiler::str(const int verbosity) const
{
  if (verbosity<0) return "";
  std::stringstream ss;
  ss << "Profiler "<<Name<<std::endl;
  for (std::map<std::string,struct times>::const_iterator result=results.begin(); result!=results.end(); ++result)
    ss << (*result).first <<" "<<(*result).second.cpu<<" "<<(*result).second.wall<<std::endl;
  return ss.str();
}
std::ostream& operator<<(std::ostream& os, Profiler const& obj)
{
  return os << obj.str();
}


#include <sys/time.h>
struct Profiler::times Profiler::getTimes()
{
  struct Profiler::times result;
  result.cpu=(double)clock()/CLOCKS_PER_SEC;
  struct timeval time;
  result.wall=(double)0;
  if (!gettimeofday(&time,NULL))
    result.wall = (double)time.tv_sec + (double)time.tv_usec * .000001;
  return result;
}
