#include <algorithm>
#include <sstream>
#include <iostream>
#include "Profiler.h"

Profiler::Profiler()
{
}

Profiler::Profiler(std::string name)
{
  reset(name);
}

void Profiler::reset(const std::string name)
{
  Name=name;
  results.clear();
  stop();
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


void Profiler::start(const std::string name)
{
  if (current!="") this->results[current] += getTimes()-startTimes;
  current=name;
  startTimes=getTimes();
}

#include <assert.h>
void Profiler::stop(const std::string name)
{
  assert(name=="" || name == current);
  start("* Other");
}

std::string Profiler::str(const int verbosity) const
{
  std::sort(results.begin(),results.end());
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

bool operator<(const struct Profiler::times & a, const struct Profiler::times & b) const
{
  return a.wall < b.wall;
}

bool operator<( const Profiler::resultMap::value_type& a, const Profiler::resultMap::value_type& b) const
{
  return a.second < b.second;
}
