#include <algorithm>
#include <sstream>
#include <iostream>
#include <deque>
#include <queue>
#include <iomanip>
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
  if (verbosity<0) return "";
  typedef std::pair<std::string,Profiler::times> data_t;
  std::priority_queue<data_t, std::deque<data_t>, compareTimes<data_t>  > q(results.begin(),results.end());
  std::stringstream ss;
  size_t maxWidth=0;
  Profiler::times totalTimes;
  for (resultMap::const_iterator s=results.begin(); s!=results.end(); ++s) {
      if ((*s).first.size() > maxWidth) maxWidth=(*s).first.size();
      totalTimes += (*s).second;
  }
  q.push(data_t("* TOTAL",totalTimes));
  ss << "Profiler "<<Name<<std::endl;
  while (! q.empty()) {
    ss.precision(2);
    ss <<std::right <<std::setw(maxWidth) << q.top().first <<": cpu="<<std::fixed<<q.top().second.cpu<<", wall="<<q.top().second.wall<<std::endl;
    q.pop();
  }
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
