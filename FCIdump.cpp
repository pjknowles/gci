#include "FCIdump.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include "gci.h" // for xout

FCIdump::FCIdump(std::string filename)
{
  fileName = filename;
  std::ifstream s;
  s.open(fileName.c_str());
  if ( (s.rdstate() & std::ifstream::failbit ) != 0 ) {
    std::cerr << "Error opening " << fileName <<std::endl;
    throw "FCIDUMP::parameter file missing";
  }
  // cache the namelist data
  std::string ss;
  s >> ss ; namelistData=","; // forget the first word of namelist
  while (s >> ss && ss != "&END" && ss != "/")
    namelistData.append(ss);
  namelistData.append(",DUMMY_KEY=,"); // dummy entry at end to simplify parsing
//  xout <<"namelistData=" <<namelistData <<std::endl;
}

#include "gciState.h"



std::vector<int> FCIdump::parameter(std::string key, std::vector<int> def) { // dirty sucking in from FCIDUMP namelist
  std::vector<int> answer;
  std::vector<std::string> strings = parameter(key,std::vector<std::string>(1," "));
//  xout <<"parameter "<<key<<std::endl;
//  xout << "parameter "<<key<<"="; for (std::vector<std::string>::iterator s=strings.begin(); s < strings.end(); s++) xout <<*s ; xout <<std::endl;
  if (strings == std::vector<std::string>(1," ")) return def;
  for (std::vector<std::string>::const_iterator s1=strings.begin(); s1!= strings.end(); s1++) {
    int i;
    std::istringstream b(*s1);
    b >>i;
//    xout <<"s1="<<*s1<<"; i="<<i<<std::endl;
    answer.push_back(i);
  }
//  xout << "parameter "<<key<<"="; for (std::vector<int>::iterator s=answer.begin(); s < answer.end(); s++) xout <<*s ; xout <<std::endl;
  return answer;
}


std::vector<double> FCIdump::parameter(std::string key, std::vector<double> def) { // dirty sucking in from FCIDUMP namelist
  std::vector<double> answer;
  std::vector<std::string> strings = parameter(key,std::vector<std::string>(1," "));
//  xout <<"parameter "<<key<<std::endl;
//  xout << "parameter "<<key<<"="; for (std::vector<std::string>::iterator s=strings.begin(); s < strings.end(); s++) xout <<*s ; xout <<std::endl;
  if (strings == std::vector<std::string>(1," ")) return def;
  for (std::vector<std::string>::const_iterator s1=strings.begin(); s1!= strings.end(); s1++) {
    double i;
    std::istringstream b(*s1);
    b >>i;
//    xout <<"s1="<<*s1<<"; i="<<i<<std::endl;
    answer.push_back(i);
  }
//  xout << "parameter "<<key<<"="; for (std::vector<double>::iterator s=answer.begin(); s < answer.end(); s++) xout <<*s ; xout <<std::endl;
  return answer;
}

std::vector<std::string> FCIdump::parameter(std::string key, std::vector<std::string> def) { // dirty sucking in from FCIDUMP namelist
  std::vector<std::string> answer;
//  xout <<"search for "<<","<<key<<"="<<std::endl;
  size_t pos = namelistData.find(","+key+"=");
  if (pos == std::string::npos) return def;
  pos=namelistData.find("=",pos)+1;
  size_t pose=namelistData.find("=",pos);
  pose=namelistData.find_last_of(",",pose)+1;
//  xout << "string parameter key="<<key<<"; values='"<<namelistData.substr(pos,pose-pos)<<"'"<<std::endl;
//  xout << "pos="<<pos<<"; pose="<<pose<<std::endl;
  for (size_t posNext=pos; posNext < pose; pos=posNext) {
    posNext=namelistData.find(",",pos)+1;
    answer.push_back(namelistData.substr(pos,posNext-pos-1));
  }
//  xout <<"answer";
//  for (std::vector<std::string>::const_iterator s=answer.begin(); s != answer.end(); s++)
//    xout <<":"<<*s;
//  xout <<std::endl;
  return answer;
}
