#ifndef GCI_H
#define GCI_H
#include "Profiler.h"

#ifndef MOLPRO
#define xout std::cout
#if defined(GA_MPI) || defined(GA_MPI2)
#define GCI_PARALLEL
#endif
#else
#include <cic/ItfCommon.h>
#include <cic/ItfFortranInt.h>
#endif

#ifdef GCI_PARALLEL
extern "C" {
#include "ppidd_c.h"
}
#endif


#include <iostream>
namespace gci {

extern Profiler profiler; // global profiler

extern int64_t parallel_rank, parallel_size;

// shared counter
extern int64_t __nextval_counter;
inline long nextval(int64_t option=parallel_size){
#ifdef GCI_PARALLEL
  int64_t value; PPIDD_Nxtval(&option,&value); //xout <<std::endl<<"@nextval("<<option<<",rank="<<parallel_rank<<")="<<value<<std::endl;
  return value;
#else
  if (option < 0) __nextval_counter=-1;
  return ++__nextval_counter;
#endif
}

// task scheduler
extern int64_t __task_granularity, __task, __my_first_task;
inline void initask(int64_t granularity) {
  __task_granularity = granularity;
  nextval(-parallel_size);
  __task=0;
  __my_first_task=-granularity-1;
}
inline bool mytask() {
  if (__my_first_task+__task_granularity <= __task)
    __my_first_task = nextval()*__task_granularity;
  return (__task++ >= __my_first_task && __task <= __my_first_task+__task_granularity) ;
}

}


#endif // GCI_H
