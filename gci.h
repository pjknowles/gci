#ifndef GCI_H
#define GCI_H
#include "Profiler.h"

#ifndef MOLPRO
#define xout std::cout
#if defined(GA_MPI) || defined(GA_MPI2)
#define GCI_PARALLEL
#endif
#else
#include "gciMolpro.h"
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
  if (option < 0) __nextval_counter=-2;
  return ++__nextval_counter;
#endif
}

// task scheduler
#ifdef MOLPRO
extern itf::FMppInt mpp;
#else
extern int64_t __task_granularity, __task, __my_first_task;
#endif
inline void DivideTasks(std::size_t ntasks, std::size_t nMinBatch = 0, std::size_t nMaxBatch = 0) {{
#ifdef MOLPRO
  mpp.DivideTasks(ntasks, nMinBatch, nMaxBatch);
#else
  // simple static LB
  __task_granularity = (ntasks-1)/parallel_size+1;
  __task_granularity = __task_granularity > 0 ? __task_granularity : 1;
  if ((int64_t)nMinBatch > 0 && __task_granularity < (int64_t)nMinBatch) __task_granularity = (int64_t)nMinBatch;
  if ((int64_t)nMaxBatch > 0 && __task_granularity > (int64_t)nMaxBatch) __task_granularity = (int64_t)nMaxBatch;
//  __task_granularity = 3;
  nextval(-parallel_size);
  __task=0;
  __my_first_task=-__task_granularity-1;
#endif
}}
inline bool NextTask() {
  {
#ifdef MOLPRO
//    return (parallel_rank==0);
  return mpp.NextTask();
#else
  if (__my_first_task+__task_granularity <= __task)
    __my_first_task = nextval()*__task_granularity;
  return (__task++ >= __my_first_task && __task <= __my_first_task+__task_granularity) ;
#endif
  }
}
inline void EndTasks() {
#ifdef MOLPRO
  mpp.EndTasks();
#else
#endif
}


}


#endif // GCI_H
