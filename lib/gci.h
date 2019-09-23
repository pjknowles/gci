#ifndef GCI_H
#define GCI_H

#include <memory>
#include <vector>
#include <cstring>
#include <stdexcept>
#include "Profiler.h"
#include <cstdint>
#include "sharedCounter.h"
#include "memory.h"

#ifndef MOLPRO
#define xout std::cout
#else
#include "molpro_config.h"
#include "gciMolpro.h"
#include "ppidd.h"
#ifdef HAVE_MPI_H
#include <mpi.h>
#define MPI_COMM_COMPUTE MPI_Comm_f2c(PPIDD_Worker_comm())
#else
#define MPI_COMM_NULL 0
#endif
#endif
#ifdef HAVE_MPI_H

#include <cstdint>
#include <climits>

#if SIZE_MAX == UCHAR_MAX
#define MPI_SIZE_T MPI_UNSIGNED_CHAR
#elif SIZE_MAX == USHRT_MAX
#define MPI_SIZE_T MPI_UNSIGNED_SHORT
#elif SIZE_MAX == UINT_MAX
#define MPI_SIZE_T MPI_UNSIGNED
#elif SIZE_MAX == ULONG_MAX
#define MPI_SIZE_T MPI_UNSIGNED_LONG
#elif SIZE_MAX == ULLONG_MAX
#define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#else
#error "what is happening here?"
#endif
#endif

#include <iostream>
#include <memory>

namespace gci {

extern std::unique_ptr<Profiler> profiler; // global profiler

extern int parallel_rank, parallel_size;

extern MPI_Comm molpro_plugin_intercomm;
extern bool molpro_plugin;

// shared counter
//extern std::unique_ptr<sharedCounter> _nextval_counter;
extern std::map<MPI_Comm, std::unique_ptr<sharedCounter>> _nextval_counter;

void create_new_counter(MPI_Comm communicator) {
    if (_nextval_counter.find(communicator) == _nextval_counter.end()) {
        _nextval_counter[communicator] = std::make_unique<sharedCounter>(sharedCounter());
    }
}

inline long nextval(MPI_Comm communicator, int64_t option = parallel_size) {
    if (option < 0) {
        if (_nextval_counter[communicator] == nullptr)
            _nextval_counter[communicator] = std::make_unique<sharedCounter>();
        _nextval_counter[communicator]->reset();
        return 0;
    }
    return _nextval_counter[communicator]->increment();
}

// task scheduler
extern std::map<MPI_Comm, long int> __task_granularity, __task, __my_first_task;

inline void DivideTasks(std::size_t ntasks, std::size_t nMinBatch, std::size_t nMaxBatch,
                        MPI_Comm communicator) {
    {
        // simple static LB
        size_t task_gran = ((ntasks - 1) / parallel_size + 1);
        task_gran = task_gran > 0 ? task_gran : 1;
        if (nMinBatch > 0 && task_gran < nMinBatch)
            task_gran = nMinBatch;
        if (nMaxBatch > 0 && task_gran > nMaxBatch)
            task_gran = nMaxBatch;
        __task_granularity[communicator] = task_gran;
        nextval(communicator, -parallel_size);
        __task[communicator] = 0;
        __my_first_task[communicator] = -__task_granularity[communicator] - 1;
    }
}

inline bool NextTask(MPI_Comm communicator) {
    {
        if (__my_first_task[communicator] + __task_granularity[communicator] <= __task[communicator])
            __my_first_task[communicator] = nextval(communicator) * __task_granularity[communicator];
        return (__task[communicator]++ >= __my_first_task[communicator] &&
                __task[communicator] <= __my_first_task[communicator] + __task_granularity[communicator]);
    }
}

inline void EndTasks(MPI_Comm communicator) {
}

// simple task distribution assembly

inline void gather_chunks(double *buffer, const size_t length, const size_t chunk, MPI_Comm communicator) {
    {
        std::vector<int> recvcounts(static_cast<unsigned long>(parallel_size)),
                displs(static_cast<unsigned long>(parallel_size));
        displs[0] = 0;
        for (size_t i = 1; i < (size_t) parallel_size; i++) {
            displs[i] = static_cast<int>((int) i * chunk);
            if (displs[i] > (int) (length)) displs[i] = (int) length;
            recvcounts[i - 1] = displs[i] - displs[i - 1];
        }
        recvcounts[parallel_size - 1] = (int) (length) - displs[parallel_size - 1];
//        xout << "nsa="<<nsa<<std::endl;
//        xout << "displ:"; for (size_t i=0; i<(size_t)parallel_size; i++) xout <<" "<<displs[i]; xout <<std::endl;
//        xout << "recvcounts:"; for (size_t i=0; i<(size_t)parallel_size; i++) xout <<" "<<recvcounts[i]; xout <<std::endl;
#ifdef HAVE_MPI_H
        MPI_Allgatherv(MPI_IN_PLACE,
                       0,
                       MPI_DATATYPE_NULL,
                       buffer,
                       &recvcounts[0],
                       &displs[0],
                       MPI_DOUBLE,
                       communicator);
#endif
    }
}

void inline gsum(double *buffer, size_t len, MPI_Comm communicator) {
#ifdef HAVE_MPI_H
    std::vector<double> result;
    if (parallel_rank == 0)
        result.resize(len);
    MPI_Reduce(buffer, &result[0], (int) len, MPI_DOUBLE, MPI_SUM, 0, communicator);
    if (parallel_rank == 0)
        std::memcpy(buffer, &result[0], sizeof(double) * len);
    MPI_Bcast(buffer, (int) len, MPI_DOUBLE, 0, communicator);
#endif
}

void inline gsum(std::map<size_t, double> &buffer, MPI_Comm communicator) {
#ifdef HAVE_MPI_H
    std::map<size_t, double> result;
    for (int rank = 0; rank < parallel_size; rank++) {
        size_t siz = buffer.size();
        MPI_Bcast(&siz, (int) 1, MPI_SIZE_T, rank, communicator);
        std::vector<size_t> addresses(siz);
        std::vector<double> values(siz);
        if (rank == parallel_rank) {
            size_t i = 0;
            for (const auto &b : buffer) {
                addresses[i] = b.first;
                values[i] = b.second;
                ++i;
            }
        }
        MPI_Bcast(addresses.data(), static_cast<int>(siz), MPI_SIZE_T, rank, communicator);
        MPI_Bcast(values.data(), static_cast<int>(siz), MPI_DOUBLE, rank, communicator);
        for (size_t i = 0; i < siz; i++) {
            result[addresses[i]] += values[i];
        }
    }
    buffer = result;
#endif
}

void inline gsum(memory::array<double> v, MPI_Comm communicator) {gsum(&v[0], v.size(), communicator);}

void inline gsum(memory::vector<double> v, MPI_Comm communicator) {gsum(&v[0], v.size(), communicator);}

void inline gsum(std::vector<double> v, MPI_Comm communicator) {gsum(&v[0], v.size(), communicator);}

}

#endif // GCI_H
