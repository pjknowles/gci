#ifndef GCI_PARALLEL_UTILS_H
#define GCI_PARALLEL_UTILS_H

namespace gci {
class Lock {
public:
    explicit Lock(int mutex = 0) : mutex(mutex) {
        GA_Lock(mutex);
    }

    ~Lock() {GA_Unlock(mutex);}

    int mutex;
};

} // namespace gci
#endif //GCI_PARALLEL_UTILS_H
