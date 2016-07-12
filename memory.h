#ifndef MEMORY_H
#define MEMORY_H
#include <stddef.h>
#include <stdint.h>
#include <iostream>
#ifdef __cplusplus
extern "C" {
#endif
  /*!
   * \brief memory_allocate Create a new array on the stack
   * \param n The number of bytes to allocate
   * \return pointer to the created array
   */
  double* memory_allocate(size_t n);
  /*!
   * \brief memory_release Release the memory allocated to a given pointer.
   * As a side-effect, the top of the stack is reset to be just below
   * the given pointer, thereby releasing memory allocated to other
   * pointers allocated later than the given one.
   * \param p
   */
  void memory_release(void* p);
  /*!
   * \brief memory_save Remember the state of the stack, with a view
   * to coming back to it via a call to \c memory_release_saved
   * \return
   */
  size_t memory_save();
  /*!
   * \brief memory_release_saved Release all stack memory allocated
   * since call to \ref memory_save
   * \param p The handle returned by a previous \ref memory_save
   */
  void memory_release_saved(size_t p);
  /*!
   * \brief memory_used Report used memory
   * \param maximum if true (default is false), report the maximum
   * memory allocated to date, in bytes.
   * \return
   */
  size_t memory_used(int maximum);
  /*!
   * \brief memory_remaining Report used memory
   * \return remaining memory in bytes
   */
  size_t memory_remaining();
  /*!
   * \brief memory_reset_maximum_stack  Reset the maximum stack used
   * statistic to the currently-used stack
   * \param level if non-negative, the stack position desired; default is
   * current stack size in bytes
   */
  void memory_reset_maximum_stack(int64_t level);

}
#ifdef __cplusplus

namespace memory {
/*!
   * An allocator suitable for use with STL, that uses Molpro's memory stack
   */
  template<typename T> class Allocator {
  public :
    typedef T value_type;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef std::size_t size_type;
    typedef std::ptrdiff_t difference_type;

  public :
    template<typename U>
      struct rebind {
        typedef Allocator<U> other;
      };

  public :
    inline explicit Allocator() {}
    inline ~Allocator() {}
    inline explicit Allocator(Allocator const&) {}
    template<typename U>
      inline explicit Allocator(Allocator<U> const&) {}

    inline pointer address(reference r) { return &r; }
    inline const_pointer address(const_reference r) { return &r; }

    inline pointer allocate(size_type cnt,
			    typename std::allocator<void>::const_pointer = 0) {
      return reinterpret_cast<pointer>(memory_allocate(cnt * sizeof (T)));
    }
    inline void deallocate(pointer p, size_type) {
      memory_release(p);
    }

    inline size_type max_size() const {
      return memory_remaining() / sizeof(T);
    }

    inline void construct(pointer p, const T& t) { new(p) T(t); }
    inline void destroy(pointer p) { p->~T(); }

    inline bool operator==(Allocator const&) { return true; }
    inline bool operator!=(Allocator const& a) { return !operator==(a); }
  };
}
inline void memory_module_test()
{
  //  std::cout << "test memory module C binding"<<std::endl;
  size_t base=memory_save();
  //  std::cout << "memory_used="<<memory_used(0)<<std::endl;
  double* dd=(double*)memory_allocate(55*sizeof(double));
  //  std::cout << "Allocated pointer at "<< (size_t) dd<<std::endl;
  //  std::cout << "memory_used="<<memory_used(0)<<std::endl;
  for (int i=0; i<55; i++) dd[i]=1;
  // memory_release(dd);
  //  std::cout << "memory_used="<<memory_used(0)<<std::endl;
  memory_release_saved(base);
  //  std::cout << "memory_used="<<memory_used(0)<<std::endl;

  {
    std::vector<double,memory::Allocator<double> > xx(55);
    for (int i=0; i<55; i++) xx[i]=1;
    //    std::cout << xx[54]<<std::endl;
  //    std::cout << "memory_used="<<memory_used(0)<<std::endl;
  }
  //  std::cout << "memory_used="<<memory_used(0)<<std::endl;
}
#endif
#endif
