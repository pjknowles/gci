#ifndef MEMORY_H
#define MEMORY_H
#include <stddef.h>
#include <stdint.h>
#include <iostream>
#include <vector>
#include <stdexcept>
#ifndef NDEBUG
#define NDEBUG
#endif
#include <assert.h>
#include <iostream>
#ifdef __cplusplus
extern "C" {
#endif
/*!
 * \brief memory_initialize Initialize the memory system
 * \param n The number of bytes to be reserved for the stack
 */
void memory_initialize(size_t n);
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
  template<typename T> class allocator {
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
        typedef allocator<U> other;
      };

  public :
    inline explicit allocator() {}
    inline ~allocator() {}
    inline explicit allocator(allocator const&) {}
    template<typename U>
      inline explicit allocator(allocator<U> const&) {}

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

    inline bool operator==(allocator const&) { return true; }
    inline bool operator!=(allocator const& a) { return !operator==(a); }
  };
}

namespace memory {


/*!
 * @brief A template for a container class like std::vector<T> but with the
 * memory provided from an existing array or std::vector, or allocate internally.
 * This allows one to optionally attach the
 * container to existing data structures, and use the same syntax as std::vector<T>
 * to access it.  Note that the array is fixed size, defined at construction, and so those methods of
 * std::vector that allocate memory
 * (for example resize()) can not be used.
 * \tparam T Type of the elements of the array.
 * \tparam Alloc Memory allocator. By default, this is the allocator in the memory namespace, which allocates from a managed fixed stack.
 * Note that with that allocator, buffers must be released in reverse order of construction, but this will normally be the case when relying
 * on the destructor of the vector class being called when objects go out of scope.
 *
 * Example:
 * \code
 * using namespace memory;
 * {
 *   vector<double> arr(5);
 *   for (size_t i=0; i<arr.size(); i++) el[i]=100*i;
 *   for (vector<double>::iterator el=arr.begin(); el!=arr.end(); el++)
 *     std::cout << "array element "<<*el<<std::endl;
 * } // arr will release its buffer and be destroyed on going out of scope
 * \endcode
 */
template<typename T, class Alloc=memory::allocator<T> > class vector {
private:
  vector<T,Alloc>() {}
  T* buffer;
  size_t length;
  std::vector<T,Alloc> * vp; // if not null, then this object owns the buffer

  template<class T1>
    class vector_iterator : std::iterator<std::bidirectional_iterator_tag, T1>
  {
    public:
      vector_iterator(T1* ptr): ptr_(ptr) { }
      vector_iterator operator++() { ptr_++; return *this;}
      vector_iterator operator++(int junk) { vector_iterator i=*this; ptr_++; return i;}
      vector_iterator operator=(const vector_iterator& other) { ptr_ = other.ptr_; return *this; }
      ptrdiff_t operator+=(const size_t incr) {ptr_+=incr;return *this;}
      ptrdiff_t operator-=(const size_t incr) {ptr_-=incr;return *this;}
      ptrdiff_t operator-(const vector_iterator& other) const {return this->ptr_-other.ptr_;}
      T1& operator*() { return *ptr_;}
      T1* operator->() { return ptr_;}
      template<class T2> friend bool operator==(vector_iterator const &lhs, vector_iterator<T2> const &rhs) { return lhs.ptr_==rhs.ptr_;}
      template<class T2> friend bool operator!=(vector_iterator const &lhs, vector_iterator<T2> const &rhs) { return lhs.ptr_!=rhs.ptr_;}
      template<class T2> friend bool operator<(vector_iterator const &lhs, vector_iterator<T2> const &rhs) { return lhs.ptr_<rhs.ptr_;}
      template<class T2> friend bool operator>(vector_iterator const &lhs, vector_iterator<T2> const &rhs) { return lhs.ptr_>rhs.ptr_;}
    private:
      T1* ptr_;
      vector_iterator(); // private for begin, end
    };

public:
  /*!
   * \brief Construct a vector of type T allocated dynamically.
   * \param length The number of elements of buffer.
   */
  inline vector<T,Alloc>(size_t const length) : length(length), vp(new std::vector<T,Alloc>(length))
  {
    this->buffer=&((*vp)[0]);
  }

  /*!
   * \brief Construct a vector of type T.
   * \param buffer Pointer to an existing array of type T.
   * \param length The number of elements of buffer.
   */
  inline vector<T,Alloc>(T * const buffer, size_t const length) : vp(0), buffer(buffer), length(length) { }
  /*!
   * \brief Construct a vector of type T, and initialise its elements.
   * \param buffer Pointer to an existing array of type T.
   * \param length The number of elements of buffer.
   * \param value The value to assign to all elements of the buffer.
   */
  inline vector<T,Alloc>(T * const buffer, size_t const length, T value) : vp(0), buffer(buffer), length(length)
  {
    for (typename vector<T,Alloc>::iterator v=this->begin(); v!=this->end(); v++) *v=value;
  }

  /*!
   * \brief Construct a vector of type T.
   * \param array Existing std::vector of type T whose buffer will be used. Any subsequent resizing of array will cause chaos.
   */
  inline vector<T,Alloc>(std::vector<T,Alloc>& array) : vp(0), buffer(&array[0]), length(array.size()) { }
  /*!
   * \brief Construct a vector of type T, and initialise its elements.
   * \param array Existing std::vector of type T whose buffer will be used. Any subsequent resizing of array will cause chaos.
   * \param value The value to assign to all elements of the buffer.
   */
  inline vector<T,Alloc>(std::vector<T,Alloc>& array, T value) : vp(0), buffer(&array[0]), length(array.size()) { assign(value); }

  /*!
   * \brief vector<T> Copy constructor
   * \param source An existing object. An element-by-element copy is made, i.e. the data buffer is allocated then copied from source.
   */
  inline vector<T,Alloc>(const vector<T,Alloc>& source) : length(source.size()), vp(new std::vector<T,Alloc>(source.size()))
  {
    this->buffer=&((*vp)[0]);
    for (typename source::const_iterator v=source.begin(); v!=source.end(); v++)
      buffer[v-source.begin()]=*v;
  }

  inline ~vector() { if (vp) delete vp; }

  T& operator[](size_t n) { assert(n<length && n>=0); return buffer[n]; }
  /*!
   * \brief Assign new contents to the vector, replacing its current contents
   * \param value The value to set all elements of the vector to.
   */
  void assign(const T& value) {
    for (typename vector<T,Alloc>::iterator v=this->begin(); v!=this->end(); v++) *v=value;
  }
  template <class InputIterator>
    void assign (InputIterator first, InputIterator last) {
    for (typename vector<T,Alloc>::iterator v=first; v!=last; v++) buffer[v-first]=*v;
  }
  inline T& at(size_t n) { if (n<length && n>=0) return buffer[n]; else throw std::out_of_range("vector:: at() out of range"); }
  inline const T& at(size_t n) const { if (n<length && n>=0) return buffer[n]; else throw std::out_of_range("vector:: at() out of range"); }
  inline T& back() { return buffer[size()-1]; }
  inline const T& back() const { return buffer[size()-1]; }
  inline size_t capacity() const { return length; }
  inline bool empty() const { return length==0; }
  inline T& front() { return buffer[0]; }
  inline const T& front() const { return buffer[0]; }
  inline size_t max_size() const { return length; }
  inline size_t size() const { return length; }
  /*!
   * \brief Exchange the content of the container by the content of x, which is another object of the same type.
   * \param x Another vector of the same type.
   */
  void swap(vector<T,Alloc>& x) {
    std::swap(this->buffer,x.buffer);
    std::swap(this->length,x.length);
    std::swap(this->vp,x.vp);
  }

  typedef vector_iterator<T> iterator;
  typedef vector_iterator<T const> const_iterator;

  iterator begin() {return iterator(buffer);}
  iterator end() {return iterator(buffer+length);}
  const_iterator begin() const {return const_iterator(buffer);}
  const_iterator end() const {return const_iterator(buffer+length);}
};


}
using namespace memory;

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

  // test memory::allocator with STL container
  {
    std::vector<double,memory::allocator<double> > xx(55);
    for (int i=0; i<55; i++) xx[i]=1;
    //    std::cout << xx[54]<<std::endl;
  //    std::cout << "memory_used="<<memory_used(0)<<std::endl;
  }
  // test memory::vector container
  {
    vector<double> xx(5);
    double val(5);
    xx.assign(val);
    for (vector<double>::iterator x=xx.begin(); x<xx.end(); x++)
      if (*x != val) throw "bad";
  }
  //  std::cout << "memory_used="<<memory_used(0)<<std::endl;
}



#endif
#endif
