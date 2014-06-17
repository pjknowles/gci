#ifndef DISTANCE_H
#define DISTANCE_H
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
namespace std{//extend std namespace for Sun which does not have regular std::distance
template <class Iterator>
  inline typename iterator_traits<Iterator>::difference_type distance( Iterator i1, Iterator i2 )
{
 typename iterator_traits<Iterator>::difference_type diff=0;
 distance(i1, i2, diff);
 return diff;
}
}
#endif
#endif // DISTANCE_H
