#ifndef GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERBASE_H
#define GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERBASE_H

namespace molpro::gci::schedule {
/*!
 * @brief Base class for shared counters defining the interface.
 */
class SharedCounterBase {
public:
  SharedCounterBase() = default;
  SharedCounterBase(const SharedCounterBase &) = delete;
  virtual ~SharedCounterBase() = default;

  /*!
   * @brief increment the counter
   * @param amount how much to increment by
   * @return current value of the counter
   */
  virtual int increment(int amount = 1) = 0;

  /*!
   * @brief reset counter value to 0
   */
  virtual void reset() = 0;
};

}; // namespace molpro::gci::schedule

#endif // GCI_SRC_MOLPRO_GCI_SCHEDULE_SHAREDCOUNTERBASE_H
