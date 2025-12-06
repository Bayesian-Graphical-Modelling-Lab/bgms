#ifndef PRINT_MUTEX_H
#define PRINT_MUTEX_H

#include <tbb/mutex.h>

inline tbb::mutex& get_print_mutex() {
  static tbb::mutex m;
  return m;
}

#endif // PRINT_MUTEX_H

// Add this header to the parallel code you wish to print from
// + the below code to print in parallel code:
//
// {
//   tbb::mutex::scoped_lock lock(get_print_mutex());
//   std::cout
//   << "print "
//   << std::endl;
// }