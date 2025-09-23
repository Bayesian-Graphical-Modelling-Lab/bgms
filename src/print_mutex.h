#ifndef PRINT_MUTEX_H
#define PRINT_MUTEX_H

#include <tbb/mutex.h>

inline tbb::mutex& get_print_mutex() {
  static tbb::mutex m;
  return m;
}

#endif // PRINT_MUTEX_H
