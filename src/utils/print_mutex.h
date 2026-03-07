#ifndef PRINT_MUTEX_H
#define PRINT_MUTEX_H

#include <tbb/mutex.h>

/**
 * Return a process-global mutex for thread-safe console output.
 *
 * Include this header and lock the returned mutex before printing
 * from parallel code:
 * @code
 * {
 *   tbb::mutex::scoped_lock lock(get_print_mutex());
 *   std::cout << "message" << std::endl;
 * }
 * @endcode
 */
inline tbb::mutex& get_print_mutex() {
  static tbb::mutex m;
  return m;
}

#endif // PRINT_MUTEX_H
// }