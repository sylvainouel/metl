#ifndef OMP_STUB_H
#define OMP_STUB_H
//#warning COMPILED WITHOUT OMP SUPPORT
inline int omp_get_num_threads() { return 1; }
inline int omp_get_thread_num() { return 0; }
inline int omp_get_max_threads() { return 1; }
inline void omp_set_num_threads(int t) {}

typedef int omp_lock_t;
#endif
