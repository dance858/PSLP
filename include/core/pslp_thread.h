#ifndef PSLP_THREAD_H
#define PSLP_THREAD_H

#if defined(_WIN32) || defined(_WIN64)

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct
{
    void *(*start_routine)(void *);
    void *arg;
    void *ret;
} ps_thread_wrapper_t;

typedef struct
{
    // dummy, no real handle needed
    ps_thread_wrapper_t wrapper;
} ps_thread_t;

// run the function immediately (no real threading)
static inline int ps_thread_create(ps_thread_t *t, void *attr,
                                   void *(*start_routine)(void *), void *arg)
{
    (void) attr;
    if (!t) return -1;

    t->wrapper.start_routine = start_routine;
    t->wrapper.arg = arg;
    t->wrapper.ret = start_routine(arg); // call immediately
    return 0;
}

#else /* POSIX */

#include <pthread.h>

typedef pthread_t ps_thread_t;

static inline int ps_thread_create(ps_thread_t *thread, void *attr,
                                   void *(*start_routine)(void *), void *arg)
{
    return pthread_create(thread, (pthread_attr_t *) attr, start_routine, arg);
}

static inline int ps_thread_join(ps_thread_t *thread, void **retval)
{
    return pthread_join(*thread, retval);
}

#endif /* Windows / POSIX */

#endif /* PSLP_THREAD_H */
