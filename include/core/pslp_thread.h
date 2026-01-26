#ifndef PSLP_THREAD_H
#define PSLP_THREAD_H

#if defined(_WIN32) || defined(_WIN64)

#include <process.h> // for _beginthreadex
#include <stdio.h>
#include <stdlib.h>
#include <windows.h>

typedef struct
{
    void *(*start_routine)(void *);
    void *arg;
    void *ret;
} ps_thread_wrapper_t;

typedef struct
{
    HANDLE handle;
    ps_thread_wrapper_t *wrapper;
} ps_thread_t;

// trampoline for _beginthreadex
static unsigned __stdcall ps_thread_trampoline(void *param)
{
    ps_thread_wrapper_t *w = (ps_thread_wrapper_t *) param;
    w->ret = w->start_routine(w->arg);
    return 0;
}

static inline int ps_thread_create(ps_thread_t *t, void *attr,
                                   void *(*start_routine)(void *), void *arg)
{
    (void) attr;

    t->wrapper = (ps_thread_wrapper_t *) calloc(1, sizeof(*t->wrapper));
    if (!t->wrapper) return -1;

    t->wrapper->start_routine = start_routine;
    t->wrapper->arg = arg;
    t->wrapper->ret = NULL;

    uintptr_t handle = _beginthreadex(NULL, // security
                                      0,    // stack size (0 = default)
                                      ps_thread_trampoline,
                                      t->wrapper, // argument
                                      0,          // create flags
                                      NULL        // thread id
    );

    if (handle == 0)
    {
        free(t->wrapper);
        t->wrapper = NULL;
        return -1;
    }

    t->handle = (HANDLE) handle;
    return 0;
}

static inline int ps_thread_join(ps_thread_t *t, void **retval)
{
    if (!t || !t->handle) return -1;

    DWORD wait_result = WaitForSingleObject(t->handle, INFINITE);
    if (wait_result != WAIT_OBJECT_0)
    {
        return -1;
    }

    if (retval && t->wrapper) *retval = t->wrapper->ret;

    CloseHandle(t->handle);
    t->handle = NULL;

    free(t->wrapper);
    t->wrapper = NULL;

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

static inline int ps_thread_cancel(ps_thread_t thread)
{
    return pthread_cancel(thread);
}

#endif /* Windows / POSIX */

#endif /* PSLP_THREAD_H */
