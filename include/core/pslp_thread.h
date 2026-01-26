#ifndef PSLP_THREAD_H
#define PSLP_THREAD_H

#if defined(_WIN32) || defined(_WIN64)

#include <stdint.h>
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

typedef struct
{
    HANDLE handle;
    ps_thread_wrapper_t *wrapper;
} ps_thread_t;

static DWORD WINAPI ps_thread_trampoline(LPVOID param)
{
    ps_thread_wrapper_t *w = param;
    w->ret = w->start_routine(w->arg);
    return 0;
}

static inline int ps_thread_create(ps_thread_t *t, void *attr,
                                   void *(*start_routine)(void *), void *arg)
{
    (void) attr;

    t->wrapper = malloc(sizeof(*t->wrapper));
    if (!t->wrapper) return -1;

    t->wrapper->start_routine = start_routine;
    t->wrapper->arg = arg;
    t->wrapper->ret = NULL;

    t->handle = CreateThread(NULL, 0, ps_thread_trampoline, t->wrapper, 0, NULL);

    if (!t->handle)
    {
        free(t->wrapper);
        return -1;
    }

    return 0;
}

static inline int ps_thread_join(ps_thread_t t, void **retval)
{
    if (WaitForSingleObject(t.handle, INFINITE) != WAIT_OBJECT_0) return -1;

    if (retval) *retval = t.wrapper->ret;

    CloseHandle(t.handle);
    free(t.wrapper);
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

static inline int ps_thread_join(ps_thread_t thread, void **retval)
{
    return pthread_join(thread, retval);
}

static inline int ps_thread_cancel(ps_thread_t thread)
{
    return pthread_cancel(thread);
}

#endif /* Windows / POSIX */

#endif /* PSLP_THREAD_H */
