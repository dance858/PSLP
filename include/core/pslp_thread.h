#ifndef PSLP_THREAD_H
#define PSLP_THREAD_H

#if defined(_WIN32) || defined(_WIN64)

/* ======================= */
/*        WINDOWS          */
/* ======================= */

#include <process.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdlib.h>
#include <windows.h>

typedef struct
{
    HANDLE handle;
    atomic_int cancel_requested;
    void *retval;
} ps_thread_t;

typedef struct
{
    void *(*fn)(void *);
    void *arg;
    ps_thread_t *thread;
} ps_thread_start_t;

static unsigned __stdcall ps_thread_trampoline(void *arg)
{
    ps_thread_start_t *s = (ps_thread_start_t *) arg;
    s->thread->retval = s->fn(s->arg);
    free(s);
    return 0;
}

static inline int ps_thread_create(ps_thread_t *thread, void *attr,
                                   void *(*start_routine)(void *), void *arg)
{
    (void) attr;

    atomic_init(&thread->cancel_requested, 0);
    thread->retval = NULL;

    ps_thread_start_t *start = malloc(sizeof(*start));
    if (!start) return -1;

    start->fn = start_routine;
    start->arg = arg;
    start->thread = thread;

    uintptr_t h = _beginthreadex(NULL, 0, ps_thread_trampoline, start, 0, NULL);

    if (!h)
    {
        free(start);
        return -1;
    }

    thread->handle = (HANDLE) h;
    return 0;
}

static inline int ps_thread_join(ps_thread_t *thread, void **retval)
{
    WaitForSingleObject(thread->handle, INFINITE);

    if (retval) *retval = thread->retval;

    CloseHandle(thread->handle);
    return 0;
}

/* Cooperative cancellation */
static inline int ps_thread_cancel(ps_thread_t *thread)
{
    atomic_store(&thread->cancel_requested, 1);
    return 0;
}

static inline int ps_thread_test_cancel(ps_thread_t *thread)
{
    return atomic_load(&thread->cancel_requested);
}

#else

/* ======================= */
/*          POSIX          */
/* ======================= */

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

#endif

#endif /* PSLP_THREAD_H */
