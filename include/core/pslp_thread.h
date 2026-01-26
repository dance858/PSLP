#ifndef PSLP_THREAD_H
#define PSLP_THREAD_H

#if defined(_WIN32) || defined(_WIN64)

#include <stdlib.h>
#include <windows.h>

typedef HANDLE ps_thread_t;

/* Wrapper to adapt POSIX-style start routine */
typedef struct
{
    void *(*start_routine)(void *);
    void *arg;
    void *ret;
} ps_thread_wrapper_t;

static DWORD WINAPI ps_thread_trampoline(LPVOID param)
{
    ps_thread_wrapper_t *w = (ps_thread_wrapper_t *) param;
    w->ret = w->start_routine(w->arg);
    return 0;
}

static inline int ps_thread_create(ps_thread_t *thread, void *attr,
                                   void *(*start_routine)(void *), void *arg)
{
    (void) attr;

    ps_thread_wrapper_t *w = malloc(sizeof(*w));
    if (!w) return -1;

    w->start_routine = start_routine;
    w->arg = arg;
    w->ret = NULL;

    *thread = CreateThread(NULL, 0, ps_thread_trampoline, w, 0, NULL);

    if (!*thread)
    {
        free(w);
        return -1;
    }

    return 0;
}

static inline int ps_thread_join(ps_thread_t thread, void **retval)
{
    if (WaitForSingleObject(thread, INFINITE) != WAIT_OBJECT_0) return -1;

    /* Recover wrapper pointer */
    ps_thread_wrapper_t *w = NULL;
    GetExitCodeThread(thread, (LPDWORD) &w);

    if (retval && w) *retval = w->ret;

    free(w);
    CloseHandle(thread);
    return 0;
}

/* Best-effort only: NOT safe */
static inline int ps_thread_cancel(ps_thread_t thread)
{
    TerminateThread(thread, 0);
    CloseHandle(thread);
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

#endif
#endif /* PSLP_THREAD_H */
