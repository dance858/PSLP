#ifndef PSLP_THREAD_H
#define PSLP_THREAD_H

#if defined(_WIN32) || defined(_WIN64)
#include <stdint.h>
#include <windows.h>
// Portable thread type
typedef HANDLE ps_thread_t;

// Internal structure to pass user function and argument
typedef struct
{
    void *(*start_routine)(void *);
    void *arg;
    void *ret;
} ps_thread_wrapper_t;

// Windows thread wrapper to safely return pointer
static DWORD WINAPI ps_thread_start_wrapper(LPVOID param)
{
    ps_thread_wrapper_t *wrapper = (ps_thread_wrapper_t *) param;
    wrapper->ret = wrapper->start_routine(wrapper->arg);
    return 0;
}

// Thread create: attr is ignored on Windows
static inline int ps_thread_create(ps_thread_t *thread, void *attr,
                                   void *(*start_routine)(void *), void *arg)
{
    (void) attr; // Not used on Windows
    ps_thread_wrapper_t *wrapper =
        (ps_thread_wrapper_t *) malloc(sizeof(ps_thread_wrapper_t));
    if (!wrapper) return -1;
    wrapper->start_routine = start_routine;
    wrapper->arg = arg;
    wrapper->ret = NULL;
    *thread = CreateThread(NULL, 0, ps_thread_start_wrapper, wrapper, 0, NULL);
    return *thread ? 0 : -1;
}

// Thread join: retrieves pointer return value if possible
static inline int ps_thread_join(ps_thread_t thread, void **retval)
{
    DWORD res = WaitForSingleObject(thread, INFINITE);
    if (res == WAIT_OBJECT_0)
    {
        // Retrieve wrapper and return value
        ps_thread_wrapper_t *wrapper = NULL;
        BOOL ok = GetExitCodeThread(thread, (LPDWORD) &wrapper);
        if (ok && wrapper)
        {
            if (retval) *retval = wrapper->ret;
            free(wrapper);
        }
        else
        {
            if (retval) *retval = NULL;
        }
        CloseHandle(thread);
        return 0;
    }
    return -1;
}

// Thread cancel: WARNING - unsafe, does not clean up resources!
// Use only if absolutely necessary. Prefer cooperative cancellation.
static inline int ps_thread_cancel(ps_thread_t thread)
{
    TerminateThread(thread, 0);
    CloseHandle(thread);
    return 0;
}
#else
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

#endif // PSLP_THREAD_H
