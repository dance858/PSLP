#ifndef PSLP_THREAD_H
#define PSLP_THREAD_H

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
typedef HANDLE ps_thread_t;
typedef DWORD(WINAPI *ps_thread_func_t)(LPVOID);
static inline int ps_thread_create(ps_thread_t *thread, void *attr,
                                   void *(*start_routine)(void *), void *arg)
{
    (void) attr;
    *thread =
        CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE) start_routine, arg, 0, NULL);
    return *thread ? 0 : -1;
}
static inline int ps_thread_join(ps_thread_t thread, void **retval)
{
    DWORD res = WaitForSingleObject(thread, INFINITE);
    if (res == WAIT_OBJECT_0)
    {
        if (retval)
        {
            DWORD code;
            GetExitCodeThread(thread, &code);
            *retval = (void *) (uintptr_t) code;
        }
        CloseHandle(thread);
        return 0;
    }
    return -1;
}
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
