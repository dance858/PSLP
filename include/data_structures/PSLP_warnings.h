#ifndef PSLP_WARNINGS_H
#define PSLP_WARNINGS_H

#if defined(__clang__)
#define PSLP_DIAG_PUSH() _Pragma("clang diagnostic push")
#define PSLP_DIAG_POP() _Pragma("clang diagnostic pop")

#define PSLP_DIAG_IGNORE_CONVERSION()                                               \
    _Pragma("clang diagnostic ignored \"-Wconversion\"")

#define PSLP_DIAG_IGNORE_SIGN_CONVERSION()                                          \
    _Pragma("clang diagnostic ignored \"-Wsign-conversion\"")

#elif defined(__GNUC__)
#define PSLP_DIAG_PUSH() _Pragma("GCC diagnostic push")
#define PSLP_DIAG_POP() _Pragma("GCC diagnostic pop")

#define PSLP_DIAG_IGNORE_CONVERSION()                                               \
    _Pragma("GCC diagnostic ignored \"-Wconversion\"")

#define PSLP_DIAG_IGNORE_SIGN_CONVERSION()                                          \
    _Pragma("GCC diagnostic ignored \"-Wsign-conversion\"")

#elif defined(_MSC_VER)
#define PSLP_DIAG_PUSH() __pragma(warning(push))
#define PSLP_DIAG_POP() __pragma(warning(pop))

/* Closest MSVC equivalents */
#define PSLP_DIAG_IGNORE_CONVERSION() __pragma(warning(disable : 4244 4267))

#define PSLP_DIAG_IGNORE_SIGN_CONVERSION() __pragma(warning(disable : 4245))

#else
#define PSLP_DIAG_PUSH()
#define PSLP_DIAG_POP()
#define PSLP_DIAG_IGNORE_CONVERSION()
#define PSLP_DIAG_IGNORE_SIGN_CONVERSION()
#endif

#endif /* PSLP_WARNINGS_H */
