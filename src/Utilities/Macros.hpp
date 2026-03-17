#pragma once

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

#if defined(__GNUC__) || defined(__clang__)
    #define ASSUME_ALIGNED(ptr, align) \
        (ptr) = static_cast<decltype(ptr)>(__builtin_assume_aligned((ptr), (align)))
#elif defined(_MSC_VER)
    #define ASSUME_ALIGNED(ptr, align) \
        __assume((reinterpret_cast<uintptr_t>(ptr) % (align)) == 0)
#else
    #define ASSUME_ALIGNED(ptr, align) ((void)0)
#endif