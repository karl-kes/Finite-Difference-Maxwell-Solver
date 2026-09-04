#pragma once

#define RESTRICT __restrict__

#define ASSUME_ALIGNED(ptr, align) \
    (ptr) = static_cast<decltype(ptr)>(__builtin_assume_aligned((ptr), (align)))
