#pragma once

namespace kestrel {

#if defined(KESTREL_FP32) && defined(KESTREL_FP64)
  #error "ERROR: Select only one floating-point precision"
#elif defined(KESTREL_FP32)
  using fp_t = float;
#else
  using fp_t = double;
#endif

}
