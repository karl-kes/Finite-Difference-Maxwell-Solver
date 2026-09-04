#pragma once

#include <concepts>
#include <type_traits>

namespace kestrel {

template <typename T>
concept arithmetic = 
  std::floating_point<T> ||
  std::integral<T>;

template <typename T>
concept enumeration = 
  std::is_enum_v<T>;

} // namespace kestrel
