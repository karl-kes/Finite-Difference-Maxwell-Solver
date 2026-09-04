#pragma once

#include "concepts.hpp"
#include <cstddef>
#include <utility>

namespace kestrel {

template <enumeration Enum>
[[nodiscard]] constexpr auto to_index(Enum value) noexcept -> std::size_t {
  return static_cast<std::size_t>(std::to_underlying(value));
}

enum class axis : std::size_t { x, y, z, count };

inline constexpr auto axis_count{std::to_underlying(axis::count)};

} // namespace kestrel
