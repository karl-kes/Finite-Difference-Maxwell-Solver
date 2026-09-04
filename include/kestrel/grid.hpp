#pragma once

#include <array>
#include "enums.hpp"
#include "precision.hpp"

namespace kestrel {

using extent = std::array<std::size_t, axis_count>;
using spacing = std::array<fp_t, axis_count>;

class grid {
private:
  extent cells_;
  spacing cell_spacing_;

public:
  explicit grid(extent cells, spacing cell_spacing) noexcept
    : cells_{cells}
    , cell_spacing_{cell_spacing}
  { }

  [[nodiscard]] constexpr auto cells() const noexcept -> extent {
    return cells_;
  }

  [[nodiscard]] constexpr auto cell_spacing() const noexcept -> spacing {
    return cell_spacing_;
  }
};

} // namespace kestrel
