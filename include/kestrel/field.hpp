#pragma once

#include <cstddef>
#include "concepts.hpp"

namespace kestrel {

template<arithmetic T, std::size_t components>
class [[nodiscard]] field {
static_assert(components > 0);

private:

public:
};

template <arithmetic T>
using scalar_field = field<T, 1uz>;

template <arithmetic T>
using vector_field = field<T, 3uz>;

} // namespace kestrel
