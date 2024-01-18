#pragma once

#include <time.h>
#include <iomanip>
#include <sstream>

namespace WHU {
    extern "C" char* strptime(const char* s,
        const char* f,
        struct tm* tm);

    template<typename T, typename U>
        requires std::floating_point<T>&& std::floating_point<U>
    bool are_equal_floating_numbers(T l, U r) {
        std::common_type_t<T, U> l_t = l, r_t = r;
        return std::abs(l_t - r_t) < 1e-6;
    }
}