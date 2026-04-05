#pragma once

#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>

namespace test {

struct TestCase {
    std::string name;
    std::string suite;
    std::function<void()> fn;
};

struct TestResult {
    std::string name;
    bool passed;
    std::string message;
    double elapsed_ms;
};

inline std::vector<TestCase>& registry() {
    static std::vector<TestCase> cases;
    return cases;
}

struct Registrar {
    Registrar( std::string suite, std::string name, std::function<void()> fn ) {
        registry().push_back( { std::move(name), std::move(suite), std::move(fn) } );
    }
};

// Exception thrown on assertion failure (caught by runner):
struct AssertionFailed : std::exception {
    std::string msg;
    AssertionFailed( std::string m ) : msg{ std::move(m) } {}
    char const* what() const noexcept override { return msg.c_str(); }
};

inline void assert_true( bool cond, char const* expr, char const* file, int line ) {
    if ( !cond ) {
        std::ostringstream ss;
        ss << file << ":" << line << ": ASSERT_TRUE(" << expr << ") failed";
        throw AssertionFailed{ ss.str() };
    }
}

inline void assert_near( double actual, double expected, double tol,
                          char const* a_str, char const* e_str, char const* file, int line ) {
    if ( std::abs( actual - expected ) > tol ) {
        std::ostringstream ss;
        ss << file << ":" << line << ": ASSERT_NEAR(" << a_str << ", " << e_str
           << ", " << tol << ") failed: actual=" << std::setprecision(12) << actual
           << " expected=" << expected << " diff=" << std::abs(actual - expected);
        throw AssertionFailed{ ss.str() };
    }
}

template <typename A, typename B>
void assert_eq( A const& a, B const& b, char const* a_str, char const* b_str, char const* file, int line ) {
    if ( !( a == b ) ) {
        std::ostringstream ss;
        ss << file << ":" << line << ": ASSERT_EQ(" << a_str << ", " << b_str
           << ") failed: " << a << " != " << b;
        throw AssertionFailed{ ss.str() };
    }
}

template <typename A, typename B>
void assert_lt( A const& a, B const& b, char const* a_str, char const* b_str, char const* file, int line ) {
    if ( !( a < b ) ) {
        std::ostringstream ss;
        ss << file << ":" << line << ": ASSERT_LT(" << a_str << ", " << b_str
           << ") failed: " << a << " >= " << b;
        throw AssertionFailed{ ss.str() };
    }
}

template <typename A, typename B>
void assert_gt( A const& a, B const& b, char const* a_str, char const* b_str, char const* file, int line ) {
    if ( !( a > b ) ) {
        std::ostringstream ss;
        ss << file << ":" << line << ": ASSERT_GT(" << a_str << ", " << b_str
           << ") failed: " << a << " <= " << b;
        throw AssertionFailed{ ss.str() };
    }
}

inline int run_all( std::string const& filter = "" ) {
    std::vector<TestResult> results;
    std::size_t passed{}, failed{};

    std::cout << "\n<--- FDTD Test Suite --->\n\n";

    std::string current_suite;

    for ( auto const& tc : registry() ) {
        if ( !filter.empty() && tc.suite.find(filter) == std::string::npos
             && tc.name.find(filter) == std::string::npos ) continue;

        if ( tc.suite != current_suite ) {
            if ( !current_suite.empty() ) std::cout << "\n";
            std::cout << "--- " << tc.suite << " ---\n";
            current_suite = tc.suite;
        }

        auto t0 = std::chrono::high_resolution_clock::now();
        try {
            tc.fn();
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            std::cout << "  PASS  " << tc.name << " (" << std::fixed << std::setprecision(1) << ms << " ms)\n";
            results.push_back( { tc.suite + "/" + tc.name, true, "", ms } );
            ++passed;
        } catch ( AssertionFailed const& e ) {
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            std::cout << "  FAIL  " << tc.name << "\n        " << e.msg << "\n";
            results.push_back( { tc.suite + "/" + tc.name, false, e.msg, ms } );
            ++failed;
        } catch ( std::exception const& e ) {
            auto t1 = std::chrono::high_resolution_clock::now();
            double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
            std::cout << "  FAIL  " << tc.name << " (exception: " << e.what() << ")\n";
            results.push_back( { tc.suite + "/" + tc.name, false, e.what(), ms } );
            ++failed;
        }
    }

    std::cout << "\n========================================\n";
    std::cout << "  Results: " << passed << " passed, " << failed << " failed, "
              << (passed + failed) << " total\n";
    std::cout << "========================================\n\n";

    return failed > 0 ? 1 : 0;
}

} // namespace test

#define TEST(suite, name)                                                      \
    void test_##suite##_##name();                                              \
    static test::Registrar reg_##suite##_##name{ #suite, #name, test_##suite##_##name }; \
    void test_##suite##_##name()

#define ASSERT_TRUE(expr)       test::assert_true((expr), #expr, __FILE__, __LINE__)
#define ASSERT_NEAR(a, e, tol)  test::assert_near((a), (e), (tol), #a, #e, __FILE__, __LINE__)
#define ASSERT_EQ(a, b)         test::assert_eq((a), (b), #a, #b, __FILE__, __LINE__)
#define ASSERT_LT(a, b)         test::assert_lt((a), (b), #a, #b, __FILE__, __LINE__)
#define ASSERT_GT(a, b)         test::assert_gt((a), (b), #a, #b, __FILE__, __LINE__)

#define ASSERT_THROW(expr, exc_type) do { bool caught = false; try { expr; } catch ( exc_type const& ) { caught = true; } if ( !caught ) { std::ostringstream ss; ss << __FILE__ << ":" << __LINE__ << ": ASSERT_THROW(" #expr ", " #exc_type ") - no exception thrown"; throw test::AssertionFailed{ ss.str() }; } } while(0)