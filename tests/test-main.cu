#include "test-framework.hpp"

// Pull in all test translation units (they self-register via static Registrar):
// (Linked at compile time via separate .cu files)

int main( int argc, char* argv[] ) {
    std::string filter;
    if ( argc > 1 ) filter = argv[1];
    return test::run_all( filter );
}