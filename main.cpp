#include <iostream>
#include <memory>
#include <cmath>
#include <fstream>
#include <string>
#include <omp.h>
#include <cassert>
#include <cctype>
#include <chrono>

class Grid {
private:
    std::size_t Nx_, Ny_, Nz_;                  // Grid Size
    double dx_, dy_, dz_;                       // Spatial Differentials
    double c_;                                  // Light Speed
    double dt_;                                 // Time Differential
    std::unique_ptr<double[]> Ex_, Ey_, Ez_;    // Electric Field
    std::unique_ptr<double[]> Bx_, By_, Bz_;    // Magnetic Field

public:
    // Constructor:
    Grid( std::size_t new_Nx = 5, std::size_t new_Ny = 5, std::size_t new_Nz = 5,
          double new_dx = 1.0, double new_dy = 1.0, double new_dz = 1.0,
          double new_c = 1.0 ):
    Nx_{ new_Nx }, Ny_{ new_Ny }, Nz_{ new_Nz },
    dx_{ new_dx }, dy_{ new_dy }, dz_{ new_dz },
    c_{ new_c },
    dt_{ 0.99 / ( c() * std::sqrt( 1.0/(dx()*dx()) + 1.0/(dy()*dy()) + 1.0/(dz()*dz()) ) ) } {
        
        std::size_t const grid_size{ Nx_ * Ny_ * Nz_ };
        Ex_ = std::make_unique<double[]>( grid_size );
        Ey_ = std::make_unique<double[]>( grid_size );
        Ez_ = std::make_unique<double[]>( grid_size );
        Bx_ = std::make_unique<double[]>( grid_size );
        By_ = std::make_unique<double[]>( grid_size );
        Bz_ = std::make_unique<double[]>( grid_size );

        for ( std::size_t i{}; i < grid_size; ++i ) {
            Ex_[i] = 0.0;
            Ey_[i] = 0.0;
            Ez_[i] = 0.0;
            Bx_[i] = 0.0;
            By_[i] = 0.0;
            Bz_[i] = 0.0;
        }
    }

    // System Simulation:
    void update_B() {
        // dB/dt = -curl( E )
        #pragma omp parallel for collapse(2)
        for ( std::size_t z = 0; z < Nz() - 1; ++z ) {
            for ( std::size_t y = 0; y < Ny() - 1; ++y ) {
                #pragma omp simd
                for ( std::size_t x = 0; x < Nx() - 1; ++x ) {
                    // ∂Bx/∂t = -(∂Ez/∂y - ∂Ey/∂z)
                    Bx_[idx(x,y,z)] -= dt() * curl_X( Ey_[idx(x,y,z)], Ey_[idx(x,y,z+1)],
                                                      Ez_[idx(x,y,z)], Ez_[idx(x,y+1,z)] );

                    // ∂By/∂t = -(∂Ex/∂z - ∂Ez/∂x)
                    By_[idx(x,y,z)] -= dt() * curl_Y( Ex_[idx(x,y,z)], Ex_[idx(x,y,z+1)],
                                                      Ez_[idx(x,y,z)], Ez_[idx(x+1,y,z)] );

                    // ∂Bz/∂t = -(∂Ey/∂x - ∂Ex/∂y)
                    Bz_[idx(x,y,z)] -= dt() * curl_Z( Ey_[idx(x,y,z)], Ey_[idx(x+1,y,z)],
                                                      Ex_[idx(x,y,z)], Ex_[idx(x,y+1,z)] );
                }
            }
        }
    }
    void update_E() {
        // dE/dt = c*c * curl(B)
        #pragma omp parallel for collapse(2)
        for ( std::size_t z = 1; z < Nz(); ++z ) {
            for ( std::size_t y = 1; y < Ny(); ++y ) {
                #pragma omp simd
                for ( std::size_t x = 1; x < Nx(); ++x ) {
                    // ∂Ex/∂t = c*c * (∂Ez/∂y - ∂Ey/∂z)
                    Ex_[idx(x,y,z)] += dt() * c_sq() * curl_X( By_[idx(x,y,z-1)], By_[idx(x,y,z)],
                                                               Bz_[idx(x,y-1,z)], Bz_[idx(x,y,z)] );

                    // ∂Ey/∂t = c*c * (∂Ex/∂z - ∂Ez/∂x)
                    Ey_[idx(x,y,z)] += dt() * c_sq() * curl_Y( Bx_[idx(x,y,z-1)], Bx_[idx(x,y,z)],
                                                               Bz_[idx(x-1,y,z)], Bz_[idx(x,y,z)] );

                    // ∂Ez/∂t = c*c * (∂Ex/∂y - ∂Ey/∂x)
                    Ez_[idx(x,y,z)] += dt() * c_sq() * curl_Z( By_[idx(x-1,y,z)], By_[idx(x,y,z)],
                                                               Bx_[idx(x,y-1,z)], Bx_[idx(x,y,z)] );
                }
            }
        }
    }
    void step() {
        update_B();
        update_E();
    }
    void inject_source( std::size_t const x,
                        std::size_t const y,
                        std::size_t const z, 
                        double const value ) {
        Ez_[idx(x,y,z)] += value;
    }
    void output_slice( std::size_t const z, std::string const &file_name ) {
        std::ofstream file( file_name );
        for ( std::size_t y{}; y < Ny(); ++ y ) {
            for ( std::size_t x{}; x< Nx(); ++x ) {
                file << Ez_[idx(x,y,z)];
                if ( x < Nx() - 1 ) {
                    file << ",";
                }
            }
            file << "\n";
        }
    }

    // Getters:
    // Differentials
    double dx() const {
        return dx_;
    }
    double dy() const {
        return dy_;
    }
    double dz() const {
        return dz_;
    }
    double dt() const {
        return dt_;
    }
    // Speed
    double c() const {
        return c_;
    }
    double c_sq() const {
        return c_*c_;
    }
    // Dimensions
    std::size_t Nx() const {
        return Nx_;
    }
    std::size_t Ny() const {
        return Ny_;
    }
    std::size_t Nz() const {
        return Nz_;
    }
    // Fields
    double get_field( char const field,
                      char const component,
                      std::size_t const x,
                      std::size_t const y,
                      std::size_t const z ) const {
        std::size_t index{ idx(x,y,z) };

        if ( std::tolower( field ) == 'e' ) {
            switch ( std::tolower( component ) ) {
                case 'x': return Ex_[index];
                case 'y': return Ey_[index];
                case 'z': return Ez_[index];
            }
        } else if ( std::tolower( field ) == 'b' ) {
            switch ( std::tolower( component ) ) {
                case 'x': return Bx_[index];
                case 'y': return By_[index];
                case 'z': return Bz_[index];
            }
        }
        std::cout << "ERROR! CHECK PARAMETERS!" << std::endl;
        return 0.0;
    }

    // Helpers:
    // Finds 3D index
    std::size_t idx( std::size_t const x, std::size_t const y, std::size_t const z ) const {
        return x + Nx() * ( y + Ny() * z );
    }
    // Curls in X, Y, Z
    double curl_X( double const Y_0, double const Y_1,
                   double const Z_0, double const Z_1 ) const {
        double dZdY{ ( Z_1 - Z_0 ) / dy() };
        double dYdZ{ ( Y_1 - Y_0 ) / dz() };

        return ( dZdY - dYdZ );
    }
    double curl_Y( double const X_0, double const X_1,
                   double const Z_0, double const Z_1 ) const {
        double dXdZ{ ( X_1 - X_0 ) / dz() };
        double dZdX{ ( Z_1 - Z_0 ) / dx() };

        return ( dXdZ - dZdX );
    }
    double curl_Z( double const Y_0, double const Y_1,
                   double const X_0, double const X_1 ) const {
        double dYdX{ ( Y_1 - Y_0 ) / dx() };
        double dXdY{ ( X_1 - X_0 ) / dy() };

        return ( dYdX - dXdY );
    }
};

int main() {
    
    /* 
        To compile and run.

        For No Parallel (NOT RECOMMENDED):
        g++ -std=c++17 main.cpp -o main.exe

        For Parallel (RECOMMENDED):
        g++ -std=c++17 main.cpp -o main.exe -fopenmp

        ./main.exe
    */

    Grid grid{ 100, 100, 5 };

    auto start = std::chrono::high_resolution_clock::now();
    for ( int t{}; t <= 1000; ++t ) {
        grid.inject_source( 5, 5, 1, std::sin( 0.1 * t ) );
        grid.step();
    }
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time: " << duration.count() << " ms" << std::endl;

    return 0;
}