#include "grid.hpp"

#if defined(__GNUC__) || defined(__clang__)
    #define RESTRICT __restrict__
#elif defined(_MSC_VER)
    #define RESTRICT __restrict
#else
    #define RESTRICT
#endif

Grid::Grid( Simulation_Config const &config )
: Nx_{ config.Nx + 1 }
, Ny_{ config.Ny + 1 }
, Nz_{ config.Nz + 1 }
, dx_{ config.dx }
, dy_{ config.dy }
, dz_{ config.dz }
, eps_{ config.eps }
, mu_{ config.mu }
, c_{ config.c }
, c_sq_{ config.c * config.c }
, dt_{ config.dt }
, pml_{ config } {
    auto allocate = [size = config.size]() {
        return std::make_unique<double[]>( size );
    };

    Ex_ = allocate(); Ey_ = allocate(); Ez_ = allocate();
    Bx_ = allocate(); By_ = allocate(); Bz_ = allocate();
    Jx_ = allocate(); Jy_ = allocate(); Jz_ = allocate();
}

void Grid::add_source( std::unique_ptr<Source> source ) {
    sources_.push_back( std::move( source ) );
}

void Grid::apply_sources( double const time_step ) {
    for ( auto const &source : sources_ ) {
        source->apply( *this, time_step );
    }
}

void Grid::update_B() {
    // ∂B/∂t = -curl( E )

    double* RESTRICT Bx{ Bx_.get() };
    double* RESTRICT By{ By_.get() };
    double* RESTRICT Bz{ Bz_.get() };
    double* RESTRICT Ex{ Ex_.get() };
    double* RESTRICT Ey{ Ey_.get() };
    double* RESTRICT Ez{ Ez_.get() };

    double const dt_local{ dt() };

    #pragma omp parallel for collapse( 2 )
    // Start at 0; stagger for Yee-Cell grid.
    for ( std::size_t z = 0; z < Nz() - 1; ++z ) {
        for ( std::size_t y = 0; y < Ny() - 1; ++y ) {
            #pragma omp simd
            for ( std::size_t x = 0; x < Nx() - 1; ++x ) {
                std::size_t const i{ idx(x,y,z) };

                // Take curl of components and apply B -= ∂B:
                double const curl_x_E{ curl_x( Ey[i], Ey[idx(x,y,z+1)], Ez[i], Ez[idx(x,y+1,z)] ) };
                double const curl_y_E{ curl_y( Ex[i], Ex[idx(x,y,z+1)], Ez[i], Ez[idx(x+1,y,z)] ) };
                double const curl_z_E{ curl_z( Ey[i], Ey[idx(x+1,y,z)], Ex[i], Ex[idx(x,y+1,z)] ) };

                // ∂B_x = ∂t * curl_x( E )
                Bx[i] -= dt_local * curl_x_E;

                // ∂B_y = ∂t * curl_y( E )
                By[i] -= dt_local * curl_y_E;

                // ∂B_z = ∂t * curl_z( E )
                Bz[i] -= dt_local * curl_z_E;
            }
        }
    }
    pml_.update_B_psi( Ex_.get(), Ey_.get(), Ez_.get(),
                       Bx_.get(), By_.get(), Bz_.get(),
                       dt(), dx(), dy(), dz() );
}

void Grid::update_E() {
    // ∂E/∂t = c*c * curl(B)

    double* RESTRICT Bx{ Bx_.get() };
    double* RESTRICT By{ By_.get() };
    double* RESTRICT Bz{ Bz_.get() };
    double* RESTRICT Ex{ Ex_.get() };
    double* RESTRICT Ey{ Ey_.get() };
    double* RESTRICT Ez{ Ez_.get() };
    double* RESTRICT Jx{ Jx_.get() };
    double* RESTRICT Jy{ Jy_.get() };
    double* RESTRICT Jz{ Jz_.get() };

    double const dt_local{ dt() };
    double const inv_eps{ 1.0 / eps() };
    double const c_sq_local{ c_sq() };

    #pragma omp parallel for collapse( 2 )
    // Start at 1; stagger for Yee-Cell grid.
    for ( std::size_t z = 1; z < Nz(); ++z ) {
        for ( std::size_t y = 1; y < Ny(); ++y ) {
            #pragma omp simd
            for ( std::size_t x = 1; x < Nx(); ++x ) {
                std::size_t const i{ idx(x,y,z) };

                // Curl of components and apply E += ∂E:
                double const curl_x_B{ curl_x( By[idx(x,y,z-1)], By[i], Bz[idx(x,y-1,z)], Bz[i] ) };
                double const curl_y_B{ curl_y( Bx[idx(x,y,z-1)], Bx[i], Bz[idx(x-1,y,z)], Bz[i] ) };
                double const curl_z_B{ curl_z( By[idx(x-1,y,z)], By[i], Bx[idx(x,y-1,z)], Bx[i] ) };

                double const jx_term{ Jx[i] * inv_eps };
                double const jy_term{ Jy[i] * inv_eps };
                double const jz_term{ Jz[i] * inv_eps };

                // ∂E_x = ∂t * c*c * (∂E_z/∂y - ∂E_y/∂z)
                Ex[i] += dt_local * ( c_sq_local * curl_x_B - jx_term );

                // ∂E_y = ∂t * c*c * (∂Ex/∂z - ∂Ez/∂x)
                Ey[i] += dt_local * ( c_sq_local * curl_y_B - jy_term );

                // ∂E_z = ∂t * c*c * (∂Ex/∂y - ∂Ey/∂x)
                Ez[i] += dt_local * ( c_sq_local * curl_z_B - jz_term );
            }
        }
    }
    pml_.update_E_psi( Ex_.get(), Ey_.get(), Ez_.get(),
                       Bx_.get(), By_.get(), Bz_.get(),
                       dt(), dx(), dy(), dz(), c_sq() );
}

void Grid::step() {
    update_B();
    update_E();
}

double Grid::field(
    Field const field,
    Component const component,
    std::size_t const x, std::size_t const y, std::size_t const z ) const {
    if ( field == Field::ELECTRIC ) {
        switch ( component ) {
            case Component::X: return Ex(x,y,z);
            case Component::Y: return Ey(x,y,z);
            case Component::Z: return Ez(x,y,z);
        }
    } else if ( field == Field::MAGNETIC ) {
        switch ( component ) {
            case Component::X: return Bx(x,y,z);
            case Component::Y: return By(x,y,z);
            case Component::Z: return Bz(x,y,z);
        }
    }
    throw std::invalid_argument{ "Invalid field or component specifier" };
}

double &Grid::field(
    Field const field,
    Component const component,
    std::size_t const x, std::size_t const y, std::size_t const z ) {
    if ( field == Field::ELECTRIC ) {
        switch ( component ) {
            case Component::X: return Ex(x,y,z);
            case Component::Y: return Ey(x,y,z);
            case Component::Z: return Ez(x,y,z);
        }
    } else if ( field == Field::MAGNETIC ) {
        switch ( component ) {
            case Component::X: return Bx(x,y,z);
            case Component::Y: return By(x,y,z);
            case Component::Z: return Bz(x,y,z);
        }
    }
    throw std::invalid_argument{ "Invalid field or component specifier" };
}

double Grid::field_magnitude(
    Field const field,
    std::size_t const x, std::size_t const y, std::size_t const z ) const {
    double Fx{ this->field( field, Component::X, x, y, z ) };
    double Fy{ this->field( field, Component::Y, x, y, z ) };
    double Fz{ this->field( field, Component::Z, x, y, z ) };

    return std::sqrt( Fx*Fx + Fy*Fy + Fz*Fz );
}

double Grid::total_energy() const {
    double energy{};
    double const dV{ dx() * dy() * dz() };

    #pragma omp parallel for collapse( 3 ) reduction( +:energy )
    for ( std::size_t z = 0; z < Nz(); ++z ) {
        for ( std::size_t y = 0; y < Ny(); ++y ) {
            for ( std::size_t x = 0; x < Nx(); ++x ) {
                    double const E_sq{ Ex(x,y,z)*Ex(x,y,z) + Ey(x,y,z)*Ey(x,y,z) + Ez(x,y,z)*Ez(x,y,z) };
                    double const B_sq{ Bx(x,y,z)*Bx(x,y,z) + By(x,y,z)*By(x,y,z) + Bz(x,y,z)*Bz(x,y,z) };

                    energy += 0.5 * ( eps() * E_sq + B_sq / mu() );
            }
        }
    }
    return energy * dV;
}

double Grid::source_power() const {
    double power{};
    double const dV{ dx() * dy() * dz() };
    
    #pragma omp parallel for collapse( 3 ) reduction( +:power )
    for ( std::size_t z = 0; z < Nz(); ++z ) {
        for ( std::size_t y = 0; y < Ny(); ++y ) {
            for ( std::size_t x = 0; x < Nx(); ++x ) {
                power -= Jx(x,y,z) * Ex(x,y,z);
                power -= Jy(x,y,z) * Ey(x,y,z);
                power -= Jz(x,y,z) * Ez(x,y,z);
            }
        }
    }
    return power * dV;
}