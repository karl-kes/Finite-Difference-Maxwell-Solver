#include "grid.hpp"

void Grid::add_source( std::unique_ptr<Source> source ) {
    sources_.push_back( std::move( source ) );
}

void Grid::apply_sources( double time_step ) {
    for ( auto const &source : sources_ ) {
        source->apply( *this, time_step );
    }
}

void Grid::update_B() {
    // ∂B/∂t = -curl( E )
    #pragma omp parallel for collapse( 2 )
    // Start at 0; stagger for Yee-Cell grid.
    for ( std::size_t z = 0; z < Nz() - 1; ++z ) {
        for ( std::size_t y = 0; y < Ny() - 1; ++y ) {
            #pragma omp simd
            for ( std::size_t x = 0; x < Nx() - 1; ++x ) {
                // Take curl of components and apply B -= ∂B:
                double const curl_x_E{ curl_x( Ey(x,y,z), Ey(x,y,z+1), Ez(x,y,z), Ez(x,y+1,z) ) };
                double const curl_y_E{ curl_y( Ex(x,y,z), Ex(x,y,z+1), Ez(x,y,z), Ez(x+1,y,z) ) };
                double const curl_z_E{ curl_z( Ey(x,y,z), Ey(x+1,y,z), Ex(x,y,z), Ex(x,y+1,z) ) };

                // ∂B_x = ∂t * curl_x( E )
                Bx(x,y,z) -= dt() * curl_x_E;

                // ∂B_y = ∂t * curl_y( E )
                By(x,y,z) -= dt() * curl_y_E;

                // ∂B_z = ∂t * curl_z( E )
                Bz(x,y,z) -= dt() * curl_z_E;
            }
        }
    }
}

void Grid::update_E() {
    // ∂E/∂t = c*c * curl(B)
    #pragma omp parallel for collapse( 2 )
    // Start at 1; stagger for Yee-Cell grid.
    for ( std::size_t z = 1; z < Nz(); ++z ) {
        for ( std::size_t y = 1; y < Ny(); ++y ) {
            #pragma omp simd
            for ( std::size_t x = 1; x < Nx(); ++x ) {
                // Curl of components and apply E += ∂E:
                double const curl_x_B{ curl_x( By(x,y,z-1), By(x,y,z), Bz(x,y-1,z), Bz(x,y,z) ) };
                double const curl_y_B{ curl_y( Bx(x,y,z-1), Bx(x,y,z), Bz(x-1,y,z), Bz(x,y,z) ) };
                double const curl_z_B{ curl_z( By(x-1,y,z), By(x,y,z), Bx(x,y-1,z), Bx(x,y,z) ) };

                double const jx_term{ Jx(x,y,z) / eps() };
                double const jy_term{ Jy(x,y,z) / eps() };
                double const jz_term{ Jz(x,y,z) / eps() };

                // ∂E_x = ∂t * c*c * (∂E_z/∂y - ∂E_y/∂z)
                Ex(x,y,z) += dt() * ( c_sq() * curl_x_B - jx_term );

                // ∂E_y = ∂t * c*c * (∂Ex/∂z - ∂Ez/∂x)
                Ey(x,y,z) += dt() * ( c_sq() * curl_y_B - jy_term );

                // ∂E_z = ∂t * c*c * (∂Ex/∂y - ∂Ey/∂x)
                Ez(x,y,z) += dt() * ( c_sq() * curl_z_B - jz_term );
            }
        }
    }
}

void Grid::step( Simulation_Config const &config, Output const &output, std::size_t curr_time ) {
    update_B();
    update_E();

    if ( ( curr_time % config.output_interval() ) == 0 ) {
        output.write_field( *this, Field::ELECTRIC, curr_time );
        output.write_field( *this, Field::MAGNETIC, curr_time );
        print_progress( curr_time, config.total_time );
    }
}

double Grid::field(
    Field field,
    Component component,
    std::size_t x, std::size_t y, std::size_t z ) const {
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
    Field field,
    Component component,
    std::size_t x, std::size_t y, std::size_t z ) {
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
    Field field,
    std::size_t x, std::size_t y, std::size_t z ) const {
    double Fx{ this->field( field, Component::X, x, y, z ) };
    double Fy{ this->field( field, Component::Y, x, y, z ) };
    double Fz{ this->field( field, Component::Z, x, y, z ) };

    return std::sqrt( Fx*Fx + Fy*Fy + Fz*Fz );
}

double Grid::total_energy() const {
    double energy{};
    double dV{ dx() * dy() * dz() };

    #pragma omp parallel for collapse( 2 ) reduction( +:energy )
    for ( std::size_t z = 0; z < Nz(); ++z ) {
        for ( std::size_t y = 0; y < Ny(); ++y ) {
            #pragma omp simd
            for ( std::size_t x = 0; x < Nx(); ++x ) {
                    double E_sq{ Ex(x,y,z)*Ex(x,y,z) + Ey(x,y,z)*Ey(x,y,z) + Ez(x,y,z)*Ez(x,y,z) };
                    double B_sq{ Bx(x,y,z)*Bx(x,y,z) + By(x,y,z)*By(x,y,z) + Bz(x,y,z)*Bz(x,y,z) };

                    energy += 0.5 * ( eps() * E_sq + B_sq / mu() );
            }
        }
    }
    return energy * dV;
}