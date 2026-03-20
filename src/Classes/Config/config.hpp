#pragma once

#include <cstddef>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

enum class Field {
    ELECTRIC,
    MAGNETIC
};

enum class Component {
    X,
    Y,
    Z
};

class Simulation_Config {
public:
    // Independent parameters defaults; overridable in config.cfg:
    std::size_t Nx{ 100 };
    std::size_t Ny{ 100 };
    std::size_t Nz{ 100 };

    double dx{ 1.0 };
    double dy{ 1.0 };
    double dz{ 1.0 };

    double mu{ 1.0 };
    double eps{ 1.0 };

    double cfl_factor{ 1.0 };
    std::size_t total_steps{ 1000 };

    bool run_validation{ true };

    bool use_pml{ true };
    int pml_order{ 3 };
    double pml_kappa_max{ 15.0 };
    double pml_alpha_max{ 0.3 };

    // Derived (computed from above):
    std::size_t size{};
    double c{};
    double dt{};
    std::size_t pml_thickness{};
    double pml_sigma_max{};

    // Default constructor uses hardcoded defaults:
    Simulation_Config() { compute_derived(); }

    // Load from file, falling back to defaults for missing keys:
    [[nodiscard]] static Simulation_Config from_file( std::string const &path ) {
        Simulation_Config cfg{};

        std::ifstream file{ path };
        if ( !file.is_open() ) {
            std::cout << "[Config] File not found: " << path
                      << "; using defaults.\n";
            return cfg;
        }

        auto map{ parse( file ) };

        read( map, "Nx",          cfg.Nx );
        read( map, "Ny",          cfg.Ny );
        read( map, "Nz",          cfg.Nz );
        read( map, "dx",          cfg.dx );
        read( map, "dy",          cfg.dy );
        read( map, "dz",          cfg.dz );
        read( map, "mu",          cfg.mu );
        read( map, "eps",         cfg.eps );
        read( map, "cfl_factor",  cfg.cfl_factor );
        read( map, "total_steps", cfg.total_steps );
        read( map, "run_validation", cfg.run_validation );
        read( map, "use_pml",     cfg.use_pml );
        read( map, "pml_order",   cfg.pml_order );
        read( map, "pml_kappa_max", cfg.pml_kappa_max );
        read( map, "pml_alpha_max", cfg.pml_alpha_max );

        if ( !map.empty() ) {
            std::cout << "[Config] Warning! unknown keys ignored:";
            for ( auto const &[key, val] : map ) {
                std::cout << " " << key;
            }
            std::cout << "\n";
        }

        cfg.compute_derived();
        return cfg;
    }

    [[nodiscard]] std::size_t output_interval() const {
        return total_steps / 100;
    }

private:
    void compute_derived() {
        size          = ( Nx + 1 ) * ( Ny + 1 ) * ( Nz + 1 );
        c             = 1.0 / std::sqrt( mu * eps );
        dt            = cfl_factor / ( c * std::sqrt( 1.0 / ( dx * dx )
                      + 1.0 / ( dy * dy ) + 1.0 / ( dz * dz ) ) );
        pml_thickness = static_cast<std::size_t>( std::cbrt(
                            static_cast<double>( size ) ) / 8.0 );
        pml_sigma_max = 0.8 * ( pml_order + 1 )
                      / ( dx * std::sqrt( mu / eps ) );
    }

    // Because C++ is so ugly:
    using Map = std::unordered_map<std::string, std::string>;

    // Parser:
    [[nodiscard]] static Map parse( std::ifstream &file ) {
        Map map{};
        std::string line{};

        while ( std::getline( file, line ) ) {
            // Strip comments:
            if ( auto pos{ line.find( '#' ) }; pos != std::string::npos ) {
                line.erase( pos );
            }
            // Skip blank lines:
            if ( line.find_first_not_of( " \t\r" ) == std::string::npos ) {
                continue;
            }
            auto eq{ line.find( '=' ) };
            if ( eq == std::string::npos ) {
                continue;
            }
            std::string key{ trim( line.substr( 0, eq ) ) };
            std::string val{ trim( line.substr( eq + 1 ) ) };

            if ( !key.empty() && !val.empty() ) {
                map[key] = val;
            }
        }
        return map;
    }

    [[nodiscard]] static std::string trim( std::string s ) {
        auto const start{ s.find_first_not_of( " \t\r" ) };
        if ( start == std::string::npos ) return {};
        auto const end{ s.find_last_not_of( " \t\r" ) };
        
        return s.substr( start, end - start + 1 );
    }

    // Type-specific readers; erase consumed keys so leftovers can be warned:

    static void read( Map &map, std::string const &key, double &out ) {
        if ( auto it{ map.find( key ) }; it != map.end() ) {
            out = std::stod( it->second );
            map.erase( it );
        }
    }

    static void read( Map &map, std::string const &key, std::size_t &out ) {
        if ( auto it{ map.find( key ) }; it != map.end() ) {
            out = static_cast<std::size_t>( std::stoull( it->second ) );
            map.erase( it );
        }
    }

    static void read( Map &map, std::string const &key, int &out ) {
        if ( auto it{ map.find( key ) }; it != map.end() ) {
            out = std::stoi( it->second );
            map.erase( it );
        }
    }

    static void read( Map &map, std::string const &key, bool &out ) {
        if ( auto it{ map.find( key ) }; it != map.end() ) {
            std::string val{ it->second };
            // Lowercase for comparison:
            for ( auto &ch : val ) ch = static_cast<char>( std::tolower( ch ) );
            out = ( val == "true" || val == "1" );
            map.erase( it );
        }
    }
};