#pragma once

#include "../Config/config.hpp"

class Grid;

// Abstract base class for all soruces:
class Source {
public:
    virtual ~Source() = default;
    virtual void apply( Grid &grid, double time_step ) = 0;
};

// Straight wire with sinusoidal current along x:
class Straight_Wire_X : public Source {
private:
    double amplitude_;
    double frequency_;
    std::size_t y_, z_;
    std::size_t x_start_, x_end_;

public:
    Straight_Wire_X(
        double amp, 
        double freq,
        std::size_t y, std::size_t z,
        std::size_t x_start, std::size_t x_end )
        : amplitude_{ amp }
        , frequency_{ freq }
        , y_{ y }
        , z_{ z }
        , x_start_{ x_start }
        , x_end_{ x_end }
        { }

    void apply( Grid &grid, double time_step ) override;
};

// Hard source injection at single point:
class Point_Source : public Source { 
private:
    double value_;
    std::size_t x_, y_, z_;

public:
    Point_Source(
        double value,
        std::size_t x, std::size_t y, std::size_t z )
        : value_{ value }
        , x_{ x }
        , y_{ y }
        , z_{ z }
        { }

    void apply( Grid &grid, double time_step ) override;
};

// Gaussian pulse current source:
class Gaussian_Pulse : public Source { 
private:
    double amplitude_;
    double t_0_;
    double width_;
    std::size_t x_, y_, z_;

public:
    Gaussian_Pulse(
        double amp,
        double t_0, 
        double width,
        std::size_t x, std::size_t y, std::size_t z )
        : amplitude_{ amp }
        , t_0_{ t_0 }
        , width_{ width }
        , x_{ x }, y_{ y }, z_{ z }
        { }

    void apply( Grid &grid, double time_step ) override;
};