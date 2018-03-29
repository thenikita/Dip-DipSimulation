//#include "MonteCarloApplication.h"
#pragma once

#include <vector>
#include <string>
#include "montecarloapplication.h"

class Particle {
public:
    Particle( double x,
              double y,
              double z,
              double mx,
              double my,
              double mz );

    static double CalculateDipoleInteractionEnergy( double lambda,
                                                    Particle first,
                                                    Particle second );

    static double CalculateInFieldEnergy( double field );

    static double ScalarProduction( std::vector<double> first,
                                    std::vector<double> second );

public:
    double x, y, z;
    double mx, my, mz;
    double d;

public:
    std::string Show( );

    static double CalculateProjection( std::vector<double> first,
                                       std::vector<double> second );

    static double CalculateVectorModule( std::vector<double> vector );

    bool CheckForErrors( double R, double L );
};