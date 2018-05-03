//#include "MonteCarloApplication.h"
#pragma once

#include <vector>
#include <string>
#include "montecarloapplication.h"

class Particle {
public:
    double x, y, z;
    double mx, my, mz;
    static double d;

public:
    Particle( double x,
              double y,
              double z,
              double mx,
              double my,
              double mz );

    static double CalculateDipoleInteractionEnergy(
            double lambda,
            double diameter,
            Particle first,
            Particle second );

    static double CalculateInFieldEnergy(std::vector<double> magneticMoment,
                                         std::vector<double> fieldVector );

    static double ProductScalars(
            std::vector<double> first,
            std::vector<double> second );

    std::string ToString( );

    static double CalculateProjection(
            std::vector<double> first,
            std::vector<double> second );

    static double CalculateVectorModule( std::vector<double> vector );

    bool CheckForErrors( double R, double L );

    static void SetDiameter( double diameter ) { d = diameter; }
};