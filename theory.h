#pragma once

class Theory {
public:
    static double CalculateMagnetizationForSingleParticle( double field );

    static double CalculateFieldMultiplier( double lambda,
                                            double volumeDensity,
                                            double magneticMoment );

    static double CalculateMagnetizationForSystem( double lambda,
                                                   double field,
                                                   double magneticMoment,
                                                   double volumeDensity );
};
