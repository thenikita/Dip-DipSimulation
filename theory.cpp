#include "theory.h"
#include <cmath>

double Theory::CalculateMagnetizationForSingleParticle( double field ) {
    return 1 / tanh( field ) - 1 / field;
}

double Theory::CalculateFieldMultiplier(
        double lambda,
        double volumeDensity,
        double magneticMoment) {

    return  4 * lambda * volumeDensity / magneticMoment;
}

double Theory::CalculateMagnetizationForSystem(
        double lambda,
        double field,
        double magneticMoment,
        double volumeDensity) {

    double M =
            CalculateMagnetizationForSingleParticle( field )
            * magneticMoment
            * volumeDensity;

    double H = field /
               CalculateFieldMultiplier(
                       lambda,
                       volumeDensity,
                       magneticMoment );

    double effectiveM = H + M / 3;

    return magneticMoment *
           volumeDensity *
           CalculateMagnetizationForSingleParticle( effectiveM );
}
