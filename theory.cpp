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
        double particleDensity) {

    double M =
            CalculateMagnetizationForSingleParticle( field )
            * magneticMoment
            * particleDensity;

    double newM = M / 3 * CalculateFieldMultiplier(lambda,
                                               particleDensity,
                                               magneticMoment);

    double effectiveM = field + newM;

    return magneticMoment *
           particleDensity *
           CalculateMagnetizationForSingleParticle( effectiveM );
}
