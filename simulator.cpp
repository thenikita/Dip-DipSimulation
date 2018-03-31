#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <thread>

#include "simulator.h"
#include "montecarloapplication.h"
#include "theory.h"

using std::cin;
using std::cout;
using std::endl;

using std::mt19937_64;
using std::uniform_real_distribution;

const double PI = 3.14;

Simulator::Simulator(
        double *lambda,
        double *targetVolumeDensity,
        double *particleMagneticMoment,
        unsigned int *particleCount,
        double *particleDiameter,
        double *aspect,
        double &targetTubeR,
        double &targetTubeL ) {

    cout << "Simulator is launched!" << endl;

    this->particleDiameter = *particleDiameter;
    this->particleCount = *particleCount;
    this->lambda = *lambda;
    this->targetVolumeDensity = *targetVolumeDensity;
    this->particleMagneticMoment = *particleMagneticMoment;
    this->aspect = *aspect;
    this->tubeLength = targetTubeR;
    this->tubeRadius = targetTubeL;

    tempRadius = tubeRadius * multiplier;
    tempLength = tubeLength * multiplier;

    cout << "Enter the field (xi) vector (x, y, z): " << endl;

    double tempField;
    for ( int i = 0; i < 3; i++ ) {
        cin >> tempField;
        field.push_back( tempField );
    }

    this->fieldModule = Particle::CalculateVectorModule(field);

    cout << "Field is set to: (";
    cout << field.at( 0 ) << ", "
                             << field.at( 1 ) << ", "
                                                 << field.at( 2 ) << ")"
                                                                     << endl;

    cout << "Field module: " << fieldModule << endl;
    cout << "\n*************************************************************\n";
}

void Simulator::Run( ) {
    GenerateTube(particleCount,
                 particleDiameter,
                 targetVolumeDensity,
                 tubeRadius,
                 tubeLength,
                 aspect);

    SetGeneratorRandomSeed( );

    GenerateParticles( this->particleCount );

    CheckSystemForErrors( );

    ShowSystem();

    MakeResizing();

    //MakeIterations(particleCount);

    ShowSystem();
}


void Simulator::ShowSystem( ) {
    cout << "\n################# SYSTEM ###################\n";
    cout << "Lambda:           " << this->lambda << endl;
    cout << "Field Module:     " << this->fieldModule << endl;
    cout << "Target Density:   " << this->targetVolumeDensity << endl;
    cout << "Current Density:  " << CalculateCurrentVolumeDensity( ) << endl;

    cout << "Particles:        " << particleCount << endl;
    cout << "  N  |"
            "            COORDINATES             |"
            "            MAGNETIC MOMENT       "
         << endl;

    for ( unsigned i = 0; i < particleCount; i++ ) {
        cout << "  " << i << "    ";
        cout << particles.at( i ).ToString( ) << endl;
    }
    cout << "\n################# SYSTEM ###################\n";
}

void Simulator::GenerateParticles( const int particleAmount ) {
    cout    << "Generating "
            << particleAmount << " particles to tube: "
            << tubeRadius << " x " << tubeLength
            << endl;

    Particle::SetDiameter(particleDiameter);
    cout << "Particle Diameter is set to "
            << particleDiameter
            << endl;

    cout << "Generating particles..." << endl;
    // TODO fix crash on here (1 particle)
    for ( int i = 0; i < particleAmount; i++ ) {
        double x, y, z;
        double R = tubeRadius;

        x = GenerateRandom( -R + 1,
                            R - 1,
                            generator );

        R = sqrt( R * R - x * x );
        y = GenerateRandom( -R + 1,
                            R - 1,
                            generator );

        z = GenerateRandom( -tubeLength / 2 + 1,
                            tubeLength / 2 - 1,
                            generator );

        Particle temp( x, y, z, 1, 0, 0 );
        if ( !CheckParticleForCollisions( temp, -1 )) {
            this->particleCount++;
            particles.push_back( temp );

            cout << "Particles ready: " << i + 1 << "\r" << std::flush;

        } else { i--; }
    }

    cout << "\nParticles've been generated!" << endl;
}


/* Generate Delta State function generates new state of i-th particle
   which is based on current position in the tube and using the
   random number generator*/
//
// IMPORTANT MEMO
// Here we calculate energy prodivided by kT so later we can
// use it as is in transition calculations without more dividings.
// 

Particle Simulator::GenerateDeltaState(
        const Particle particle,
        bool ifNeedResize ) {

    double newX = 0;
    double newY = 0;
    double newZ = 0;

    if ( ifNeedResize ) {

        if ( particle.x < 0 ) {

            newX = particle.x +
                   GenerateRandom( 0, 1, generator ) * deltaCoordinate;

        } else if ( particle.x > 0 ) {

            newX = particle.x +
                   GenerateRandom( -1, 0, generator ) * deltaCoordinate;

        }

        if ( particle.y < 0 ) {

            newY = particle.y +
                   GenerateRandom( 0, 1, generator ) * deltaCoordinate;

        } else if ( particle.y > 0 ) {

            newY = particle.y +
                   GenerateRandom( -1, 0, generator ) * deltaCoordinate;

        }

        if ( particle.z < 0 ) {

            newZ = particle.z +
                   GenerateRandom( 0, 1, generator ) * deltaCoordinate;

        } else if ( particle.z > 0 ) {

            newZ = particle.z +
                   GenerateRandom( -1, 0, generator ) * deltaCoordinate;

        }

    } else {

        newX = particle.x +
               GenerateRandom( -1, 1, generator ) * deltaCoordinate;

        newY = particle.y +
               GenerateRandom( -1, 1, generator ) * deltaCoordinate;

        newZ = particle.z +
               GenerateRandom( -1, 1, generator ) * deltaCoordinate;

    }

    double left = 1;

    double newMZ = GenerateRandom( -1, 1, generator );

    left -= newMZ * newMZ;
    double newMY = GenerateRandom( -sqrt( left ), sqrt( left ), generator );

    left -= newMY * newMY;
    //TODO: check whether we need +- x generation
    double newMX = sqrt( left );

    auto *newState = new Particle(
            newX,
            newY,
            newZ,
            newMX,
            newMY,
            newMZ);

    return *newState;
}


/* Make Iterations function is making iterations moving particles
   inside the tube. Use after resizing the tube to relevant size
   if wanna get results accounting demagnetisation and other shit*/

//TODO: add deltaCoordinate decreasing in future iterations to immprove speed
void Simulator::MakeIterations( int particleAmount ) {
    cout << "\n****************ITERATIONS******************\n";

    double averageProjection = 0;

    //check it
    for ( unsigned i = 0; i < stepsAmount; i++ ) {
        cout << 100 * i / stepsAmount << " %  Tube: " << tubeRadius << " " << tubeLength << "   \r"
             << std::flush;

        for ( unsigned j = 0; j < particleAmount; j++ ) {
            Particle temp = GenerateDeltaState( particles.at( j ), false );

            if ( !CheckParticleForCollisions( temp, j )) {
                double energy =
                        CalculateParticleEnergy( temp, true ) -
                        CalculateParticleEnergy( particles.at( j ), true );

                if ( exp( energy ) > GenerateRandom( 0, 1, generator )) {
                    particles.at( j ) = temp;
                }

                averageProjection += particles.at( j ).mz;
            }
        }
    }

    cout
            << "100 %  Tube: " << tubeRadius << " "
            << tubeLength << "   \r"

            << "\nIterations've finished!" << endl

            << "Average is:           "
            << averageProjection / stepsAmount / particleAmount
            << endl

            << "Theory for isolated:  "
            << Theory::CalculateMagnetizationForSingleParticle( fieldModule )
            << endl

            << "Theory for system:    "
            << Theory::CalculateMagnetizationForSystem( lambda,
                                                        fieldModule,
                                                        particleMagneticMoment,
                                                        targetVolumeDensity )
            << endl;

    CheckSystemForErrors( );

    cout << "****************ITERATIONS******************\n";
}


/* Make Resizing function performs resizing of the tube to get particles'
   volume density be equal to target*/
void Simulator::MakeResizing( ) {
    cout << "\n*****************RESIZING*******************\n";

    bool goOn = true;

    //TODO: finish the resizing procedure
    while ( goOn ) {

        tempRadius = tubeRadius;
        tempLength = tubeLength;

        cout
                << " Tube: "
                << tubeRadius << " "
                << tubeLength << " "
                << " Max: "
                << tempRadius << " "
                << tempLength << " " << endl;

        goOn = CheckIfTubeNeedResize( );

        for ( unsigned j = 0; j < particleCount; j++ ) {
            Particle temp = GenerateDeltaState( particles.at( j ), goOn );

            if ( !CheckParticleForCollisions( temp, j )) {
                double energy =
                        CalculateParticleEnergy( temp, true ) -
                        CalculateParticleEnergy( particles.at( j ), true );

                if ( exp( energy ) > GenerateRandom( 0, 1, generator )) {
                    particles.at( j ) = temp;
                }
            }
        }

        ResizeTubeIfPossible( );

        /*
        k++;

        cout
            << " Tube: "
            << tubeRadius << " "
            << tubeLength << " "
            << " Max: "
            << tempRadius << " "
            << tempLength << "          \n" << std::flush;

        for (int i = 0; i < particleCount; i++)
        {
            Particle temp = GenerateDeltaState(particles.at(i));

            if (!CheckParticleForCollisions(temp, i))
            {
                double energy =
                    CalculateParticleEnergy(temp, false) -
                    CalculateParticleEnergy(particles.at(i), false);

                if (exp(energy) > GenerateRandom(0, 1, generator))
                {
                    particles.at(i) = temp;
                    CollectTubeSizes(i);
                }
            }

        }
*/
    }

    cout << "Resize finished! Now will check if there's any errors..." << endl;
    CheckSystemForErrors( );

    cout << "\n*****************RESIZING*******************\n";
}


/* Calculate Particle Energy function calculates energy of particles 
   inside the outer magnetic fieldModule and in dipole-dipole interaction */
double Simulator::CalculateParticleEnergy(
        const Particle particle,
        bool mode ) {

    //here's the fieldModule part of energy
    std::vector<double> point;
    point.push_back( particle.x );
    point.push_back( particle.y );
    point.push_back( particle.z );

    std::vector<double> outerField = GetVectorField( point, mode );

    std::vector<double> moment;
    moment.push_back( particle.mx );
    moment.push_back( particle.my );
    moment.push_back( particle.mz );

    double outerFieldEnergy =
            Particle::CalculateProjection( moment, outerField );

    double dipoleInteractionEnergy = 0;


    //TODO: calculate dipole interaction energy
    //Already have dip-dip function in particle obj

    double totalEnergy = dipoleInteractionEnergy + outerFieldEnergy;

    return totalEnergy;
}

double Simulator::CalculateCurrentVolumeDensity( ) {
    double particleVolume =
            double(particleCount) *
            4 / 3 *
            PI *
            pow(( particleDiameter / 2 ), 3 );

    double tubeVolume = PI * tubeRadius * tubeRadius * tubeLength;

    return particleVolume / tubeVolume;
}


/* Check For Collisions function checks if particle is collided with
   one of unignored particles or walls */
//TODO: IMPORTANT checking for collisions: /walls /particles
bool Simulator::CheckParticleForCollisions(
        const Particle particle,
        int ignored ) {

    //w walls
    double R =
            particle.x * particle.x +
            particle.y * particle.y;

    if (( sqrt( R ) + particleDiameter / 2 ) > tubeRadius ) return true;
    if ( fabs( particle.z ) + particleDiameter / 2 > tubeLength / 2 ) return true;

    //w parts
    for ( unsigned i = 0; i < this->particleCount; i++ ) {
        if (( i != ignored )) {
            double distance =
                    ( particles.at( i ).x - particle.x ) * ( particles.at( i ).x - particle.x ) +
                    ( particles.at( i ).y - particle.y ) * ( particles.at( i ).y - particle.y ) +
                    ( particles.at( i ).z - particle.z ) * ( particles.at( i ).z - particle.z );

            if ( distance < particleDiameter * particleDiameter ) return true;
        }
    }

    return false;
}


/* Generate Random function just generates random double value */
double Simulator::GenerateRandom( double min, double max, mt19937_64 &generator ) {
    uniform_real_distribution<double> distr( min, max );
    return distr( generator );
}


/* Set Random Seed function sets seed for random generator */
void Simulator::SetGeneratorRandomSeed( ) {
    generator.seed( 1 );

    std::chrono::steady_clock::time_point now =
            std::chrono::steady_clock::now( );

    auto timePast = int( GenerateRandom( 0, 1000000000, generator ));

    std::this_thread::sleep_for( std::chrono::nanoseconds( timePast ));

    std::chrono::steady_clock::time_point then =
            std::chrono::steady_clock::now( );

    unsigned seed = static_cast<unsigned int>(
            std::chrono::duration_cast<std::chrono::nanoseconds>( then - now ).
                    count( ));

    cout << "RANDOM SEED IS: " << seed
         << " In time: " << timePast << "nsecs."
         << endl;

    generator.seed( seed );
}


/* Get Vector Field function returns outer fieldModule vector in two
   modes: normal and radial */
std::vector<double> Simulator::GetVectorField( std::vector<double> point, bool mode ) {
    // true = normal to z oriented fieldModule
    // false = center oriented fieldModule
    if ( mode ) {
        return field;
    } else {
        std::vector<double> temp;
        temp.push_back( -point.at( 0 ));
        temp.push_back( -point.at( 1 ));
        temp.push_back( -point.at( 2 ));

        return temp;
    }
}


/* Check System For Errors is backup function which checks system
   just one more time for collisions and maybe other errors */
bool Simulator::CheckSystemForErrors( ) {
    cout << "Checking system for errors in " << this->particleCount << endl;
    for ( unsigned i = 0; i < this->particleCount; i++ ) {
        //cout << i << std::flush;
        for ( unsigned j = i + 1; j < this->particleCount; j++ ) {
            //cout << j << std::flush;

            double distance =
                    pow( particles.at( i ).x - particles.at( j ).x, 2 ) +
                    pow( particles.at( i ).y - particles.at( j ).y, 2 ) +
                    pow( particles.at( i ).z - particles.at( j ).z, 2 );

            if ( distance < particleDiameter * particleDiameter ) {
                cout << "\nERROR on " << i << ", " << j << endl;
                return true;
            }
        }
        //cout << endl;
    }
    cout << "NO ERRORS" << endl;
    return false;
}


/* Collect Tube Sizes function checks if it's possible to remember
   new smallest tube's parameters */
void Simulator::CollectTubeSizes( int num ) {
    double R = sqrt(
            particles.at( num ).x * particles.at( num ).x +
            particles.at( num ).y * particles.at( num ).y );

    double L = abs( particles.at( num ).z ) * 2;

    if ( R < this->tempRadius ) this->tempRadius = R;
    if ( L < this->tempLength ) this->tempLength = L;
}


void Simulator::ResizeTubeIfPossible( ) {
    if ( tempLength < tubeLength ) tubeLength = tempLength;
    if ( tempRadius < tubeRadius ) tubeRadius = tempRadius;
}


// TODO not sure need resize checker works correct
bool Simulator::CheckIfTubeNeedResize( ) {
    double R = tempRadius;
    double A = aspect;

    return !(( R > tubeRadius ) && ( R * A > tubeLength ));
}

void Simulator::GenerateTube(
        unsigned int particleCount,
        double particleDiameter,
        double targetVolumeDensity,
        double &targetTubeR,
        double &targetTubeL,
        double aspect ) {

    cout << "\n*************************************************************\n";
    cout << "Now going to generate tube..." << endl;

    // N * 4/3 pi r^3
    double allParticlesVolume =
            particleCount * 4.0 / 3.0 * PI * pow( particleDiameter, 3 ) / 8.0;
    double tubeVolume = allParticlesVolume / targetVolumeDensity;

    // V = pi r^2 * L = pi r^3 * aspect
    targetTubeR = pow( tubeVolume / PI / aspect, 0.3333333 );

    if ( targetTubeR < particleDiameter ) {
        targetTubeR = 0.55 * particleDiameter;
    }

    targetTubeL = targetTubeR * aspect;

    cout    << "New generated tube: "
            << "R = " << targetTubeR
            << " and L = " << targetTubeL
            << endl;
}