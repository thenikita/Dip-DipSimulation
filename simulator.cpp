#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>
#include <thread>

#include "simulator.h"
#include "theory.h"

using std::cin;
using std::cout;
using std::endl;
using std::flush;

using std::mt19937_64;
using std::uniform_real_distribution;

const double PI = 3.14;

Simulator::Simulator(
        double *lambda,
        double *targetVolumeDensity,
        double *particleMagneticMoment,
        unsigned int *particleCount,
        double *particleDiameter,
        double *aspect ) {

    cout << "Simulator is launched!" << endl;

    this->particleDiameter = *particleDiameter;
    this->particleCount = *particleCount;
    this->lambda = *lambda;
    this->targetVolumeDensity = *targetVolumeDensity;
    this->particleMagneticMoment = *particleMagneticMoment;
    this->aspect = *aspect;

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

    ShowSystem();

    //MakeIterations(particleCount);

    //ShowSystem();
}


void Simulator::ShowSystem( ) {
    cout << "\n################# SYSTEM ###################\n";
    cout << "Lambda:           " << this->lambda << endl;
    cout << "Field Module:     " << this->fieldModule << endl;
    cout << "Target Density:   " << this->targetVolumeDensity << endl;
    cout << "Current Density:  " << CalculateCurrentVolumeDensity( ) << endl;

    cout << "Current Length:   " << this->tempLength << endl;
    cout << "Current Radius:   " << this->tempRadius << endl;

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

void Simulator::GenerateParticles( const int particleCount ) {
    // increasing tube size to generate particles in it
    tempRadius = tubeRadius * multiplier;
    tempLength = tubeLength * multiplier;

    double R = tempRadius;
    double L = tempLength;
    double d = particleDiameter;

    int tries = 1;

    cout    << "Generating "
            << particleCount << " particles to tube: "
            << R << " x " << L
            << " with d: " << d
            << endl;

    for ( unsigned i = 0; i < particleCount; i++ ) {
        double x, y, z;

        x = GenerateRandom( -R + d/2,
                            R - d/2,
                            generator );

        R = sqrt( R * R - x * x );
        y = GenerateRandom( -R + d/2,
                            R - d/2,
                            generator );

        z = GenerateRandom( -L / 2.0 + d/2,
                            L / 2.0 - d/2,
                            generator );


        cout << x << " " << y << " " << z << endl;
        Particle temp( x, y, z, 1, 0, 0 );
        if ( !CheckParticleForCollisions( temp, -1 )) {
            this->currentParticles++;
            particles.push_back( temp );

            cout
                 << "Particles ready: "
                 << i + 1 << " in " << tries
                 << particles[i].ToString()
                 << "\n" << std::flush;

            tries = 1;

        } else {
            i--;
            tries++;
            if (tries % 1000 == 0) {

                cout
                     << "Particles ready: "
                     << i << " in " << tries
                     << "\r" << std::flush;
            }
        }
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
        double minX,
        double maxX,
        double minY,
        double maxY,
        double minZ,
        double maxZ) {

    double newX = 0;
    double newY = 0;
    double newZ = 0;

    newX = particle.x +
           GenerateRandom( minX, maxX, generator ) * deltaCoordinate;

    newY = particle.y +
           GenerateRandom( minY, maxY, generator ) * deltaCoordinate;

    newZ = particle.z +
           GenerateRandom( minZ, maxZ, generator ) * deltaCoordinate;

    double left = 1;

    double newMZ = GenerateRandom( -1.0, 1.0, generator );

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
            Particle temp = GenerateDeltaState( particles.at( j ) );

            if ( !CheckParticleForCollisions( temp, j )) {
                double energy =
                        CalculateParticleEnergy( temp, true ) -
                        CalculateParticleEnergy( particles.at( j ), true );

                if ( exp( energy ) > GenerateRandom( 0.0, 1.0, generator )) {
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

    cout
            << "Target:  "
            << tubeRadius << " x "
            << tubeLength << endl;

    bool goOn = true;
    SetStartingTubeSize( );

    int i = 0;
    //TODO: finish the resizing procedure
    while ( goOn ) {
        i++;
        cout
                << "Current: "
                << tempRadius << " x "
                << tempLength << "\r" << flush;

        goOn = CheckIfTubeNeedResize( );

        for ( unsigned j = 0; j < particleCount; j++ ) {
            bool notJumped = true;

            int l = 0;

            /*
            double R = particles[j].x * particles[j].x
                       + particles[j].y * particles[j].y;
            cout << "\nR: " << sqrt( R ) << endl;
            */

            double minX, minY, minZ;
            double maxX, maxY, maxZ;

            if ( tempLength > 0.9 * tubeLength ) {

                if ( particles[j].z > 0 ) {

                    minZ = -1.0;
                    maxZ = 0.0;

                } else {

                    minZ = 0.0;
                    maxZ = 1.0;

                }

            } else {

                minZ = -1.0;
                maxZ = 1.0;

            }

            if ( tempRadius > 0.9 * tubeRadius ) {

                if ( particles[j].x > 0 ) {

                    minX = -1.0;
                    maxX = 0.0;

                } else {

                    minX = 0.0;
                    maxX = 1.0;

                }

                if ( particles[j].y > 0 ) {

                    minY = -1.0;
                    maxY = 0.0;

                } else {

                    minY = 0.0;
                    maxY = 1.0;

                }

            } else {

                minX = -1.0;
                minY = -1.0;

                maxX = 1.0;
                maxY = 1.0;

            }

            while ( notJumped ){

                Particle temp = GenerateDeltaState( particles.at( j ),
                                                    minX,
                                                    maxX,
                                                    minY,
                                                    maxY,
                                                    minZ,
                                                    maxZ );

                //cout << "\ntemp:     " << temp.ToString();

                if ( !CheckParticleForCollisions( temp, j )) {
                    particles.at( j ) = temp;
                    notJumped = false;
                    l = 0;
                }
                l++;

                if ( l > 10 ){
                    notJumped = false;
                }
            }
        }
        /*
        for ( int k = 0; k < particleCount; k++ ) {
            cout << "\njumped to " << particles[k].ToString() << endl;
        }
        cout << endl;
        */

        bool resizeR = true;
        bool resizeL = true;
        if ( tempRadius < tubeRadius ) resizeR = false;
        if ( tempLength < tubeLength ) resizeL = false;
        ResizeTubeIfPossible( resizeR, resizeL );

        /*
        if ( i > 1000 ) {
            goOn = false;
        }
         */
    }

    cout << "\nResize finished! Now will check if there's any errors..." << endl;
    cout << "Finally: " << tempRadius << " x " << tempLength << endl;

    CheckSystemForErrors( );

    if ( tempRadius < tubeRadius ) {

        tempRadius = tubeRadius;
    } else { cout << "ERROR AFTER RESIZING" << endl; }
    if ( tempLength < tubeLength ) {

        tempLength = tubeLength;
    } else { cout << "ERROR AFTER RESIZING" << endl; }

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
            pow(( particleDiameter / 2 ), 3 );

    double tubeVolume = tempRadius * tempRadius * tempLength;

    return particleVolume / tubeVolume;
}


/* Check For Collisions function checks if particle is collided with
   one of unignored particles or walls */
//TODO: IMPORTANT checking for collisions: /walls /particles
bool Simulator::CheckParticleForCollisions(
        const Particle particle,
        int ignored ) {

    //cout << "\nchecking collisions: ";
    //w walls
    double R =
            particle.x * particle.x +
            particle.y * particle.y;

    if (( sqrt( R ) + particleDiameter / 2 ) > tempRadius ){
        /*
        cout <<
             "wall 1: " <<
             sqrt( R ) + particleDiameter / 2 <<
             " " << tempRadius <<
             "\n" << std::flush;
         */
        return true;
    }

    if ( fabs( particle.z ) + particleDiameter / 2 > tempLength / 2 ) {
        /*
        cout <<
             "wall 2: " <<
             fabs( particle.z ) + particleDiameter / 2 <<
             " " << tempLength / 2 <<
             "\n" << std::flush;
        */
        return true;
    }

    //w parts
    for ( unsigned i = 0; i < this->currentParticles; i++ ) {
        if (( i != ignored )) {
            double distance =
                    ( particles.at( i ).x - particle.x ) *
                    ( particles.at( i ).x - particle.x ) +

                    ( particles.at( i ).y - particle.y ) *
                    ( particles.at( i ).y - particle.y ) +

                    ( particles.at( i ).z - particle.z ) *
                    ( particles.at( i ).z - particle.z );

            if ( distance < particleDiameter * particleDiameter ) {
                //cout << "particle\n" << std::flush;
                return true;
            }
        }
    }
    //cout << "none\n" << std::flush;
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

    auto timePast = int( GenerateRandom( 0.0, 1000000000.0, generator ));

    std::this_thread::sleep_for( std::chrono::nanoseconds( timePast ));

    std::chrono::steady_clock::time_point then =
            std::chrono::steady_clock::now( );

    unsigned seed = static_cast<unsigned int>(
            std::chrono::duration_cast<std::chrono::nanoseconds>( then - now ).
                    count( ));

    cout << "RANDOM SEED IS: " << seed
         << " In time: " << timePast << " nsec"
         << endl;

    generator.seed( seed );

    cout <<  "Test sequence:" << endl;
    for ( int i = 0; i < 10; ++i ) {
        cout << GenerateRandom(0, 1, generator) << endl;
    }
    cout << endl;
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
void Simulator::CollectTubeSizes( unsigned num ) {
    double R = sqrt(
            particles.at( num ).x * particles.at( num ).x +
            particles.at( num ).y * particles.at( num ).y );

    double L = fabs( particles.at( num ).z ) * 2;

    if ( R < this->tempRadius ) this->tempRadius = R;
    if ( L < this->tempLength ) this->tempLength = L;
}


void Simulator::ResizeTubeIfPossible( bool resizeR = true,
                                      bool resizeL = true) {
    double R = 0.0;
    double L = 0.0;

    for ( int i = 0; i < particleCount; i++ ) {
        double r = particles[i].x * particles[i].x
                   + particles[i].y * particles[i].y;

        r = sqrt( r );

        if ( fabs( r ) > R ) {

            R = fabs( r );
        }

        if ( fabs( particles[i].z ) > L / 2 ) {

            L = fabs( particles[i].z * 2 );
        }
    }

    R += particleDiameter / 2;
    L += particleDiameter;

    //cout << "L: " << L << " R: " << R << endl;
    if (( L < tempLength ) and resizeL ) {
        if ( L > 0.9 * tubeLength ) {
            /*
            cout << "resize L to " <<
                 L <<
                 " from " <<
                 tempLength << endl;
            */
            tempLength = L;
        }
    }

    if (( R < tempRadius ) and resizeR ) {
        if ( R > 0.9 * tubeRadius ) {
            /*
            cout << "resize R to " <<
                 R <<
                 " from " <<
                 tempRadius << endl;
            */
            tempRadius = R;
        }
    }
}


// TODO not sure need resize checker works correct
bool Simulator::CheckIfTubeNeedResize( ) {
    double R = tempRadius;
    double L = tempLength;

    return (( R > tubeRadius ) || ( L > tubeLength ));
}

void Simulator::GenerateTube(
        unsigned int particleCount,
        double particleDiameter,
        double targetVolumeDensity,
        double &tubeRadius,
        double &tubeLength,
        double aspect ) {

    cout << "\n*************************************************************\n";
    cout << "Now going to generate tube..." << endl;

    Particle::SetDiameter(particleDiameter);
    cout << "Particle Diameter is set to "
         << particleDiameter
         << endl;

    // N * 4/3 pi r^3
    double allParticlesVolume =
            particleCount * 4.0 / 3.0 * PI * pow( particleDiameter, 3 ) / 8.0;
    double tubeVolume = allParticlesVolume / targetVolumeDensity;

    // V = pi r^2 * L = pi r^3 * aspect
    tubeRadius = pow( tubeVolume / PI / aspect, 0.3333333 );

    if ( tubeRadius < particleDiameter ) {
        tubeRadius = 0.55 * particleDiameter;
    }

    tubeLength = tubeRadius * aspect;

    cout    << "New generated tube: "
            << "R = " << tubeRadius
            << " and L = " << tubeLength
            << "\n"
            << endl;
}

void Simulator::SetStartingTubeSize( ) {

    double L = 0;
    double R = 0;

    for ( int i = 0; i < particleCount; i++ ) {

        double l = 2 * fabs(particles[i].z) + particleDiameter;

        double r = particles[i].x * particles[i].x +
                   particles[i].y * particles[i].y;

        r = sqrt( r ) + particleDiameter / 2;

        if ( l > L ) L = l;
        if ( r > R ) R = r;
    }

    tempRadius = R;
    tempLength = L;

}
