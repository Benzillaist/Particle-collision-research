#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

//declaration of tweakable variables, keep tickrate very high and don't add too many particles or increase the radius by too much
double std_acceleration = 1;
int numberOfParticles = 10;
double minVelocity = 5;
double maxVelocity = 10;
double trackLength = 10;
double tickRate = 10000;
double particleRadius = 0.02;
double tickFrequency = 1/tickRate;

//Standard particle class for storing particle data
class Particle {
    public:

    //each particle consists of the ideal velocity and starting position, which will stay static, and then the current position and velocity, which will be changed during the simulation
        double idealvelocity;
        double initialposition;
        double currentvelocity;
        double currentposition;
        Particle(double idealv, double initialp) { //constructor function
            idealvelocity = idealv;
            initialposition = initialp;
        }
};

vector<Particle> particleList;

//integer modulo function that handles negatives >:(
double nmod(double dividend, double divisor) {
    if(dividend < 0) {
        return (dividend + divisor);
    }
    return fmod(dividend, divisor);
};

//finds the initial spacing between particles
double findSpacing() {
    return trackLength/numberOfParticles;
};

//returns a random velocity between the min and max velocities specified above
double randVelocity() {
    double cutRand = (double)rand() / 32767;
    return (double)minVelocity + (((double)maxVelocity - (double)minVelocity) * cutRand);
};

//returns a single particle object
Particle generateParticle(double idealvelocity, double initialposition) {
    Particle particle = Particle(idealvelocity, initialposition);
    return particle;
};

//generates the particles and adds them to the particle list
void loadParticles() {
    for(int i = 0; i < numberOfParticles; i++) {
        double rV = randVelocity();
        particleList.push_back(generateParticle(rV, i*findSpacing()));
        particleList[i].currentposition = particleList[i].initialposition;
        particleList[i].currentvelocity = particleList[i].idealvelocity;
    }
};


//outdated collision management method
void didCollide(int i) {
    
    //estimates where the particle will be after one tick and checks if it steps into the bounding box of another particle
    double pf0 = fmod(particleList[nmod((i-1),numberOfParticles)].currentposition + (double)(particleList[nmod((i-1),numberOfParticles)].currentvelocity*tickFrequency + particleRadius), trackLength);
    double pf1l = fmod(particleList[i].currentposition + (double)(particleList[i].currentvelocity*tickFrequency - particleRadius), trackLength);
    double pf1m = fmod(particleList[i].currentposition + (double)(particleList[i].currentvelocity*tickFrequency + particleRadius), trackLength);
    double pf2 = fmod(particleList[nmod((i+1),numberOfParticles)].currentposition + (double)(particleList[nmod((i+1),numberOfParticles)].currentvelocity*tickFrequency - particleRadius), trackLength);
    
    //compares positions
    bool comparePiAndMore = pf1m < pf2;
    bool comparePiAndLess = pf1l < pf0;
    bool comparePMoreAndLess = pf2 < pf0;

    //checks if compared positions are correct
    if(!(((comparePiAndMore) && (comparePiAndLess) && (comparePMoreAndLess)) || ((comparePiAndMore) && (!comparePiAndLess) && (!comparePMoreAndLess)) || ((!comparePiAndMore) && (!comparePiAndLess) && (comparePMoreAndLess)))) {
        //printf("yup: %d\n", i);
        double pv1 = particleList[i].currentvelocity;
        double pv2 = particleList[nmod((i+1),numberOfParticles)].currentvelocity;
        particleList[i].currentvelocity = pv2;
        particleList[nmod((i+1),numberOfParticles)].currentvelocity = pv1;
        //didCollide(nmod((i-1),numberOfParticles));
        //didCollide(nmod((i+1),numberOfParticles));
    };
};

//updated collision management function
void didCollideV2() {
    double x0, x1, distance, v0, v1;
    for(int i = 0; i < numberOfParticles; i++) { //for each particle and the succeeding one:

        //get both particle's positions
        x0 = particleList[i].currentposition;
        x1 = particleList[nmod(i+1, numberOfParticles)].currentposition;

        //finds the distance between
        distance = x1 - x0; 

        //get each particles velocity
        v0 = particleList[i].currentvelocity;
        v1 = particleList[nmod(i+1, numberOfParticles)].currentvelocity;

        if(distance <= (-2 * particleRadius)) {
            distance += trackLength; //if the distance is negative, it means that the succeeding particle has looped around, so i add on to the distance to account for that change
        }
        if((abs(distance) < 2 * particleRadius) && (v0 > v1)) { //if the particles are overlapping and are coming closer together, the velocities switch
            //swaps velocities of the two particles
            particleList[i].currentvelocity = v1;
            particleList[nmod((i+1),numberOfParticles)].currentvelocity = v0;
        }
    }
}

//function that accelerates particles back to their ideal velocity if the particle is not at previously states ideal velocity
void accelerateParticles() {
    for(int i = 0; i < numberOfParticles; i++) {

        //updates position of each of the particles
        particleList[i].currentposition = fmod(particleList[i].currentposition + (particleList[i].currentvelocity*tickFrequency), trackLength);

        //if the difference in between the current and ideal velocities is equal to or lower than the change that would occur if the particle is accelerated, the current velocity is set equal to the ideal velocity
        if(abs(particleList[i].idealvelocity - particleList[i].currentvelocity) <= (std_acceleration*tickFrequency)) {
            particleList[i].currentvelocity = particleList[i].idealvelocity;
        }

        //if the current velocity is greater than the ideal velocity, then the current velocity is decreased by the acceleration every tick
        else if((particleList[i].idealvelocity - particleList[i].currentvelocity) < 0) {
            particleList[i].currentvelocity -= (std_acceleration*tickFrequency);
        }

        //if the current velocity is less than the ideal velocity, then the current velocity is increased by the acceleration every tick
        else if((particleList[i].idealvelocity - particleList[i].currentvelocity) > 0) {
            particleList[i].currentvelocity += (std_acceleration*tickFrequency);
        };
    }
};

//old method that advanced the particles one tick forwards
void iterate() {

    for(int i = 0; i < numberOfParticles; i++) {
        didCollide(i);
    };
    accelerateParticles();
}

//figures out if any of the particles are in the wrong position   
void areParticlesMessedUp() {
    int pf0, pf1, pf2;
    for(int i = 0; i < numberOfParticles; i++) {

        //gets position of each particle, the one in front of it, and the one behind it
        pf0 = particleList[nmod((i-1),numberOfParticles)].currentposition;
        pf1 = particleList[i].currentposition;
        pf2 = particleList[nmod((i+1),numberOfParticles)].currentposition;

        //compares the positions of each of the particles with each other
        bool comparePiAndMore = pf1 < pf2;
        bool comparePiAndLess = pf1 < pf0;
        bool comparePMoreAndLess = pf2 < pf0;

        //checks if the combination of the above comparisons is one of the three legal permutations, and if not, makes you aware of it
        if(!(((comparePiAndMore) && (comparePiAndLess) && (comparePMoreAndLess)) || ((comparePiAndMore) && (!comparePiAndLess) && (!comparePMoreAndLess)) || ((!comparePiAndMore) && (!comparePiAndLess) && (comparePMoreAndLess)))) {
            printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH: %d\n", i);
        };
    };
};

//DATA ANALYSIS METHODS

double radiusTolerance = 5;

//searches through the particle list and identifies the one with the Lowest Ideal Velocity
int lowestIV() {
    double lV = particleList[0].idealvelocity;
    double tV;
    int li = 0;
    for(int i = 1; i < numberOfParticles; i++) {
        tV = particleList[i].idealvelocity;
        if(tV<lV) {
            lV = tV;
            li = i;
        }
    }
    return li;
};

//finds the front of the chain of particles if the entire collection of particles is together
int findFront() {
    int i = lowestIV(); //first index
    while(((particleList[nmod(i+1, numberOfParticles)].currentposition-particleList[i].currentposition) < (radiusTolerance * particleRadius)) && ((particleList[nmod(i+1, numberOfParticles)].currentposition-particleList[i].currentposition)>=0)) {
        i = nmod(i+1, numberOfParticles);
    }
    return i;
};

//tests if every particle is in the same chain
bool fullHouse() {
    int i = nmod(findFront()+1, numberOfParticles); //first index
    int nP = 0; //number of particles in loop
    while(((particleList[nmod(i+1, numberOfParticles)].currentposition-particleList[i].currentposition) < (radiusTolerance * particleRadius)) && ((particleList[nmod(i+1, numberOfParticles)].currentposition-particleList[i].currentposition)>=0)) {
        i = nmod(i+1, numberOfParticles);
        nP++;
    }
    if(nP == numberOfParticles-1) {
        return 1;
    }
    return 0;
}

//finds the mean velocity of all of the particles
double meanVelocity() {
    double v_sum = 0;

    for(int i = 0; i < numberOfParticles; i++) {
        v_sum += particleList[i].currentvelocity;
    };

    return (v_sum / numberOfParticles);
};

//finds the difference from the ideal velocity of each particle to the mean velocity of all of the particles current velocities
void meanResiduals() {
    double meanvelocity = meanVelocity();
    double lowestResidual = abs(particleList[0].idealvelocity - meanvelocity);
    int lowestIndex;
    double individualResidual;
    printf("\nMean residuals:\n");
    for(int i = 0; i < numberOfParticles; i++) {
        individualResidual = abs(particleList[i].idealvelocity - meanvelocity);
        if (individualResidual < lowestResidual) {
            lowestResidual = individualResidual;
            lowestIndex = i;
        };
        printf("%d: %lf\n", i, individualResidual);
    };
    printf("Lowest residual: %d: %lf\n", lowestIndex, lowestResidual);
}

//prints the current stats on all of the particles
void printCurrentParticleData() {
    for(int i = 0; i < numberOfParticles; i++) {
        printf("%d: Current position: %lf, Current velocity: %lf\n", i, particleList[i].currentposition, particleList[i].currentvelocity);
    };
};

//prints the intial stats of all the particles
void printInitialParticleData() {
    for(int i = 0; i < numberOfParticles; i++) {
        printf("%d: Initial position: %lf, Ideal velocity: %lf\n", i, particleList[i].initialposition, particleList[i].idealvelocity);
    };
};



//max number of iterations for each run, the total seconds that elapse is equal to iterationNumber/tickRate
int iterationNumber = 5000000;

//number of runs that will be ran
int numberOfRuns = 1;

int main(void) {
    int isFrontRun = 0;
    int completedRuns = 0;

    //runs the number of previously specified simulations
    for(int j = 0; j < numberOfRuns; j++) {

        //resets the particle list from any past runs
        particleList.clear();

        //resets the random number generator "seed" so I don't get the same particles each time
        srand(time(0));

        //prints the "seed", so each run can be re-ran at a later time
        printf("Seed: %d\n",rand());

        //generates all the particles based off the previous random number generation
        loadParticles();

        //prints initial particle data
        printCurrentParticleData();

        //runs for a maximum of cycles that was previously defined
        for(int i = 0; i <= iterationNumber; i++) {

            //manages and collisions that would have happened for each particle
            didCollideV2();

            //moves and accelerates each particle
            accelerateParticles();

            //checks every 1000 ticks to see if the particles have all bunched up, if so, it prints the particle data and then leaves the loop
            if((i%1000) == 0) {


                if(fullHouse() == 1) {
                    printf("\nIteration number: %d\n", i);
                    printCurrentParticleData();
                    completedRuns++;
                    break;    
                }
            }

            //since printing out each tick would clog the console and not be useful for debugging or analyzing, I only print out every 10000 ticks
            if((i%1000) == 0) {
                printf("\nIteration number: %d\n", i);
                printCurrentParticleData();
            }
        }

        //general methods that are useful for analyzing data and debugging
        printf("\nInitial particle data: \n");
        printInitialParticleData();
        printf("\nMean velocity: %lf", meanVelocity());
        meanResiduals();
        printf("Lowest ideal velocity: %d: %lf\n", lowestIV(), particleList[lowestIV()].idealvelocity);
        printf("Start of chain: %d\n", findFront());

        //if the particle with the lowest velocity is the one at the front of the chain, it adds one to the counter
        if((findFront() == lowestIV()) && (fullHouse() == 1)) {
            isFrontRun++;
        }
    }

    //prints the number of successful runs to the number of runs where the particle with the lowest velocity was at the front
    printf("Total completed runs: %d, of which, %d had the particle with the lowest ideal velocity at the front\n", completedRuns, isFrontRun);
};
