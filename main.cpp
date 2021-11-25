#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

double std_acceleration = 1;
int numberOfParticles = 16;
double minVelocity = 5;
double maxVelocity = 10;
double trackLength = 10;
double tickRate = 1000;
double particleRadius = 0.01;
double tickFrequency = 1/tickRate;

//Standard particle class for storing particle data
class Particle {
    public:
        double idealvelocity;
        double initialposition;
        double currentvelocity;
        double currentposition;
        Particle(double idealv, double initialp) {
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
        printf("yup: %d\n", i);
        double pv1 = particleList[i].currentvelocity;
        double pv2 = particleList[nmod((i+1),numberOfParticles)].currentvelocity;
        particleList[i].currentvelocity = pv2;
        particleList[nmod((i+1),numberOfParticles)].currentvelocity = pv1;
        //didCollide(nmod((i-1),numberOfParticles));
        //didCollide(nmod((i+1),numberOfParticles));
    };
};

void accelerateParticles() {
    for(int i = 0; i < numberOfParticles; i++) {
        particleList[i].currentposition = fmod(particleList[i].currentposition + (particleList[i].currentvelocity*tickFrequency), trackLength);
        if(abs(particleList[i].idealvelocity - particleList[i].currentvelocity) <= (std_acceleration*tickFrequency)) {
            particleList[i].currentvelocity = particleList[i].idealvelocity;
        }
        else if((particleList[i].idealvelocity - particleList[i].currentvelocity) < 0) {
            particleList[i].currentvelocity -= (std_acceleration*tickFrequency);
        }
        else if((particleList[i].idealvelocity - particleList[i].currentvelocity) > 0) {
            particleList[i].currentvelocity += (std_acceleration*tickFrequency);
        };
    }
};

//advances the particles one tick forwards
void iterate() {

    for(int i = 0; i < numberOfParticles; i++) {
        didCollide(i);
    };
    accelerateParticles();
}

//figures out if any of the particles are in the wrong position   
void areParticlesFucked() {
    int pf0, pf1, pf2;
    for(int i = 0; i < numberOfParticles; i++) {

        pf0 = particleList[nmod((i-1),numberOfParticles)].currentposition;
        pf1 = particleList[i].currentposition;
        pf2 = particleList[nmod((i+1),numberOfParticles)].currentposition;

        bool comparePiAndMore = pf1 < pf2;
        bool comparePiAndLess = pf1 < pf0;
        bool comparePMoreAndLess = pf2 < pf0;

        if(!(((comparePiAndMore) && (comparePiAndLess) && (comparePMoreAndLess)) || ((comparePiAndMore) && (!comparePiAndLess) && (!comparePMoreAndLess)) || ((!comparePiAndMore) && (!comparePiAndLess) && (comparePMoreAndLess)))) {
            printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH: %d\n", i);
        };
    };
};

//Data analysis methods

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


//number of iterations
int iterationNumber = 20000;


int main(void) {
    printf("nmod: %lf", nmod(-0.99, 10));
    srand(time(0));
    printf("Seed: %d\n",rand());
    loadParticles();
    printCurrentParticleData();
    for(int i = 0; i < iterationNumber; i++) {
        iterate();
        printf("\nIteration number: %d\n", i + 1);
        printCurrentParticleData();
        //areParticlesFucked();
    }
    printf("\nInitial particle data: \n");
    //printCurrentParticleData();
    printInitialParticleData();
    printf("\nMean velocity: %lf", meanVelocity());
    meanResiduals();
};
