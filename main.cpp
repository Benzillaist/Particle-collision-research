#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

double std_acceleration = 1;
int numberOfParticles = 10;
double minVelocity = 5;
double maxVelocity = 10;
int trackLength = 10;
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

//finds the initial spacing between particles
double findSpacing() {
    return trackLength/numberOfParticles;
};

//returns a random velocity between the min and max velocities specified above
double randVelocity() {
    int randToInt = (1000000 * (rand() / 32767));
    double randToDouble = minVelocity + ((maxVelocity - minVelocity) * (((double)randToInt) / 1000000));
    return randToDouble;
};

//returns a single particle object
Particle generateParticle(double idealvelocity, double initialposition) {
    Particle particle = Particle(idealvelocity, initialposition);
    return particle;
};

//generates the particles and adds them to the particle list
void loadParticles() {
    for(int i = 0; i < numberOfParticles; i++) {
        particleList.push_back(generateParticle(randVelocity(), i*findSpacing()));
        particleList[i].currentposition = particleList[i].initialposition;
        particleList[i].currentvelocity = particleList[i].idealvelocity;
    }
}

//manages the velocities of the particle after the collision
void Collision(Particle particle1, Particle particle2) {
    double pv1 = particle1.currentvelocity;
    double pv2 = particle2.currentvelocity;
    particle1.currentvelocity = pv2;
    particle2.currentvelocity = pv1;
};

bool didCollide(Particle particle1, Particle particle2) {
    double pi1 = particle1.currentposition;
    double pi2 = particle2.currentposition;
    if(pi1>pi2) {
        pi2 += trackLength;
    };
    
    //estimates where the particle will be after one tick and checks if it steps into the bounding box of another particle
    double pf1 = fmod(pf1 + (double)(particle1.currentvelocity*tickFrequency + particleRadius), trackLength);
    double pf2 = fmod(pf2 + (double)(particle2.currentvelocity*tickFrequency + particleRadius), trackLength);
    
    if(pf1>=pf2) {
        return 1;
    }
    else {
        return 0;
    };
}

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
        if(didCollide(particleList[i], particleList[(i+1)%numberOfParticles])) {
            Collision(particleList[i], particleList[(i+1)%numberOfParticles]);
        }
    };
    accelerateParticles();
}

//figures out if any of the particles are in the wrong position   
void areParticlesFucked() {
    for(int i = 0; i < numberOfParticles; i++) {
        if(!(((particleList[i].currentposition < particleList[(i+1)%numberOfParticles].currentposition) && (particleList[i].currentposition < particleList[(i-1)%numberOfParticles].currentposition) && (particleList[(i+1)%numberOfParticles].currentposition < particleList[(i-1)%numberOfParticles].currentposition)) || ((particleList[i].currentposition < particleList[(i+1)%numberOfParticles].currentposition) && (particleList[i].currentposition > particleList[(i-1)%numberOfParticles].currentposition) && (particleList[(i+1)%numberOfParticles].currentposition > particleList[(i-1)%numberOfParticles].currentposition)) || ((particleList[i].currentposition > particleList[(i+1)%numberOfParticles].currentposition) && (particleList[i].currentposition > particleList[(i-1)%numberOfParticles].currentposition) && (particleList[(i+1)%numberOfParticles].currentposition < particleList[(i-1)%numberOfParticles].currentposition)))) {
            printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH: %d", i);
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
    printf("\nMean residuals:");
    for(int i = 0; i < numberOfParticles; i++) {
        individualResidual = abs(particleList[i].idealvelocity - meanvelocity);
        if (individualResidual < lowestResidual) {
            lowestResidual = individualResidual;
            lowestIndex = i;
        };
        printf("%d: %d", i, individualResidual);
    };
    printf("Lowest residual: %d: %d", lowestIndex, lowestResidual);
}

//prints the current stats on all of the particles
void printCurrentParticleData() {
    for(int i = 0; i < numberOfParticles; i++) {
        printf("%d: Current position: %d, Current velocity: %d", i, particleList[i].currentposition, particleList[i].currentvelocity);
    };
};

//prints the intial stats of all the particles
void printInitialParticleData() {
    for(int i = 0; i < numberOfParticles; i++) {
        printf("%d: Initial position: %d, Ideal velocity: %d", i, particleList[i].initialposition, particleList[i].idealvelocity);
    };
};


//number of iterations
int iterationNumber = 100000;


int main(void) {
    srand(time(0));
    printf("Seed: %d\n",rand());
    loadParticles();
    printCurrentParticleData();
    for(int i = 0; i < numberOfParticles; i++) {
        iterate();
        printf("\nIteration number: %d", i + 1);
        printCurrentParticleData();
        areParticlesFucked();
    }
    printf("\nInitial particle data: ");
    printInitialParticleData();
    printf("Mean velocity: %d", meanVelocity());
    meanResiduals();
};