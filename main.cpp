#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

//declaration of tweakable variables, keep tickrate very high and don't add too many particles or increase the radius by too much
double std_acceleration = 1;
int numParticles = 100;
int particleListSize = numParticles;
double minVelocity = 0;
double maxVelocity = 10;
double trackLength = 10;
double tickRate = 10000;
double particleRadius = 0.02;
double tickFrequency = 1/tickRate;

//max number of iterations for each run, the total seconds that elapse is equal to iterationNumber/tickRate
int iterationNumber = 5000000;

//collision type, for elastic put true, for inelastic put false
bool elasticCollision = false;

//number of runs that will be ran
int numberOfRuns = 100000;

//Lists:

//Standard particle class for storing particle data
class Particle {
    public:

    //each particle consists of the ideal velocity and starting position, which will stay static, and then the current position and velocity, which will be changed during the simulation
        double idealvelocity;
        double initialposition;
        double currentvelocity;
        double currentposition;
        int currentmass;
        int leadParticle;
        double leadParticleVelocity;
        Particle(double idealv, double initialp, int ID) { //constructor function
            idealvelocity = idealv;
            initialposition = initialp;
            currentmass = 1;
            leadParticle = ID;
            leadParticleVelocity = idealv;
        }
};

vector<Particle> particleList;
vector<double> averageMassList (numParticles, 0);
vector<double> finalMassList (numParticles, 0);
vector<double> deltaTimeList;
vector<double> tempTimeList;
vector<double> stdDevList;
vector<double> meanVelocityList;
vector<double> leadVelocityList;
double tempTime = 0;

//METHODS

//integer modulo function that handles negatives >:(
double nmod(double dividend, double divisor) {
    if(dividend < 0) {
        return (dividend + divisor);
    }
    return fmod(dividend, divisor);
};

//finds the initial spacing between particles
double findSpacing() {
    return trackLength/particleListSize;
};

//returns a random velocity between the min and max velocities specified above
double randVelocity() {
    double cutRand = (double)rand() / 32767;
    return (double)minVelocity + (((double)maxVelocity - (double)minVelocity) * cutRand);
};

//returns a single particle object
Particle generateParticle(double idealvelocity, double initialposition) {
    Particle particle = Particle(idealvelocity, initialposition, particleList.size());
    return particle;
};

//generates the particles and adds them to the particle list
void loadParticles() {
    for(int i = 0; i < particleListSize; i++) {
        double rV = randVelocity();
        particleList.push_back(generateParticle(rV, i*findSpacing()));
        particleList[i].currentposition = particleList[i].initialposition;
        particleList[i].currentvelocity = particleList[i].idealvelocity;
    }
};


//outdated collision management method
void E_didCollide(int i) {
    
    //estimates where the particle will be after one tick and checks if it steps into the bounding box of another particle
    double pf0 = fmod(particleList[nmod((i-1),particleListSize)].currentposition + (double)(particleList[nmod((i-1),particleListSize)].currentvelocity*tickFrequency + particleRadius), trackLength);
    double pf1l = fmod(particleList[i].currentposition + (double)(particleList[i].currentvelocity*tickFrequency - particleRadius), trackLength);
    double pf1m = fmod(particleList[i].currentposition + (double)(particleList[i].currentvelocity*tickFrequency + particleRadius), trackLength);
    double pf2 = fmod(particleList[nmod((i+1),particleListSize)].currentposition + (double)(particleList[nmod((i+1),particleListSize)].currentvelocity*tickFrequency - particleRadius), trackLength);
    
    //compares positions
    bool comparePiAndMore = pf1m < pf2;
    bool comparePiAndLess = pf1l < pf0;
    bool comparePMoreAndLess = pf2 < pf0;

    //checks if compared positions are correct
    if(!(((comparePiAndMore) && (comparePiAndLess) && (comparePMoreAndLess)) || ((comparePiAndMore) && (!comparePiAndLess) && (!comparePMoreAndLess)) || ((!comparePiAndMore) && (!comparePiAndLess) && (comparePMoreAndLess)))) {
        //printf("yup: %d\n", i);
        double pv1 = particleList[i].currentvelocity;
        double pv2 = particleList[nmod((i+1),particleListSize)].currentvelocity;
        particleList[i].currentvelocity = pv2;
        particleList[nmod((i+1),particleListSize)].currentvelocity = pv1;
        //didCollide(nmod((i-1),particleListSize));
        //didCollide(nmod((i+1),particleListSize));
    };
};

//updated collision management function
void E_didCollideV2() {
    double x0, x1, distance, v0, v1;
    for(int i = 0; i < particleListSize; i++) { //for each particle and the succeeding one:

        //get both particle's positions
        x0 = particleList[i].currentposition;
        x1 = particleList[nmod(i+1, particleListSize)].currentposition;

        //finds the distance between
        distance = x1 - x0; 

        //get each particles velocity
        v0 = particleList[i].currentvelocity;
        v1 = particleList[nmod(i+1, particleListSize)].currentvelocity;

        if(distance <= (-2 * particleRadius)) {
            distance += trackLength; //if the distance is negative, it means that the succeeding particle has looped around, so i add on to the distance to account for that change
        }
        if((abs(distance) < 2 * particleRadius) && (v0 > v1)) { //if the particles are overlapping and are coming closer together, the velocities switch
            //swaps velocities of the two particles
            particleList[i].currentvelocity = v1;
            particleList[nmod((i+1),particleListSize)].currentvelocity = v0;
        }
    }
}



//function that accelerates particles back to their ideal velocity if the particle is not at previously states ideal velocity
void E_accelerateParticles() {
    for(int i = 0; i < particleListSize; i++) {

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

//processes the next inelastic collision
void I_processNextCollision() {
    double minCollisionTime = -1;
    int minCollisionIndex, p0Mass, p1Mass, totalMass;
    double deltaVelocity, deltaPosition, deltaTime;

    for(int i = 0; i < particleListSize; i++)
    {
        deltaVelocity = particleList[i].currentvelocity-particleList[nmod(i+1,particleListSize)].currentvelocity;
        deltaPosition = particleList[nmod(i+1,particleListSize)].currentposition-particleList[i].currentposition;

        if(deltaPosition<0) {
            deltaPosition += trackLength;
        }

        deltaTime = deltaPosition / deltaVelocity;
        //printf("\nIndex: %d tempminCollisionTime: %lf", i, deltaTime);
        //printf("\ndeltaVelocity: %lf deltaPosition: %lf", deltaVelocity, deltaPosition);
        if(deltaVelocity>0) { 
            if(((deltaTime)<minCollisionTime)||(minCollisionTime<0)) {
                minCollisionIndex = i;
                minCollisionTime = deltaTime;
            }
        }
    }

    //printf("\nminCollisionIndex: %d minCollisionTime: %lf", minCollisionIndex, minCollisionTime);

    for(int i = 0; i < particleListSize; i++) {
        particleList[i].currentposition = nmod(particleList[i].currentposition + (particleList[i].idealvelocity * minCollisionTime), trackLength);
    }

    if(minCollisionTime == -1) {
        printf("no elegible particle was removed\n");
    }
    else {
        tempTime += minCollisionTime;
        p0Mass = particleList[minCollisionIndex].currentmass;
        p1Mass = particleList[nmod(minCollisionIndex + 1, particleListSize)].currentmass;
        totalMass = p0Mass + p1Mass;
        particleList[minCollisionIndex].currentvelocity = ((particleList[minCollisionIndex].currentvelocity * (double)p0Mass) + (particleList[nmod(minCollisionIndex + 1, particleListSize)].currentvelocity * (double)p1Mass))/ (double)totalMass;
        particleList[minCollisionIndex].currentmass = totalMass;
        //printf("\nvelocity %d: %lf velocity %d: %lf", minCollisionIndex, particleList[minCollisionIndex].currentvelocity, (int)nmod(minCollisionIndex + 1, particleListSize), particleList[nmod(minCollisionIndex + 1, particleListSize)].currentvelocity);
        //printf("\nleadVelocity %d: %lf leadVelocity %d: %lf", minCollisionIndex, particleList[minCollisionIndex].leadParticleVelocity, (int)nmod(minCollisionIndex + 1, particleListSize), particleList[nmod(minCollisionIndex + 1, particleListSize)].leadParticleVelocity);
        particleList[minCollisionIndex].leadParticle = particleList[nmod(minCollisionIndex + 1, particleListSize)].leadParticle;
        particleList[minCollisionIndex].leadParticleVelocity = particleList[nmod(minCollisionIndex + 1, particleListSize)].leadParticleVelocity;
        particleList.erase(particleList.begin() + nmod(minCollisionIndex + 1, particleListSize));
        //printf("Particle number %d was removed\n", particleList.begin() + nmod(minCollisionIndex + 1, particleListSize));
    }

}



//figures out if any of the particles are in the wrong position   
void areParticlesMessedUp() {
    int pf0, pf1, pf2;
    for(int i = 0; i < particleListSize; i++) {

        //gets position of each particle, the one in front of it, and the one behind it
        pf0 = particleList[nmod((i-1),particleListSize)].currentposition;
        pf1 = particleList[i].currentposition;
        pf2 = particleList[nmod((i+1),particleListSize)].currentposition;

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
int lowestInitVIndex() {
    double lV = particleList[0].idealvelocity;
    double tV;
    int li = 0;
    for(int i = 1; i < particleListSize; i++) {
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
    int i = lowestInitVIndex(); //first index
    while(((particleList[nmod(i+1, particleListSize)].currentposition-particleList[i].currentposition) < (radiusTolerance * particleRadius)) && ((particleList[nmod(i+1, particleListSize)].currentposition-particleList[i].currentposition)>=0)) {
        i = nmod(i+1, particleListSize);
    }
    return i;
};

//tests if every particle is in the same chain
bool fullHouse() {
    int i = nmod(findFront()+1, particleListSize); //first index
    int nP = 0; //number of particles in loop
    while(((particleList[nmod(i+1, particleListSize)].currentposition-particleList[i].currentposition) < (radiusTolerance * particleRadius)) && ((particleList[nmod(i+1, particleListSize)].currentposition-particleList[i].currentposition)>=0)) {
        i = nmod(i+1, particleListSize);
        nP++;
    }
    if(nP == particleListSize-1) {
        return 1;
    }
    return 0;
}

//finds the mean velocity of all of the particles
double E_meanVelocity() {
    double v_sum = 0;

    for(int i = 0; i < particleListSize; i++) {
        v_sum += particleList[i].currentvelocity;
    };

    return (v_sum / particleListSize);
};

//finds the mean velocity of all of the particles
double I_meanVelocity() {
    double v_sum = 0;
    int mass_sum = 0;
    int tempmass = 0;

    for(int i = 0; i < particleListSize-1; i++) {
        tempmass = particleList[i].currentmass;
        v_sum += (particleList[i].currentvelocity * tempmass);
        mass_sum += tempmass;
    };

    return (v_sum / mass_sum);
};

//finds the difference from the ideal velocity of each particle to the mean velocity of all of the particles current velocities
void meanResiduals() {
    double meanvelocity = E_meanVelocity();
    double lowestResidual = abs(particleList[0].idealvelocity - meanvelocity);
    int lowestIndex;
    double individualResidual;
    double residualSum = 0;
    vector<double> tempResidualList;
    printf("\nMean residuals:\n");
    for(int i = 0; i < particleListSize; i++) {
        tempResidualList.push_back(particleList[i].idealvelocity - meanvelocity);
        individualResidual = abs(particleList[i].idealvelocity - meanvelocity);
        if (individualResidual < lowestResidual) {
            lowestResidual = individualResidual;
            lowestIndex = i;
        };
        printf("%d: %lf\n", i, individualResidual);
    };
    for(int i = 0; i < particleListSize; i++) {
        residualSum += tempResidualList[i] * tempResidualList[i];
        printf("\nTemp residual: %lf", tempResidualList[i]);
    }
    stdDevList.push_back(sqrt(residualSum/(double)tempResidualList.size()));
    meanVelocityList.push_back(meanvelocity);
    printf("Lowest residual: %d: %lf\n", lowestIndex, lowestResidual);
}

//prints the current stats on all of the particles
void E_printCurrentParticleData() {
    for(int i = 0; i < particleListSize; i++) {
        printf("%d: Current position: %lf, Current velocity: %lf\n", i, particleList[i].currentposition, particleList[i].currentvelocity);
    };
};

//prints the current stats on all of the particles
void I_printCurrentParticleData() {
    for(int i = 0; i < particleListSize-1; i++) {
        printf("%d: Current position: %lf, Current velocity: %lf, Current mass: %d\n", i, particleList[i].currentposition, particleList[i].currentvelocity, particleList[i].currentmass);
    };
};

//prints the intial stats of all the particles
void printInitialParticleData() {
    for(int i = 0; i < particleListSize; i++) {
        printf("%d: Initial position: %lf, Ideal velocity: %lf\n", i, particleList[i].initialposition, particleList[i].idealvelocity);
    };
};


//returns the largest mass out of all the particles
int E_largestMass() {
    double lM = particleList[0].currentmass;
    double tM;
    int li = 0;
    for(int i = 1; i < particleListSize; i++) {
        tM = particleList[i].currentmass;
        if(tM>lM) {
            lM = tM;
        }
    }
    return lM;
};

//returns average mass of all the particles
double E_averageMass() {
    double sum = 0;
    for(int i = 0; i < particleListSize; i++) {
        sum += (double)particleList[i].currentmass;
    }
    return sum/(double)particleListSize;
}

//prints the average masses from all the runs
void E_printRunEndData() {
    printf("\nAverage largest mass after each collision:");
    for(int i = 0; i < averageMassList.size(); i++) {
        printf("\n%d: %lf", i, averageMassList[i] / double(numberOfRuns));
        finalMassList[i] = averageMassList[i] / double(numberOfRuns);
    }

    printf("\nLead particle ID: %d", particleList[0].leadParticle);
    printf("\nLead particle velocity: %lf", particleList[0].leadParticleVelocity);

    leadVelocityList.push_back(particleList[0].leadParticleVelocity);
}

//prints lead particle data
void leadParticleVelocity() {

}


//Main method, runs code :D

int main(void) {
    int isFrontRun = 0;
    int completedRuns = 0;

    //runs the number of previously specified simulations
    for(int j = 0; j < numberOfRuns; j++) {
        //printf("J: %d\n", j);

        //resets the particle list from any past runs
        particleList.clear();

        //resets the random number generator "seed" so I don't get the same particles each time
        srand(time(0)+j);

        //prints the "seed", so each run can be re-ran at a later time
        printf("\nSeed: %d\n",rand());

        //generates all the particles based off the previous random number generation
        loadParticles();

        if(elasticCollision == true) {

            //prints initial particle data
            E_printCurrentParticleData();

            //runs for a maximum of cycles that was previously defined
            for(int i = 0; i <= iterationNumber; i++) {

                //manages and collisions that would have happened for each particle
                E_didCollideV2();

                //moves and accelerates each particle
                E_accelerateParticles();

                //checks every 1000 ticks to see if the particles have all bunched up, if so, it prints the particle data and then leaves the loop
                if((i%1000) == 0) {

                    if(fullHouse() == 1) {
                        printf("\nIteration number: %d\n", i);
                        E_printCurrentParticleData();
                        completedRuns++;
                        break;    
                    }
                }

                //since printing out each tick would clog the console and not be useful for debugging or analyzing, I only print out every 10000 ticks
                if((i%10000) == 0) {
                    printf("\nIteration number: %d\n", i);
                    E_printCurrentParticleData();
                }
            }
            
            //general methods that are useful for analyzing data and debugging
            printf("\nInitial particle data: \n");
            printInitialParticleData();
            //printf("\nMean velocity: %lf", E_meanVelocity());
            meanResiduals();
            printf("Lowest ideal velocity: %d: %lf\n", lowestInitVIndex(), particleList[lowestInitVIndex()].idealvelocity);
            printf("Start of chain: %d\n", findFront());

            //if the particle with the lowest velocity is the one at the front of the chain, it adds one to the counter
            if((findFront() == lowestInitVIndex()) && (fullHouse() == 1)) {
                isFrontRun++;
            }
        }
        else if(elasticCollision == false) {

            //averageMassList.assign(10,0);
            tempTime = 0;

            printf("\nInitial particle data: \n");
            printInitialParticleData();
            printf("\nMean velocity: %lf", I_meanVelocity());
            meanResiduals();
            printf("Lowest ideal velocity: %d: %lf\n", lowestInitVIndex(), particleList[lowestInitVIndex()].idealvelocity);
            printf("Start of chain: %d\n", findFront());

            while(particleListSize>1) {
                averageMassList[(int)(numParticles-particleListSize)] += E_largestMass();
                //cout << "\n" << particleListSize << "\n";
                I_processNextCollision();
                //I_printCurrentParticleData();
                printf("\nMean velocity: %lf\n", I_meanVelocity());
                particleListSize--;
            }

            deltaTimeList.push_back(tempTime);
            averageMassList[(int)(numParticles-particleListSize)] += E_largestMass();
            E_printRunEndData();

            particleListSize = numParticles;
        }
    }

    //prints the number of successful runs to the number of runs where the particle with the lowest velocity was at the front
    double sum = 0;
    if(elasticCollision == false) {
        for(int i = 0; i < deltaTimeList.size(); i++) {
            sum += deltaTimeList[i];
        }

        ofstream timeFile ("timeFile.txt");

        timeFile << "times{";
        for(int i = 0; i < deltaTimeList.size(); i ++){
            timeFile << deltaTimeList[i] << " " ;
        }
        timeFile << "}stdDevs{";
        for(int i = 0; i < stdDevList.size(); i ++) {
            timeFile <<  stdDevList[i] << " ";
        }
        timeFile << "}meanVelocities{";
        for(int i = 0; i < meanVelocityList.size(); i ++) {
            timeFile <<  meanVelocityList[i] << " ";
        }
        timeFile << "}leadVelocities{";
        for(int i = 0; i < leadVelocityList.size(); i ++) {
            timeFile << leadVelocityList[i] << " ";
        }
        timeFile << "}";

        ofstream massFile ("massFile.txt");

        for(int i = 0; i < finalMassList.size(); i ++){
            massFile << finalMassList[i] << " " ;
        }
    }

    if(elasticCollision == true) {
        printf("Total completed runs: %d, of which, %d had the particle with the lowest ideal velocity at the front\n", completedRuns, isFrontRun);
    }
};