import math
import random


#GLOBAL VARIABLES
global stdAcceleration
global numberOfParticles
global minVelocity
global maxVelocity
global trackLength

std_acceleration = 1
numberOfParticles = 10
minVelocity = 5
maxVelocity = 10
trackLength = 10
tickRate = 1000
particleRadius = 0.01
particleList = []
tickFrequency = 1/tickRate

#Standard particle class for storing particle data
class Particle:
    def __init__(self, idealvelocity, initialposition):
        self.idealvelocity = idealvelocity
        self.initialposition = initialposition
        self.currentvelocity = 0
        self.currentposition = 0

#finds the initial spacing between particles
def findSpacing():
    return trackLength/numberOfParticles

#returns a random velocity between the min and max velocities specified above
def randVelocity():
    return round(random.uniform(minVelocity, maxVelocity),6)

#returns a single particle object
def generateParticle(idealvelocity, initialposition):
    particle = Particle(idealvelocity, initialposition)
    return particle

#generates the particles and adds them to the particle list
def loadParticles():
    for i in range(numberOfParticles):
        particleList.append(generateParticle(randVelocity(), i*findSpacing()))
        particleList[i].currentposition = particleList[i].initialposition
        particleList[i].currentvelocity = particleList[i].idealvelocity

#manages the velocities of the particle after the collision
def Collision(particle1, particle2):
    pv1 = particle1.currentvelocity
    pv2 = particle2.currentvelocity
    particle1.currentvelocity = pv2
    particle2.currentvelocity = pv1

def didCollide(particle1, particle2):
    pi1 = particle1.currentposition
    pi2 = particle2.currentposition
    if(pi1>pi2):
        pi2 += trackLength
    
    #estimates where the particle will be after one tick and checks if it steps into the bounding box of another particle
    pf1 = (pi1 + (particle1.currentvelocity*tickFrequency) + particleRadius)%trackLength
    pf2 = (pi2 + (particle2.currentvelocity*tickFrequency) - particleRadius)%trackLength
    if(pf1>=pf2):
        return True
    else:
        return False

def accelerateParticles():
    for i in range(numberOfParticles):
        particleList[i].currentposition = (particleList[i].currentposition + (particleList[i].currentvelocity*tickFrequency))%10
        if(abs(particleList[i].idealvelocity - particleList[i].currentvelocity) <= (std_acceleration*tickFrequency)):
            particleList[i].currentvelocity = particleList[i].idealvelocity
        elif((particleList[i].idealvelocity - particleList[i].currentvelocity) < 0):
            particleList[i].currentvelocity -= (std_acceleration*tickFrequency)
        elif((particleList[i].idealvelocity - particleList[i].currentvelocity) > 0):
            particleList[i].currentvelocity += (std_acceleration*tickFrequency)

#advances the particles one tick forwards
def iterate():
    for i in range(numberOfParticles):
        if(didCollide(particleList[i], particleList[(i+1)%numberOfParticles])):
            Collision(particleList[i], particleList[(i+1)%numberOfParticles])
    accelerateParticles()

#figures out if any of the particles are in the wrong position   
def areParticlesFucked():
    for i in range(numberOfParticles):
        if(not (((particleList[i].currentposition < particleList[(i+1)%numberOfParticles].currentposition) and (particleList[i].currentposition < particleList[(i-1)%numberOfParticles].currentposition) and (particleList[(i+1)%numberOfParticles].currentposition < particleList[(i-1)%numberOfParticles].currentposition)) or ((particleList[i].currentposition < particleList[(i+1)%numberOfParticles].currentposition) and (particleList[i].currentposition > particleList[(i-1)%numberOfParticles].currentposition) and (particleList[(i+1)%numberOfParticles].currentposition > particleList[(i-1)%numberOfParticles].currentposition)) or ((particleList[i].currentposition > particleList[(i+1)%numberOfParticles].currentposition) and (particleList[i].currentposition > particleList[(i-1)%numberOfParticles].currentposition) and (particleList[(i+1)%numberOfParticles].currentposition < particleList[(i-1)%numberOfParticles].currentposition)))):
            print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAH: " + str(i))

#Data analysis methods

#finds the mean velocity of all of the particles
def meanVelocity():
    v_sum = 0
    for i in range(numberOfParticles):
        v_sum += particleList[i].currentvelocity
    return (v_sum / numberOfParticles)

#finds the difference from the ideal velocity of each particle to the mean velocity of all of the particles current velocities
def meanResiduals():
    meanvelocity = meanVelocity()
    lowestResidual = abs(particleList[0].idealvelocity - meanvelocity)
    lowestIndex = 0
    print("\nMean residuals:")
    for i in range(numberOfParticles):
        individualResidual = abs(particleList[i].idealvelocity - meanvelocity)
        if (individualResidual < lowestResidual):
            lowestResidual = individualResidual
            lowestIndex = i
        print(str(i) + ": " + str(individualResidual))
    print("Lowest residual: " + str(lowestIndex) + ": " + str(lowestResidual))

#prints the current stats on all of the particles
def printCurrentParticleData():
    for i in range(numberOfParticles):
        print(str(i) + ": Current position: " + str(particleList[i].currentposition) + ", Current velocity: " + str(particleList[i].currentvelocity))

#prints the intial stats of all the particles
def printInitialParticleData():
    for i in range(numberOfParticles):
        print(str(i) + ": Initial position: " + str(particleList[i].initialposition) + ", Ideal velocity: " + str(particleList[i].idealvelocity))

def print_hi(name):
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
    for i in range(numberOfParticles):
        print(i)
        print(randVelocity())

#number of iterations
iterationNumber = 100000

if __name__ == '__main__':
    #print_hi('PyCharm')
    loadParticles()
    printCurrentParticleData()
    for i in range(iterationNumber):
        iterate()
        print("\nIteration number: " + str(i + 1))
        printCurrentParticleData()
        areParticlesFucked()
    print("\nInitial particle data: ")
    printInitialParticleData()
    print("Mean velocity: " + str(meanVelocity()))
    meanResiduals()