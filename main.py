import math
import random


#GLOBAL VARIABLES
global stdAcceleration
global numberOfParticles
global minVelocity
global maxVelocity
global trackLength
global particleList

std_acceleration = 1
numberOfParticles = 5
minVelocity = 5
maxVelocity = 10
trackLength = 1000


particle = {
    "idealVelocity": 0

}

#Standard particle class for storing particle data
class Particle:
    def __init__(self, idealvelocity, initialposition):
        self.idealvelocity = idealvelocity
        self.initialposition = initialposition
        self.currentposition
        self.currentvelocity

    def get_IdealVelocity(self):
        return self.idealvelocity

    def get_InitialPosition(self):
        return self.initialposition

    def get_CurrentPosition(self):
        return self.currentposition

    def get_currentvelocity(self):
        return self.currentvelocity

    def set_CurrentPosition(self, position):
        currentposition = position

    def set_CurrentVelocity(self, velocity):
        currentvelocity = velocity;

#manages the velocities of the particle after the collision
def Collision(particle1, particle2):
    E_neti = (particle1.get_currentvelocity())^2 + (particle2.get_currentvelocity())^2
    p_neti = particle1.get_currentvelocity() + particle2.get_currentvelocity()

    v_1 = 1/2 * (math.sqrt((2 * E_neti) - p_neti**2) + p_neti)
    v_2 = 1/2 * (p_neti - math.sqrt((2 * E_neti) - p_neti**2))
    if(v_1 != particle1.get_currentvelocity()):
        particle1.set_CurrentVelocity(v_1)
        particle2.set_CurrentVelocity(v_2)
    else:
        particle1.set_CurrentVelocity(v_2)
        particle2.set_CurrentVelocity(v_1)

#finds the initial spacing between particles
def findSpacing():
    return trackLength/numberOfParticles

#returns a single partcile object
def generateParticle(idealvelocity, initialposition):
    particle = Particle(idealvelocity, initialposition)
    return particle

#generates the particles and adds them to the particle list
def loadParticles():
    for i in range(numberOfParticles):
        particleList[i] = generateParticle(randVelocity(), i*findSpacing())

def randVelocity():
    return round(random.uniform(minVelocity, maxVelocity),6)

def print_hi(name):
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.
    for i in range(numberOfParticles):
        print(i)
        print(randVelocity())

if __name__ == '__main__':
    print_hi('PyCharm')

