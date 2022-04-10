import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

startIndex = 0
endIndex = 0
runDataList = [];
removedCounter = 0
timeCutoff = 1

#loading data into the code
def getDataList(fileName, string):
    with open(fileName) as f:
        fileContents = f.read()
    startIndex = fileContents.find(string) + 2
    startIndex = fileContents.find("{", startIndex) + 1
    endIndex = fileContents.find("}", startIndex) - 1
    strList = fileContents[startIndex:endIndex].split(" ")
    for i in range(0, len(strList)):
        strList[i] = float(strList[i])
    return strList

class runData:
    def __init__(self, fileName):
        self.particleNum = np.array(fileName[fileName.find("P") + 1:fileName.find("L")])
        self.particleNum = self.particleNum.astype(np.float64)
        self.trackLength = np.array(fileName[fileName.find("L") + 1:fileName.find(".")])
        self.trackLength = self.trackLength.astype(np.float64)
        self.times = np.array(getDataList(fileName, "times"))
        self.times = self.times.astype(np.float64)
        self.times = np.sort(self.times)
        self.stdDevs = np.array(getDataList(fileName, "stdDevs"))
        self.stdDevs = self.stdDevs.astype(np.float64)
        self.stdDevs = np.sort(self.stdDevs)
        self.meanVelocities = np.array(getDataList(fileName, "meanVelocities"))
        self.meanVelocities = self.meanVelocities.astype(np.float64)
        self.meanVelocities = np.sort(self.meanVelocities)
        self.leadVelocities = np.array(getDataList(fileName, "leadVelocities"))
        self.leadVelocities = self.leadVelocities.astype(np.float64)
        self.leadVelocities = np.sort(self.leadVelocities)
        self.avgMasses = np.array(getDataList(fileName, "avgMasses"))
        self.avgMasses = self.avgMasses.astype(np.float64)
        self.avgMasses = np.sort(self.avgMasses)
        self.individualCollisionTimes = np.array(getDataList(fileName, "individualCollisionTimes"))
        self.individualCollisionTimes = self.individualCollisionTimes.astype(np.float64)
        self.individualCollisionTimes = np.sort(self.individualCollisionTimes)

#loads files into program
for i in os.listdir("runFiles"):
    f = os.path.join("runFiles", i)
    if os.path.isfile(f):
        runDataList.append(runData(f))
dataListLen = len(runDataList)


#eliminates statistically extreme values
for j in range(len(runDataList)):
    i = 0
    Q1 = runDataList[j].times[int(len(runDataList[j].times) / 4)]
    Q3 = runDataList[j].times[int(len(runDataList[j].times) * (3 / 4))]
    IQR = Q3 - Q1

    while i < len(runDataList[j].times):
        if((runDataList[j].times[i] < (Q1 - (3 * IQR))) or (runDataList[j].times[i] > (Q3 + (3 * IQR)))):
            runDataList[j].times = np.delete(runDataList[j].times, i)
            runDataList[j].stdDevs = np.delete(runDataList[j].stdDevs, i)
            runDataList[j].meanVelocities = np.delete(runDataList[j].meanVelocities, i)
            runDataList[j].leadVelocities = np.delete(runDataList[j].leadVelocities, i)
        else:
            i += 1

    revTimesMax = np.max(runDataList[j].times)
    revTimesMean = np.mean(runDataList[j].times)
    revTimesStd = np.std(runDataList[j].times)
    revTimesMin = np.min(runDataList[j].times)

    print(f'P: {runDataList[j].particleNum} L: {runDataList[j].trackLength} stdDev: {revTimesStd}')
    print(f'Mean of revised1 times: {revTimesMean}')
    print(f'Min of revised1 times: {revTimesMin}')
    print(f'Max of revised1 times: {revTimesMax}')
    print(f'len revTimes: {len(runDataList[j].times)}')

def meanOfEachComponent(list):
    copyList = np.empty(len(list))
    for i in range(len(list)):
        copyList[i] = np.mean(list[i])
    return copyList

particleNumL = np.empty(dataListLen, dtype=object)
trackLengthL = np.empty(dataListLen, dtype=object)
timesL = np.empty(dataListLen, dtype=object)
stdDevsL = np.empty(dataListLen, dtype=object)
meanVelocitiesL = np.empty(dataListLen, dtype=object)
leadVelocitiesL = np.empty(dataListLen, dtype=object)
massesL = np.empty(dataListLen, dtype=object)
individualCollisionTimesL = np.empty(dataListLen, dtype=object)

for i in range(len(runDataList)):
    particleNumL[i] = runDataList[i].particleNum
    trackLengthL[i] = runDataList[i].trackLength
    timesL[i] = runDataList[i].times
    stdDevsL[i] = runDataList[i].stdDevs
    meanVelocitiesL[i] = runDataList[i].meanVelocities
    leadVelocitiesL[i] = runDataList[i].leadVelocities
    massesL[i] = runDataList[i].avgMasses
    individualCollisionTimesL[i] = runDataList[i].individualCollisionTimes

MparticleNumL = meanOfEachComponent(particleNumL)
MtrackLengthL = meanOfEachComponent(trackLengthL)
MtimesL = meanOfEachComponent(timesL)
MstdDevsL = meanOfEachComponent(stdDevsL)
MmeanVelocitiesL = meanOfEachComponent(meanVelocitiesL)
MleadVelocitiesL = meanOfEachComponent(leadVelocitiesL)
MmassesL = meanOfEachComponent(massesL)
MindividualCollisionTimesL = meanOfEachComponent(individualCollisionTimesL)

times = runDataList[8].times
stdDevs = runDataList[8].stdDevs
meanVelocities = runDataList[8].meanVelocities
leadVelocities = runDataList[8].leadVelocities
masses = runDataList[8].avgMasses
individualCollisionTimes = runDataList[8].individualCollisionTimes
logInvCollisionTimes = np.log(np.array(individualCollisionTimes))  
timesS = np.sort(times)
timesMean = np.mean(timesS)
timesStd = np.std(timesS)
timesMax = np.max(timesS)
timesMin = np.min(timesS)

#logTimes = np.log(times)

revTimes = times


#creation of histogram polygon
time_bins = np.linspace(1, revTimesMax, 30)
timeBinCenters = 0.5*(time_bins[1:]+ time_bins[:-1])
timesHistPolygon,edges = np.histogram(revTimes, time_bins)
timesPolygonMaxIndex = timesHistPolygon.tolist().index(max(timesHistPolygon.tolist()))
timesMaxBinLow = time_bins[timesPolygonMaxIndex]

fitTimes = timesHistPolygon[timesPolygonMaxIndex:]
timeFitCenters = timeBinCenters[timesPolygonMaxIndex:]


#creation of bins and binning of times for stddevs and mean velocities
#print(f'stdDevs min: {min(stdDevs)} stdDevs max: {max(stdDevs)}')
numBins = 30

stdDevs_bins = np.linspace(min(stdDevs), max(stdDevs), numBins)
stdDevsRange = np.empty([0])
meanVelocities_bins = np.linspace(min(meanVelocities), max(meanVelocities), numBins)
meanVelocitiesRange = np.empty([0])
stdDevsTimes = np.empty([0])

individualCollisionX = np.arange(0, len(individualCollisionTimes), 1)


for i in range(0, len(stdDevs_bins) - 1):
    sum = 0
    count = 0
    for j in range(0, len(stdDevs)):
            if(stdDevs[j] >= stdDevs_bins[i] and stdDevs[j] < stdDevs_bins[i + 1]):
                sum += revTimes[j]
                count += 1
    if(count == 0):
        stdDevsTimes = np.append(stdDevsTimes, 0)
    else:
        stdDevsTimes = np.append(stdDevsTimes, sum / count)

for i in range(0, len(stdDevs_bins) - 1):
    stdDevsRange = np.append(stdDevsRange, (stdDevs_bins[i] + stdDevs_bins[i + 1]) / 2)

meanVelocitiesTimes = np.empty([0])
for i in range(0, len(meanVelocities_bins) - 1):
    sum = 0
    count = 0
    for j in range(0, len(meanVelocities)):
            if(meanVelocities[j] >= meanVelocities_bins[i] and meanVelocities[j] < meanVelocities_bins[i + 1]):
                sum += revTimes[j]
                count += 1
    if(count == 0):
        meanVelocitiesTimes = np.append(meanVelocitiesTimes, 0)
    else:
        meanVelocitiesTimes = np.append(meanVelocitiesTimes, sum / count)
    
for i in range(0, len(meanVelocities_bins) - 1):
    meanVelocitiesRange = np.append(meanVelocitiesRange, (meanVelocities_bins[i] + meanVelocities_bins[i + 1]) / 2)

meanVelocitiesLead = np.empty([0])
for i in range(0, len(meanVelocities_bins) - 1):
    sum = 0
    count = 0
    for j in range(0, len(meanVelocities)):
            if(meanVelocities[j] >= meanVelocities_bins[i] and meanVelocities[j] < meanVelocities_bins[i + 1]):
                sum += leadVelocities[j]
                count += 1
    if(count == 0):
        meanVelocitiesLead = np.append(meanVelocitiesLead, 0)
    else:
        meanVelocitiesLead = np.append(meanVelocitiesLead, (sum / count) / meanVelocitiesRange[i])

meanVelocitiesLeadBins = np.linspace(0, 1, 30)
meanVelocitiesLeadDistro = np.empty([0])
for i in range(0, len(meanVelocitiesLeadBins) - 1):
    sum = 0
    count = 0
    for j in range(0, len(meanVelocities)):
            if(leadVelocities[j] / meanVelocities[j] >= meanVelocitiesLeadBins[i] and leadVelocities[j] / meanVelocities[j] < meanVelocitiesLeadBins[i + 1]):
                count += 1
    if(count == 0):
        meanVelocitiesLeadDistro = np.append(meanVelocitiesLeadDistro, 0)
    else:
        meanVelocitiesLeadDistro = np.append(meanVelocitiesLeadDistro, count)

meanVelLeadDistroRange = np.empty(0)
for i in range(0, len(meanVelocitiesLeadBins) - 1):
    meanVelLeadDistroRange = np.append(meanVelLeadDistroRange, (meanVelocitiesLeadBins[i] + meanVelocitiesLeadBins[i + 1]) / 2)

#chi squared
def chiSquared(expectedArr, actualArr):
    chiSquaredSum = 0
    for k in range(timesHistPolygon.tolist().index(max(timesHistPolygon.tolist())) + 1, len(expectedArr)):
        ciNom = (actualArr[k] - expectedArr[k])**2
        ciDenom = expectedArr[k]
        chiSquaredSum += abs(ciNom / ciDenom)
    #print(f'chiSquaredSum: {chiSquaredSum} ciNom: {ciNom} ciDenom: {ciDenom}')

    return chiSquaredSum

#==========================================
#--------------BEST FIT LINES--------------
#==========================================

#Weibull equation fit for times
def timeHistWeibullBFL(x, a, b, c, d, e):
    return x
    #return (a * d / b) * np.power((x - e) / (c * b), a - 1) * np.exp(-np.power((x - e) / (c * b), a))

guesses = (1.5, 0.05, 460, 1600, 8)

(a, b, c, d, e), cc = curve_fit(timeHistWeibullBFL, timeBinCenters, timesHistPolygon, p0 = guesses)
(ua, ub, uc, ud, ue) = np.sqrt(np.diag(cc))
(a, b, c, d, e) = (1.5, 0.05, 460, 1600, 8)
xTimesWeibullBFLPlot = np.linspace(timeBinCenters[0], timeBinCenters[len(timeBinCenters) - 1], 100)
timesWeibullBFLPlot = timeHistWeibullBFL(xTimesWeibullBFLPlot, a, b, c, d, e)

print(f'Weibull times run best fit coefficients: a: {a}+-{ua} b: {b}+-{ub} c: {c}+-{uc} d: {d}+-{ud} e: {e}+-{ue}')

#1/x^2 best fit for times
def timeHistBFL(x, a, b, c):
    return (a / ((x - b)**2)) + c

guesses = (1, 1, 1000)

(a, b, c), cc = curve_fit(timeHistBFL, timeFitCenters, fitTimes, p0 = guesses)
(ua, ub, uc) = np.sqrt(np.diag(cc))

xTimesBFLPlot = np.linspace(timeFitCenters[0], timeFitCenters[len(timeFitCenters) - 1], 100)
timesBFLPlot = timeHistBFL(xTimesBFLPlot, a, b, c)

print(f'Run best fit coefficients: a: {a:.3f}+-{ua:.3f} b: {b:.3f}+-{ub:.3f} c: {c:.3f}+-{uc:.3f}')

#1/x^2 best fit for times
def timeHistBFL2(x, a, n):
    return a / (x**2)**n

guesses = (20000000, 1)
weights = np.linspace(1, 10, len(timeFitCenters))

print(timeFitCenters)

(a, n), cc = curve_fit(timeHistBFL2, timeFitCenters, fitTimes, p0 = guesses)
(ua, un) = np.sqrt(np.diag(cc))

xTimesBFLPlot = np.linspace(timeFitCenters[0], timeFitCenters[len(timeFitCenters) - 1], 100)
timesBFLPlot = timeHistBFL2(xTimesBFLPlot, a, n)

print(f'Run best fit coefficients: a: {a} n: {n}')

#time per collision best fit
def timeCollisionBFL(x, a, b, c):
    return (a * np.exp(x / b)) + c

guesses = (1e-17, 4.4, 0)

print(f'individualCollisionX: {individualCollisionX}')
print(f'individualCollisionTimes: {individualCollisionTimes}')

#(a, b, c), cc = curve_fit(timeCollisionBFL, individualCollisionX, individualCollisionTimes, p0 = guesses)
#(ua, ub, u) = np.sqrt(np.diag(cc))

#xinvCollisionBFLPlot = np.linspace(individualCollisionX[0], individualCollisionX[len(individualCollisionX) - 1], 100)
#invCollisionBFLPlot = timeCollisionBFL(xinvCollisionBFLPlot, a, b, c)

print(f'Individual collision times coefficients: a: {a} b: {b} c: {c}')

#common fit functions
def linFitHelper(x, a, b):
    return (a * x) + b

def bestFitLinearPlot(x, y):
    guesses = (1, 1)
    (a, b), cc = curve_fit(linFitHelper, x, y, p0 = guesses)
    return a, b

def logFitHelper(x, a, b, c):
    return (a * np.log(b * x)) + c

def bestFitLogPlot(x, y):
    guesses = (1, 1, 1)
    (a, b, c), cc = curve_fit(logFitHelper, x, y, p0 = guesses)
    return a, b, c

def powerFitHelper(x, a, b, c):
    return (a * (x**b)) + c

def bestFitPowerPlot(x, y):
    guesses = (1, 1, 1)
    (a, b, c), cc = curve_fit(powerFitHelper, x, y, p0 = guesses)
    return a, b, c

revTimesHist = np.histogram(revTimes, bins = time_bins)
revTimesHistMax = max(revTimesHist[0])

#Run length histograms
def runLengthHistogram():
    plt.hist(revTimes, bins = time_bins)
    plt.plot(timeBinCenters, timesHistPolygon,'-*')
    plt.plot(xTimesBFLPlot, timesBFLPlot, "m")
    plt.plot(xTimesWeibullBFLPlot, timesWeibullBFLPlot, color = 'green')
    plt.ylim([0, revTimesHistMax * 1.5])
    plt.xlabel("Run times")
    plt.ylabel("Frequency of run times")
    plt.title("Run time frequencies")
    plt.show()

    plt.plot(xTimesWeibullBFLPlot, timesWeibullBFLPlot, color = 'green')
    plt.show()

def fitRunLengthHist():
    revTimesHist = np.histogram(revTimes, bins = time_bins)
    revTimesHistMax = max(revTimesHist[0])
    plt.hist(revTimes, bins = time_bins)
    plt.plot(timeBinCenters, timesHistPolygon,'-*')
    plt.plot(xTimesBFLPlot, timesBFLPlot, "m")
    plt.ylim([0, revTimesHistMax * 1.5])
    plt.xlabel("Run times")
    plt.ylabel("Frequency of run times")
    plt.title("Run time frequencies")
    plt.show()

def runLengthHistLog():
    plt.hist(revTimes, bins = time_bins)
    f = 1000000 / np.power((time_bins - 10) / 100, 2)
    plt.plot(xTimesBFLPlot, timesBFLPlot, "m")
    #plt.ylim([0, revTimesHistMax * 1.5])
    plt.yscale("log")
    #plt.xscale("log")
    plt.xlabel("Run times")
    plt.ylabel("Frequency of run times")
    plt.title("Run time frequencies")
    plt.show()

def runLengthHistPolygon():
    plt.plot(timeBinCenters, timesHistPolygon,'-*')
    #plt.plot(xTimesBFLPlot, timesBFLPlot, "m")
    plt.ylim([0, revTimesHistMax * 1.5])
    #plt.plot(timeBinCenters,timeBestFit,'-')
    plt.xlabel("Run times")
    plt.ylabel("Frequency of run times")
    plt.title("Polygon of run time frequencies")
    plt.show()

#3D scatterplot showing the track length, number of particles, and average run lengths
def scatterPLRL():
    fig = plt.figure()
    ax = fig.add_subplot(projection = '3d')
    print(f'MtimesL: {MtimesL}')
    ax.scatter(MparticleNumL, MtrackLengthL, MtimesL)
    ax.set_xlabel("Number of particles")
    ax.set_ylabel("Track length")
    ax.set_zlabel("Average run time")
    plt.show()

#scatterplot showing the number of particles vs mean run lengths
def scatterPRL():
    #linear fit
    a, b = bestFitLinearPlot(MparticleNumL, MtimesL)
    print(f'Best linear fit line for number of particles vs mean run times: T_avg = ({a}*x) + {b}')
    scatterLinPRLPlot = linFitHelper(MparticleNumL, a, b)
    plt.plot(MparticleNumL, scatterLinPRLPlot)

    #power fit
    a, b, c = bestFitPowerPlot(MparticleNumL, MtimesL)
    print(f'Best power fit line for number of particles vs mean run times: T_avg = ({a}*x^{b}) + {c}')
    MPNL = MparticleNumL
    scatterPowerPRLPlot = powerFitHelper(np.sort(MPNL), a, b, c)
    plt.plot(np.sort(MPNL), scatterPowerPRLPlot)

    #actual data
    plt.scatter(MparticleNumL, MtimesL)
    plt.xlabel("Number of particles")
    plt.ylabel("Mean time to end of run")
    plt.show()

#3D stddev, mean velocity, and run length scatterplot
def scatterStdMVRL():
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(stdDevs, meanVelocities, revTimes)
    ax.set_xlabel('Standard deviation of velocities')
    ax.set_ylabel('Mean of velocities')
    ax.set_zlabel('Z Label')
    plt.show()

#scatterplot of stddev and run length
def scatterStdRL():
    plt.scatter(stdDevs, revTimes)
    plt.xlabel("Standard deviation of velocities across runs")
    plt.ylabel("Run times")
    plt.title("Correlation between the standard deviation of velocities and run times")
    plt.show()

#plot of stddev and run length
def plotStdRL():
    plt.plot(stdDevsRange, stdDevsTimes)
    plt.xlabel("Standard deviation of velocities across runs")
    plt.ylabel("Run times")
    plt.title("Correlation between the standard deviation of velocities and run times")
    plt.show()

#scatterplot of mean velocities and run length
def scatterMVRL():
    plt.scatter(meanVelocities, revTimes)
    plt.xlabel("Mean of velocities across runs")
    plt.ylabel("Run times")
    plt.title("Correlation between mean velocities and run times")
    plt.show()

#plot of mean velocities vs run length
def plotMVRL():
    plt.plot(meanVelocitiesRange, meanVelocitiesTimes)
    plt.xlabel("Mean of velocities across runs")
    plt.ylabel("Run times")
    plt.title("Correlation between mean velocities and run times")
    plt.show()

#scatterplot of lead velocities vs mean velocities
def scatterLVMV():
    plt.scatter(meanVelocities, leadVelocities)
    plt.xlabel("Mean of velocities across runs")
    plt.ylabel("Lead velocities across runs")
    plt.title("Correlation between lead and mean velocities")
    plt.show()

#line showing average of lead velocities vs mean velocities with respect to the top curve
def plotNormalLVMV():
    ohPtFive = np.linspace(0.5, 0.5, 100)
    xOhPtFive = np.linspace(min(meanVelocitiesRange), max(meanVelocitiesRange), 100)
    plt.plot(meanVelocitiesRange, meanVelocitiesLead)
    plt.plot(xOhPtFive, ohPtFive)
    plt.xlabel("Mean of velocities across runs")
    plt.ylabel("Means of lead velocities over expected slope")
    plt.title("Correlation between lead and mean velocities")
    plt.show()

#distrobution of lead velocities with respect to the mean
def plotNormalLV():
    plt.plot(meanVelLeadDistroRange, meanVelocitiesLeadDistro)
    plt.xlabel("Lead velocities / mean velocities")
    plt.ylabel("Count")
    plt.title("Distribution of lead velocities over mean velocities")
    plt.show()

def plotMassPerCollision():
    mass_X = np.arange(0, len(masses), 1)
    #average mass per collision
    plt.plot(mass_X, masses)
    plt.xlabel("Collision number")
    plt.ylabel("Average mass")
    plt.title("Average mass at given collision number")
    plt.show()

#time between each collision
def plotTimePerCollision():
    plt.plot(individualCollisionX, individualCollisionTimes)
    #plt.plot(xinvCollisionBFLPlot, invCollisionBFLPlot)
    #plt.plot(xinvCollisionBFLPlot, timeCollisionBFL(xinvCollisionBFLPlot, 2.14834340373572e-12, 1.8605148946198493, 2.5667979668948684, 0))
    plt.xlabel("Collision number")
    plt.ylabel("Average time between each collision")
    plt.title("Average time between each collision per collision number")
    plt.show()

    plt.plot(individualCollisionX, individualCollisionTimes)
    #plt.plot(xinvCollisionBFLPlot, invCollisionBFLPlot)
    #plt.plot(xinvCollisionBFLPlot, timeCollisionBFL(xinvCollisionBFLPlot, 2.14834340373572e-12, 1.8605148946198493, 3.5667979668948684, 0))
    plt.yscale("log")
    plt.xlabel("Collision number")
    plt.ylabel("Average time between each collision")
    plt.title("Average time between each collision per collision number")
    plt.show()

plotTimePerCollision()
scatterPLRL()
scatterPRL()

#print(times)