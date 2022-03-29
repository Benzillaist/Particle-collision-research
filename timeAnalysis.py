import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

startIndex = 0
endIndex = 0
runDataList = [];

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
        self.particleNum = fileName[fileName.find("P") + 1:fileName.find("L")]
        self.trackLength = fileName[fileName.find("L") + 1:fileName.find(".")]
        self.times = getDataList(fileName, "times")
        self.stdDevs = getDataList(fileName, "stdDevs")
        self.meanVelocities = getDataList(fileName, "meanVelocities")
        self.leadVelocities = getDataList(fileName, "leadVelocities")
        self.avgMasses = getDataList(fileName, "avgMasses")
        self.individualCollisionTimes = getDataList(fileName, "individualCollisionTimes")

for i in os.listdir("runFiles"):
    f = os.path.join("runFiles", i)
    if os.path.isfile(f):
        runDataList.append(runData(f))

times = runDataList[8].times
stdDevs = runDataList[8].stdDevs
meanVelocities = runDataList[8].meanVelocities
leadVelocities = runDataList[8].leadVelocities
masses = runDataList[8].avgMasses
individualCollisionTimes = runDataList[8].individualCollisionTimes
print(masses[-1])
timesS = np.sort(times)
timesMean = np.mean(timesS)
timesStd = np.std(timesS)
timesMax = np.max(timesS)
timesMin = np.min(timesS)

print(max(times))

#logTimes = np.log(times)

revTimes = times

removedCounter = 0
i = 0
while i < len(revTimes):
    if i >= len(times):
        break
    if (times[i] >= ((timesMean + timesStd) / 10)):
        revTimes.pop(i)
        stdDevs.pop(i)
        meanVelocities.pop(i)
        leadVelocities.pop(i)
        removedCounter += 1
    else:
        i += 1

revTimesMax = np.max(revTimes)
revTimesMean = np.mean(revTimes)
revTimesStd = np.std(revTimes)
revTimesMin = np.min(revTimes)

print(f'Standard deviation of revised times: {revTimesStd}')
print(f'Mean of revised times: {revTimesMean}')
print(f'Min of revised times: {revTimesMin}')
print(f'Max of revised times: {revTimesMax}')

i = 0
while i < len(revTimes):
    if i >= len(revTimes):
        break
    if (revTimes[i] >= ((revTimesMean + revTimesStd) / 10)):
        revTimes.pop(i)
        stdDevs.pop(i)
        meanVelocities.pop(i)
        leadVelocities.pop(i)
        removedCounter += 1
    else:
        i += 1



revTimesMean = np.mean(revTimes)
revTimesStd = np.std(revTimes)
revTimesMax = np.max(revTimes)
revTimesMin = np.min(revTimes)

print(f'Standard deviation of times: {timesStd}')
print(f'Mean of times: {timesMean}')
print(f'Min of times: {timesMin}')
print(f'Max of times: {timesMax}')

print(f'Standard deviation of revised times: {revTimesStd}')
print(f'Mean of revised times: {revTimesMean}')
print(f'Min of revised times: {revTimesMin}')
print(f'Max of revised times: {revTimesMax}')
#print(times)

print(f'Removed Counter: {removedCounter}')

#creation of histogram polygon
time_bins = np.linspace(1, revTimesMax, 30)
print(time_bins[0])
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
    return (a * d / b) * np.power((x - e) / (c * b), a - 1) * np.exp(-np.power((x - e) / (c * b), a))

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
def timeHistBFL(x, a, n):
    return a / (x**2)**n

guesses = (20000000, 1)
weights = np.linspace(1, 10, len(timeFitCenters))

(a, n), cc = curve_fit(timeHistBFL, timeFitCenters, fitTimes, p0 = guesses)
(ua, un) = np.sqrt(np.diag(cc))

xTimesBFLPlot = np.linspace(timeFitCenters[0], timeFitCenters[len(timeFitCenters) - 1], 100)
timesBFLPlot = timeHistBFL(xTimesBFLPlot, a, n)

print(f'Run best fit coefficients: a: {a} n: {n}')

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
    individualCollisionX = np.arange(0, len(individualCollisionTimes), 1)
    plt.plot(individualCollisionX, individualCollisionTimes)
    plt.yscale("log")
    plt.xlabel("Collision number")
    plt.ylabel("Average time between each collision")
    plt.title("Average time between each collision per collision number")
    plt.show()

plotTimePerCollision()

#print(times)