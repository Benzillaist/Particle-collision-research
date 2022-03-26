import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

startIndex = 0
endIndex = 0

def getDataList(fileName, string):
    with open(fileName) as f:
        fileContents = f.read()
    startIndex = fileContents.find(string) + 2
    startIndex = fileContents.find("{", startIndex) + 1
    endIndex = fileContents.find("}", startIndex) - 2    
    strList = fileContents[startIndex:endIndex].split(" ")
    for i in range(0, len(strList)):
        strList[i] = float(strList[i])
    return strList

times = getDataList("timeFile.txt", "times")
timesS = np.sort(times)
timesMean = np.mean(timesS)
timesStd = np.std(timesS)
timesMax = np.max(timesS)
timesMin = np.min(timesS)

print(max(times))

#logTimes = np.log(times)

stdDevs = getDataList("timeFile.txt", "stdDevs")
meanVelocities = getDataList("timeFile.txt", "meanVelocities")
leadVelocities = getDataList("timeFile.txt", "leadVelocities")

masses = np.genfromtxt("massFile.txt", delimiter = " ", usemask = True)

revTimes = times

removedCounter = 0
i = 0
while i < len(revTimes):
    if i >= len(times):
        break
    if (times[i] >= (timesMean + timesStd)):
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

i = 0
while i < len(revTimes):
    if i >= len(revTimes):
        break
    if (revTimes[i] >= (revTimesMean + revTimesStd)):
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
time_bins = np.linspace(1, np.amax(revTimesMax), 30)
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
    return (a * d / b) * np.power(((x / e) - c) / b, a - 1) * np.exp(-np.power(((x / e) - c) / b, a))

guesses = (2, 4.3, 6.4, 100000, 10)

(a, b, c, d, e), cc = curve_fit(timeHistWeibullBFL, timeBinCenters, timesHistPolygon, p0 = guesses)
(ua, ub, uc, ud, ue) = np.sqrt(np.diag(cc))
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

#Run length histograms
revTimesHist = np.histogram(revTimes, bins = time_bins)
revTimesHistMax = max(revTimesHist[0])
plt.plot(xTimesWeibullBFLPlot, timesWeibullBFLPlot, color = 'green')
plt.hist(revTimes, bins = time_bins)
plt.plot(timeBinCenters, timesHistPolygon,'-*')
plt.plot(xTimesBFLPlot, timesBFLPlot, "m")
plt.ylim([0, revTimesHistMax * 1.5])
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Run time frequencies")
plt.show()

plt.plot(xTimesWeibullBFLPlot, timesWeibullBFLPlot, color = 'green')
plt.show()
"""
#Run length fit histograms
revTimesHist = np.histogram(revTimes, bins = time_bins)
revTimesHistMax = max(revTimesHist[0])

plt.hist(revTimes, bins = time_bins)
plt.plot(timeBinCenters, timesHistPolygon,'-*')
#plt.plot(xTimesBFLPlot, timesBFLPlot, "m")
plt.ylim([0, revTimesHistMax * 1.5])
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Run time frequencies")
plt.show()
"""
#Run length histograms
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

#Polygon of run length histogram
plt.plot(timeBinCenters, timesHistPolygon,'-*')
#plt.plot(xTimesBFLPlot, timesBFLPlot, "m")
plt.ylim([0, revTimesHistMax * 1.5])
#plt.plot(timeBinCenters,timeBestFit,'-')
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Polygon of run time frequencies")
plt.show()

#3D stddev, mean velocity, and run length scatterplot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(stdDevs, meanVelocities, revTimes)
ax.set_xlabel('Standard deviation of velocities')
ax.set_ylabel('Mean of velocities')
ax.set_zlabel('Z Label')
plt.show()

#scatterplot of stddev and run length
plt.scatter(stdDevs, revTimes)
plt.xlabel("Standard deviation of velocities across runs")
plt.ylabel("Run times")
plt.title("Correlation between the standard deviation of velocities and run times")
plt.show()

#plot of stddev and run length
plt.plot(stdDevsRange, stdDevsTimes)
plt.xlabel("Standard deviation of velocities across runs")
plt.ylabel("Run times")
plt.title("Correlation between the standard deviation of velocities and run times")
plt.show()

#scatterplot of mean velocities and run length
plt.scatter(meanVelocities, revTimes)
plt.xlabel("Mean of velocities across runs")
plt.ylabel("Run times")
plt.title("Correlation between mean velocities and run times")
plt.show()

#plot of mean velocities vs run length
plt.plot(meanVelocitiesRange, meanVelocitiesTimes)
plt.xlabel("Mean of velocities across runs")
plt.ylabel("Run times")
plt.title("Correlation between mean velocities and run times")
plt.show()

#scatterplot of lead velocities vs mean velocities
plt.scatter(meanVelocities, leadVelocities)
plt.xlabel("Mean of velocities across runs")
plt.ylabel("Lead velocities across runs")
plt.title("Correlation between lead and mean velocities")
plt.show()

#line showing average of lead velocities vs mean velocities with respect to the top curve
ohPtFive = np.linspace(0.5, 0.5, 100)
xOhPtFive = np.linspace(min(meanVelocitiesRange), max(meanVelocitiesRange), 100)
plt.plot(meanVelocitiesRange, meanVelocitiesLead)
plt.plot(xOhPtFive, ohPtFive)
plt.xlabel("Mean of velocities across runs")
plt.ylabel("Means of lead velocities over expected slope")
plt.title("Correlation between lead and mean velocities")
plt.show()

#distrobution of lead velocities with respect to the mean
plt.plot(meanVelLeadDistroRange, meanVelocitiesLeadDistro)
plt.xlabel("Lead velocities / mean velocities")
plt.ylabel("Count")
plt.title("Distrobution of lead velocities over mean velocities")
plt.show()

#best fit time line
"""
plt.plot(timeBinCenters,timeBestFit,'-')
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Best fit of run time frequencies")
plt.show()
"""

mass_X = np.arange(0, masses.size, 1)

#average mass per collition
plt.plot(mass_X, masses)
plt.xlabel("Collision number")
plt.ylabel("Average mass")
plt.title("Average mass at given collision number")
plt.show()


#print(times)