import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy

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
y,edges = np.histogram(revTimes, time_bins)

#creation of best-fit line
timeBestFitCoeffs = np.polyfit(timeBinCenters, y, 2)
timeBestFit = timeBestFitCoeffs[2] + timeBestFitCoeffs[1] * pow(timeBinCenters, 1) + timeBestFitCoeffs[0] * pow(timeBinCenters, 2)

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
    

#finding best fit lines
def chiSquared(expectedArr, actualArr):
    chiSquaredSum = 0
    for k in range(y.tolist().index(max(y.tolist())) + 1, len(expectedArr)):
        ciNom = (actualArr[k] - expectedArr[k])**2
        ciDenom = expectedArr[k]
        chiSquaredSum += abs(ciNom / ciDenom)
    #print(f'chiSquaredSum: {chiSquaredSum} ciNom: {ciNom} ciDenom: {ciDenom}')

    return chiSquaredSum

x = timeBinCenters
runLengthBFL = 0
iB = -1
jB = -1
kB = -1
chiSquaredBest = -1
for i in range(800,900, 10):
    for j in range(-50, 1, 1):
        for k in range(300, 351, 1):
            runLengthBFL = i / np.power((timeBinCenters + j) / k, 2)
            chiSquaredTemp = chiSquared(runLengthBFL, y)
            if(chiSquaredBest < 0 or chiSquaredTemp < chiSquaredBest):
                print(chiSquaredTemp)
                chiSquaredBest = chiSquaredTemp
                iB = i
                jB = j
                kB = k
runLengthBFL = iB / np.power((timeBinCenters + jB) / kB, 2)
print(runLengthBFL)
print(f'iB: {iB} jB: {jB} kB: {kB}')

print(f'arr max index: {y.tolist().index(max(y.tolist()))}')
print(y)

#Run length histograms
revTimesHist = np.histogram(revTimes, bins = time_bins)
revTimesHistMax = max(revTimesHist[0])

plt.hist(revTimes, bins = time_bins)
plt.plot(timeBinCenters,y,'-*')
plt.plot(x, runLengthBFL, "m")
plt.ylim([0, revTimesHistMax * 1.5])
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Run time frequencies")
plt.show()

#Run length histograms
plt.hist(revTimes, bins = time_bins)
plt.plot(x, runLengthBFL)
f = 1000000 / np.power((time_bins - 10) / 100, 2)
plt.plot(x, runLengthBFL, "m")
plt.ylim([0, revTimesHistMax * 1.5])
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Run time frequencies")
plt.show()

#Polygon of run length histogram
plt.plot(timeBinCenters,y,'-*')
plt.plot(x, runLengthBFL, "m")
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