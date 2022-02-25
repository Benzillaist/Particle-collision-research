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

#logTimes = np.log(times)

stdDevs = getDataList("timeFile.txt", "stdDevs")

meanVelocities = getDataList("timeFile.txt", "meanVelocities")

leadVelocities = getDataList("timeFile.txt", "leadVelocities")

masses = np.genfromtxt("massFile.txt", delimiter = " ", usemask = True)

revTimes = times

removedCounter = 0
for i in range(0, len(times)):
    if i >= len(times):
        break
    if not(times[i] <= (timesMean + timesStd)):
        revTimes.pop(i)
        stdDevs.pop(i)
        meanVelocities.pop(i)
        leadVelocities.pop(i)
        removedCounter += 1
        

revTimesMean = np.mean(revTimes)
revTimesStd = np.std(revTimes)
revTimesMax = np.max(revTimes)
revTimesMin = np.min(revTimes)

print(f'Standard deviation of times: {timesStd}')
print(f'Mean of times: {timesMean}')
print(f'Min of times: {timesMin}')

print(f'Standard deviation of revised times: {revTimesStd}')
print(f'Mean of revised times: {revTimesMean}')
print(f'Min of revised times: {revTimesMin}')
#print(times)

print(f'Removed Counter: {removedCounter}')

#creation of histogram polygon
time_bins = np.linspace(0, np.amax(revTimesMean + timesStd), 15)
timeBinCenters = 0.5*(time_bins[1:]+ time_bins[:-1])
y,edges = np.histogram(revTimes, time_bins)

#creation of best-fit line
timeBestFitCoeffs = np.polyfit(timeBinCenters, y, 2)
timeBestFit = timeBestFitCoeffs[2] + timeBestFitCoeffs[1] * pow(timeBinCenters, 1) + timeBestFitCoeffs[0] * pow(timeBinCenters, 2)



#Run length histograms
plt.hist(revTimes, bins = time_bins)
plt.yscale("log")
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Run time frequencies")
plt.show()
"""
#log of run length histograms
logTime_bins = np.arange(0, max(times))
logTimes = np.histogram(times)
plt.hist(logTimes, bins = logTime_bins)
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Run time frequencies")
plt.show()
"""
#Polygon of run length histogram
plt.plot(timeBinCenters,y,'-*')
#plt.plot(timeBinCenters,timeBestFit,'-')
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Polygon of run time frequencies")
plt.show()

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
                sum += times[j]
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
                sum += times[j]
                count += 1
    if(count == 0):
        meanVelocitiesTimes = np.append(meanVelocitiesTimes, 0)
    else:
        meanVelocitiesTimes = np.append(meanVelocitiesTimes, sum / count)
    
for i in range(0, len(meanVelocities_bins) - 1):
    meanVelocitiesRange = np.append(meanVelocitiesRange, (meanVelocities_bins[i] + meanVelocities_bins[i + 1]) / 2)
    

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

#scatterplot of mean velocities and run length
plt.scatter(meanVelocities, revTimes)
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

#plot of stddev and run length
plt.plot(stdDevsRange, stdDevsTimes)
plt.show()

#plot of mean velocities vs run length
plt.plot(meanVelocitiesRange, meanVelocitiesTimes)
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