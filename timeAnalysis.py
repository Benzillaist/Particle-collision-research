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

meanResiduals = getDataList("timeFile.txt", "stdDevs")

masses = np.genfromtxt("massFile.txt", delimiter = " ", usemask = True)

revTimes = np.empty(shape = 0)

for i in timesS:
    if i <= (500):
        revTimes = np.append(revTimes, i)

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
time_bins = np.linspace(0, np.amax(revTimesMean + revTimesStd), 15)
timeBinCenters = 0.5*(time_bins[1:]+ time_bins[:-1])
y,edges = np.histogram(revTimes, time_bins)

timeBestFitCoeffs = np.polyfit(timeBinCenters, y, 2)
timeBestFit = timeBestFitCoeffs[2] + timeBestFitCoeffs[1] * pow(timeBinCenters, 1) + timeBestFitCoeffs[0] * pow(timeBinCenters, 2)

plt.hist(revTimes, bins = time_bins)
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Run time frequencies")
plt.show()


plt.plot(timeBinCenters,y,'-*')
#plt.plot(timeBinCenters,timeBestFit,'-')
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Polygon of run time frequencies")
plt.show()

print(meanResiduals)
print(times)

plt.scatter(meanResiduals, times)
plt.show()

"""
plt.plot(timeBinCenters,timeBestFit,'-')
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Best fit of run time frequencies")
plt.show()
"""

mass_X = np.arange(0, masses.size, 1)

plt.plot(mass_X, masses)
plt.xlabel("Collision number")
plt.ylabel("Average mass")
plt.title("Average mass at given collision number")
plt.show()


#print(times)