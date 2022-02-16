import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import scipy

count = 0
with open('timeFile.txt') as fp:
    for line in fp:
        if line.strip():
            count += 1

print(count)

times = np.genfromtxt("timeFile.txt", delimiter = " ", usemask = True, max_rows = 1)
times = np.sort(times)
timesMean = np.mean(times)
timesStd = np.std(times)
timesMax = np.max(times)
timesMin = np.min(times)

meanResiduals = times = np.genfromtxt("timeFile.txt", delimiter = " ", usemask = True, skip_header = 1, skip_footer = 0)

masses = np.genfromtxt("massFile.txt", delimiter = " ", usemask = True)

revTimes = np.empty(shape = 0)

for i in times:
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

"""
plt.plot(timeBinCenters,y,'-*')
plt.plot(timeBinCenters,timeBestFit,'-')
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Polygon of run time frequencies")
plt.show()

print(timeBinCenters)
print(timeBestFit)

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