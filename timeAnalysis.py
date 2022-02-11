import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

times = np.genfromtxt("timeFile.txt", delimiter = " ", usemask = True)
times = np.sort(times)
timesMean = np.mean(times)
timesStd = np.std(times)

masses = np.genfromtxt("timeFile.txt", delimiter = " ", usemask = True)

revTimes = np.empty(shape = 0)

for i in times:
    if i <= (timesMean + (timesStd)):
        revTimes = np.append(revTimes, i) 

#print(revTimes)
time_bins = np.linspace(0, np.amax(revTimes), 20)

plt.hist(revTimes, bins = time_bins)
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.title("Run time frequencies")
plt.show()

mass_bins = np.arange(0, masses.size, 1)

plt.hist(masses, bins = mass_bins)
plt.xlabel("Collision number")
plt.ylabel("Average mass")
plt.title("Average mass at given collision number")
plt.show()

print(f'Standard deviation of times: {timesStd}')
print(f'Mean of times: {timesMean}')

#print(times)