import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt

times = np.genfromtxt("timeFile.txt", delimiter = " ", usemask = True)
times = np.sort(times)
timesMean = np.mean(times)
timesStd = np.std(times)

revTimes = np.empty(shape = 0)

for i in times:
    if i <= (timesMean + (timesStd)):
        revTimes = np.append(revTimes, i) 

#print(revTimes)
time_bins = np.linspace(0, np.amax(revTimes), 20)

plt.hist(revTimes, bins = time_bins)
plt.xlabel("Run times")
plt.ylabel("Frequency of run times")
plt.show()

print(f'Standard deviation of times: {timesStd}')
print(f'Mean of times: {timesMean}')

#print(times)