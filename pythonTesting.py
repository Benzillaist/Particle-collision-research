import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma

ogArr = np.array([1,2,6,8,4,9,6,7,3.5,7.8,1.9])
print(ogArr)

mask1 = (ogArr > 5) 
mask2 = (ogArr >= 8)
mask = mask1
print(mask)

upMask = ma.make_mask(mask)
print(upMask)

print(ogArr[mask])
print(ogArr[upMask])

arrt = np.empty([12])
arrt[:] = 10
print(arrt)


arrp = np.empty([12, 4])
arrp[:, 1] = arrt

print(arrp)