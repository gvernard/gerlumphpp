import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import sys

sampling = 10
emap_width = 760  # two times the emap offset
emap_height = 760 # two times the emap offset

fig = plt.figure(figsize=(20,8))

img = mpimg.imread('data/convolved_map.png')
imgplot = plt.imshow(img)
plt.xlim([0,emap_width])
plt.ylim([0,emap_height])

x1,y1,x2,y2 = np.loadtxt('data/lc_locs.dat',unpack=True)
for i in range(0,len(x1)):
    plt.plot([x1[i]/sampling,x2[i]/sampling],[y1[i]/sampling,y2[i]/sampling],color="black")


plt.show()
