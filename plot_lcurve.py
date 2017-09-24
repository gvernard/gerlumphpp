import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import sys
import os.path


index = int(sys.argv[1])
sampling = 10
emap_width = 760  # two times the emap offset
emap_height = 760 # two times the emap offset

fig = plt.figure(figsize=(20,8))



x1,y1,x2,y2 = np.loadtxt('data/lc_locs.dat',unpack=True)





xmin = min(x1[index],x2[index])/sampling
xmax = max(x1[index],x2[index])/sampling
ymin = min(y1[index],y2[index])/sampling
ymax = max(y1[index],y2[index])/sampling
Dx = xmax - xmin
Dy = ymax - ymin
D = max(Dx,Dy)
dD = 10

ax1 = fig.add_subplot(121)
ax1.set_xlim([xmin-dD,xmin+D+dD])
ax1.set_ylim([ymin-dD,ymin+D+dD])
img = mpimg.imread('data/convolved_map.png')
imgplot = ax1.imshow(img)

plt.plot([x1[index]/sampling,x2[index]/sampling],[y1[index]/sampling,y2[index]/sampling])




ax2 = fig.add_subplot(122)
ax2.set_ylim([0,10])


filename = 'data/lcdata_full_'+str(index)+'.dat'
if os.path.isfile(filename):
    x,y,e = np.loadtxt(filename,unpack=True)
    ax2.fill_between(x,y-e,y+e,alpha=0.3,lw=0,color="black")
    ax2.plot(x,y,color="black")


filename = 'data/lcdata_sampled_'+str(index)+'.dat'
if os.path.isfile(filename):
    x,y,e = np.loadtxt(filename,unpack=True)
    ax2.errorbar(x,y,yerr=e,color="blue")


filename = 'data/lcdata_strategy_'+str(index)+'.dat'
if os.path.isfile(filename):
    x,y,e = np.loadtxt(filename,unpack=True)
    ax2.errorbar(x,y,yerr=e,color="red")






plt.show()
