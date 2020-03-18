import matplotlib.animation as animation
import numpy as np
import h5py
from pylab import *

def ani_frame():
    tight_layout()
    plt.axis('off')


def update_img(key):
    dataset=f[key]
    im.set_data(dataset)
    return im

    #legend(loc=0)
f = h5py.File('data.h5', 'r');
keys=list(f.keys())

fig = plt.figure()
ax = fig.add_subplot(111)
im = ax.imshow(rand(1024,1024),interpolation='nearest')
ani = animation.FuncAnimation(fig,update_img,keys,init_func=ani_frame,interval=30)
writer = animation.writers['ffmpeg'](fps=30)

ani.save('demo.mkv',writer=writer)
