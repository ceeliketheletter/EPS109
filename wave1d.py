
# coding: utf-8

# In[23]:


## %matplotlib inline

get_ipython().run_line_magic('matplotlib', 'osx')


# In[24]:


import matplotlib.pyplot as plt
import numpy as np


# In[25]:


from matplotlib.animation import FFMpegWriter
metadata = dict(title='My first animation in 2D', artist='Matplotlib',comment='Wakanda is coming.')
writer = FFMpegWriter(fps=15, metadata=metadata)
fig = plt.figure()

with writer.saving(fig, "wave1d.mp4", dpi=200):
    nf = 100
    for it in range(nf):
        if (it%10==0): print(it,end='')
        print('.',end='')

        n = 50
        y = np.zeros(n)
        f = 2.0*np.pi/n
        for i in range(n):
            y[i] = np.cos(f*(i+it)) + np.sin(f*it)*np.cos(3*f*(i+it))
        plt.clf()
        plt.plot(y, 'ro-',mfc='w')
        plt.show()
        plt.draw()
        plt.pause(0.05)
        writer.grab_frame()

