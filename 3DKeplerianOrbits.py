
# coding: utf-8

# In[199]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
#%matplotlib osx    
import orbital
import math


# In[172]:


def Newtons_Equation(t, X):
    r = X[0:3]
    v = X[3:6]

    
    F = -G * ms * mp /np.linalg.norm(r)**3 * r
    a    = F / mp
    
    return np.concatenate((v, a))


# In[225]:


ms = 3e+6
mp = 1
M = ms * mp
G  = 1e-6
mu = G * M

r_moon = 385000 #km # 385,000 / 6,371 (earth radius) = 60 times earth's radius
v_moon = 1.022 #km/s #3679.2 #km/h 

# //// Fienga et al 2016 ///
i_p9 = 30 #degrees
Omega_p9 = 113 #degrees
argp_p9 = 150 #degrees
nu_p9 = 118 #degrees


# In[226]:


x = 1 
y = 1
z = 1
vx = 0
vy = 0.5
vz = 0

r = np.array([x,y,z])
v = np.array([vx,vy,vz])

print('r = ', r)
print('v = ', v)

#/// Calculate orbital elements /// 
#follows from the method laid out in Fundamentals of Astrodynamics and Applications, by Vallado, 2007

h=np.cross(r,v)         # angular momentum
K = np.array([0, 0, 1])
nhat=np.cross(K,h)       #node vector

#eccentricity vector
evec = ((np.linalg.norm(v)**2 - mu/np.linalg.norm(r))*r - np.dot(r,v)*v)/ mu  
e = np.linalg.norm(evec)

if e == 0:
    orbit = 'circular orbit'
if e > 0 and e < 1:
    orbit = 'elliptical orbit'
if e == 1:
    orbit ='parabolic orbit'
if e > 1:
    orbit = 'hyperbolic orbit'
print('e = ', e, orbit)

energy = np.linalg.norm(v)**2/2 - mu/np.linalg.norm(r) #Specific mechanical energy
print('E = ', energy)

if abs(e) != 1: 
   a = -mu/(2*energy)
   #p = a*(1-e**2)
else:
   #p = np.linalg.norm(h)**2/mu
   a = 'inf'
    
print('a = ', a)         #semi major axis

p = np.sqrt(abs(a)**3)
print('p = ', p)        #period

i = np.arccos(h[2]/np.linalg.norm(h))  #inclination
print('i = ', i)        

Omega = np.arccos(nhat[0]/np.linalg.norm(nhat))  #swivel: the angle from the principal direction to the ascending node, 0° ≤ Ω < 360°
if nhat[1] < 0:
   Omega = 360-Omega
print('Omega = ', Omega)

argp = np.arccos(np.dot(nhat,evec)/(np.linalg.norm(nhat)*e)) # argument of perigee, ω: 0° ≤ ω < 360
if evec[2]<0:
   argp = 360-argp
print('argp = ', argp)

nu_0 = np.arccos(np.dot(evec,r)/(e*np.linalg.norm(r)))  # True anomaly, ν, is the angle along the orbital path from perigee to the spacecraft’s position vector r.  0 <= v < 360°
if np.dot(r,v)<0:
   nu = 360 - nu_0
print('initial nu = ', nu_0)  # changes with time, location of the spacecraft in its orbit


# In[233]:


#Calculate r using orbital elements
U, V, W = orbital.utilities.uvw_from_elements(math.radians(i_p9), raan=math.radians(Omega_p9), arg_pe=math.radians(argp_p9), f=nu_p9)
x = np.linalg.norm(U)
y = np.linalg.norm(V)
z = np.linalg.norm(W)
r = np.array([x,y,z])
print(r)


# In[228]:


dt = .001
time = np.arange(0,10.001, dt)
tmax = 10


# In[229]:


xt = []
yt = []
zt = []
vx = []
vy = []
vz = []

position = []

tt = []

t = 0


X = np.concatenate((r,v))

get_ipython().run_line_magic('matplotlib', 'osx')


while True:

    r = X[0:3]
    v = X[3:6]
    
    
    xt.append(X[0])
    yt.append(X[1])
    zt.append(X[2])
    
    
    #X = X + Newtons_Equation(X) * h    #EULERS METHOD
    
    f1 = Newtons_Equation(t       ,X          )    #RUNGE KUTTA METHOD
    f2 = Newtons_Equation(t+dt/2.0,X+f1*dt/2.0)
    f3 = Newtons_Equation(t+dt/2.0,X+f2*dt/2.0)
    f4 = Newtons_Equation(t+dt    ,X+f3*dt    )

    X = X + (f1 + 2.0*f2 + 2.0*f3 + f4) / 6.0 * dt
    
    nu = np.arccos(np.dot(evec,r)/(e*np.linalg.norm(r)))  # True anomaly, ν, is the angle along the orbital path from perigee to the spacecraft’s position vector r.  0 <= v < 360°
    if np.dot(r,v)<0:
        nu = 360 - nu
    #print('nu = ', nu)
    if nu *10000 //1 / 10000 == (nu_0-.001)*10000 //1 / 10000:
        break
    
    tt.append(t)

    if t >= tmax:
        break
    t += dt


# In[166]:


print(t)
print(p)
print(p/dt)


# In[155]:


#time plot
plt.plot(tt,yt, 'r-')
plt.title('Height of Y over time')
plt.show()


# In[152]:


#time plot
plt.plot(tt,zt, 'r-')
plt.title('Height of Z over time')
plt.show()


# In[ ]:





# In[168]:


#2D orbit plot
plt.plot(xt, yt, 'r-')
plt.xlabel('X axis')
plt.ylabel('Y axis')
plt.plot(0,0,'*',mfc='w',ms=10)
plt.gca().set_aspect('equal', adjustable='box')  
plt.show()


# In[169]:


#2D orbit plot
plt.plot(xt, zt, 'r-')
plt.xlabel('X axis')
plt.ylabel('Z axis')
plt.plot(0,0,'*',mfc='w',ms=10)
plt.gca().set_aspect('equal', adjustable='box')  
plt.show()


# In[14]:


#time plot
plt.plot(tt,zt, 'r-')
plt.title('Height of Z over time')
plt.show()


# In[230]:


#3D orbit plot
XT = np.array(xt)
YT = np.array(yt)
ZT = np.array(zt)


plt.rcParams['legend.fontsize'] = 10

fig = plt.figure(figsize=(5,5))
ax = fig.gca(projection='3d')
ax.set_aspect("equal")


ax.plot(XT, YT, ZT, label='Orbit', antialiased=False)     #orbit
#ax.text(0, 0, 0, '*', zdir=None)


u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:100j]   #sphere
x = 0.1 * np.cos(u)*np.sin(v)
y = 0.1 * np.sin(u)*np.sin(v)
z = 0.1 * np.cos(v)
ax.plot_surface(x, y, z, color="b")


ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-1.5, 1.5)
ax.set_zlim(-1.5, 1.5)
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

plt.draw()
plt.pause(0.005)
plt.show()

