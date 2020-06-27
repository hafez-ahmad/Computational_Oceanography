# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 16:50:07 2020

@author: hafez
"""


#https://www.codeproject.com/Articles/1087025/Using-Python-to-Solve-Computational-Physics-Proble

import numpy as np
import matplotlib.pyplot as plt
# Set maximum iteration
maxIter = 500

# Set Dimension and delta
lenX = lenY = 20 #we set it rectangular
delta = 1

# Boundary condition
Ttop = 0 #100,0
Tbottom = 100
Tleft = 30
Tright = 30 #0,30

# Initial guess of interior grid
Tguess = 30
# Set colour interpolation and colour map.
# You can try set it to 10, or 100 to see the difference
# You can also try: colourMap = plt.cm.coolwarm
colorinterpolation = 50
colourMap = plt.cm.jet

#set mesh

X,Y=np.meshgrid(np.arange(0,lenX),np.arange(0,lenY))

# Set array size and set the interior value with Tguess
T = np.empty((lenX, lenY))
T.fill(Tguess)

# Set Boundary condition
T[(lenY-1):, :] = Ttop # bottom row 100
T[:1, :] = Tbottom  #top row 0
T[:, (lenX-1):] = Tright
T[:, :1] = Tleft

# Iteration (We assume that the iteration is convergence in maxIter = 500)
print("Please wait for a moment")
for k in range(0,maxIter):
    for i in range(1,lenX-1,delta):
        for j in range(1,lenY-1,delta):
            T[i,j]=0.25*(T[i+1][j]+T[i-1][j]+T[i][j+1]+T[i][j-1])
        
print('finisted')
#Configure the contour
plt.title("Contour of Temperature")
plt.contourf(X,Y,T, colorinterpolation, cmap=colourMap)

# Set Colorbar
plt.colorbar()

# Show the result in the plot window
plt.show()

#http://kitchingroup.cheme.cmu.edu/pycse/pycse.html

import numpy as np
from scipy.integrate import odeint

Ca0 = 2     # Entering concentration
vo = 2      # volumetric flow rate
volume = 20 # total volume of reactor, spacetime = 10
k = 1       # reaction rate constant

N = 100     # number of points to discretize the reactor volume on

init = np.zeros(N)    # Concentration in reactor at t = 0
init[0] = Ca0         # concentration at entrance

V = np.linspace(0, volume, N) # discretized volume elements
tspan = np.linspace(0, 25)    # time span to integrate over

def method_of_lines(C, t):
    'coupled ODES at each node point'
    D = -vo * np.diff(C) / np.diff(V) - k * C[1:]**2
    return np.concatenate([[0],D])

sol = odeint(method_of_lines, init, tspan)

    
# steady state solution
def pfr(C, V):
    return 1.0 / vo * (-k * C**2)

ssol = odeint(pfr, Ca0, V)

plt.plot(tspan,sol[:,1])
plt.show()
plt.plot(V, ssol, label='Steady state')
plt.plot(V, sol[-1], label='t = {}'.format(tspan[-1]))
plt.xlabel('Volume')
plt.ylabel('$C_A$')
plt.legend(loc='best')


from matplotlib import animation 

fig=plt.figure()
ax=plt.axes()
line,=ax.plot(V,init,lw=2)
plt.show()

def anim(i):
    line.set_xdata(V)
    line.set_ydata(sol[i])
    ax.set_title('t = {0}'.format(tspan[i]))
    ax.figure.canvas.draw()
    return line,

anim=animation.FuncAnimation(fig,anim,frames=50,blit=True)
anim
plt.show()

anim.save('transient_pfr.mp4', fps=10)
