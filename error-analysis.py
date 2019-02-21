########################################################################
# Team Rabbits: Austin Czyzewski | Megan Davis | Molly Janasik | David Pham
# AST 304, Fall 2018
# Michigan State University
########################################################################

"""
<Computes the error in total energy from Euler Integration along with 2nd and 4th order Runge-Kutta.>
"""
#imports
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import matplotlib
from ode import fEuler, rk2, rk4
from kepler import *

#######################################################################
# Setting up our problem
#######################################################################

m = 1
a = 1
e = 0.5
z0, eps0, Tperiod = set_initial_conditions(a, m, e) #set initial conditions

real_TE = total_energy(z0,m) #Store total energy for use later

#######################################################################
# Setting up the different time steps to work with
#Store both step sizes and names
#these stepsizes are getting progressively smaller and smaller
#######################################################################

h0 = 0.1*Tperiod
h = []
h_names = []
i = 0
while 2**i < 1025: #2,4,8,16,32,64,128,256,512,1024
  h.append(h0/2**i)
  h_names.append("$h_0$ / %s" % (2**i))
  i += 1
#print(h)

#######################################################################
# Testing each integration for Total Energy
#######################################################################

Total_Energy = [] #Store last TE value for each int method
for method in integration_methods:
  for step in h:
    ts, Xs, Ys, KEs, PEs, TEs = integrate_orbit(z0=z0,m=m,tend=Tperiod*3,h=step,method=method)
    Total_Energy.append(TEs[-1]) # This will give us a list of all errors, each should fill up a third of this list
    
#######################################################################
# Plotting
#######################################################################

Euler = abs(Total_Energy[0:11] - real_TE) #Retrieve relevant data for Euler
RK2 = abs(Total_Energy[11:22] - real_TE)
RK4 = abs(Total_Energy[22:33] - real_TE)

plt.figure(figsize = (12,8)) #Make figures
ax = plt.subplot()
plt.plot(h,Euler,'r',label = "Euler") #plot int method TE error as function of H
plt.plot(h,RK2,'g',label = "RK2")
plt.plot(h,RK4,'b',label = "RK4")
plt.yscale('log')
plt.xscale('log')


#Sorry, old stuff we don't want to delete
#for i in range(len(Total_Energy)):
#  if i <= len(Total_Energy)/3:
#    plt.plot(i,Total_Energy[i],'r') #Plotting the Euler method
#  elif len(Total_Energy)/3 < i < 2*len(Total_Energy)/3:
#    plt.plot(i,Total_Energy[i],'g') #Plotting the RK2 method
#  else:
#    plt.plot(i,Total_Energy[i],'b') #Plotting the RK4 method


#Format Plots
plt.title("Total Energy Error as a function of time-step size")
plt.xlabel("Time-Step size in terms of $h_0$")
plt.ylabel("Error")
plt.legend()
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter()) #adjust xaxis
plt.xticks(fontsize = 10)
plt.xticks(h,h_names) # A label for each of the ticks for the appropriate h value
plt.savefig("Error_analysis.jpg", overwrite=True) #Set to false if you don't want to overwrite
plt.show() #Make sure savefig before show, otherwise the figure is blank
