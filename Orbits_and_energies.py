########################################################################
# Team Rabbits: Austin Czyzewski | Megan Davis | Molly Janasik | David Pham
# AST 304, Fall 2018
# Michigan State University
########################################################################

"""
This code will create orbit plots and energy plots.

Please note: plots will be saved in a folder unless you comment out the plt.savefig code noted.

"""
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
from ode import fEuler, rk2, rk4
from kepler import *
from matplotlib import rcParams



print("Hi, the output (plots) of this script will be saved. Just a warning!")
#######################################################################
# Setting up our problem
#######################################################################

m = 1
a = 1
e = -m/(2*a)
z0, eps0, Tperiod = set_initial_conditions(a, m, e)

######################################################################
#Step_size creation
######################################################################
h0 = 0.1*Tperiod
step_sizes = []
h_names = []
i = 0
while 2**i < 1025:
  step_sizes.append(h0/2**i)
  h_names.append("$h_0$ / %s" % (2**i))
  i += 1

def orbit_plotting(tend1= 3*Tperiod, h1=h0, h2=h0/2, method1='Euler'):
    """
    A function to plot the orbits with 2 different step_sizes with any given integration method.
    
    Arguments
    
        tend1: endtime, default: 3*Tperiod
    
        h1: stepsize, default: h0
    
        h2: stepsize, default: h0/2
    
        method1: integration method, default: Euler
    
    
    Returns: no defined return statement. Will plot of one method with two different stepsizes
    
    """
    ts1, Xs1, Ys1, KEs1, PEs1, TEs1= integrate_orbit(z0=z0,m=m,tend=tend1,h=h1,method= method1)
    ts2, Xs2, Ys2, KEs2, PEs2, TEs2= integrate_orbit(z0=z0,m=m,tend=tend1,h=h2,method= method1)
    
    ax1 = plt.subplot(211)
    ax1.plot(Xs1, Ys1)
    ax1.set_title('Method:{} Stepsize:{}'.format(method1, h_names[0]))
    ax1.set_xlabel('X Position')
    ax1.set_ylabel('Y Position')
    
    ax2 = plt.subplot(212)
    ax2.plot(Xs2, Ys2)
    ax2.set_title('Method:{} Stepsize:{}'.format(method1, h_names[-1]))
    ax2.set_xlabel('X Position')
    ax2.set_ylabel('Y Position')
    
    plt.tight_layout()
    plt.savefig("Two_orbits_{}.jpg".format(method1), overwrite=True) #comment out if you don't want to save them
    plt.show()
   
def energy_plotting(tend1= 3*Tperiod, h1=h0, h2=h0/2, method1='Euler'):
    """
    A function to plot the kinetic, potential, and total energies
    with 2 different step_sizes with any given integration method.
    
    Arguments
    
        tend1: endtime, default: 3*Tperiod
    
        h1: stepsize, default: h0
    
        h2: stepsize, default: h0/2
    
        method1: integration method, default: Euler
    
    
    Returns: No defined return statement. Will plot of one method with two different stepsizes
    
    """
    ts1, Xs1, Ys1, KEs1, PEs1, TEs1= integrate_orbit(z0=z0,m=m,tend=tend1,h=h1,method= method1)
    ts2, Xs2, Ys2, KEs2, PEs2, TEs2= integrate_orbit(z0=z0,m=m,tend=tend1,h=h2,method= method1)
    
    
    
    fig = plt.figure(figsize = (8,6))
    ax = plt.subplot(111)

    ax.plot(ts1, KEs1,'r--' , label='KE, Step: '+h_names[0])
    ax.plot(ts1, PEs1,'b--', label='PE, Step: '+h_names[0])
    ax.plot(ts1, TEs1,'g--', label='TE, Step: '+h_names[0])
    ax.plot(ts2, KEs2,'r', label='KE, Step: '+h_names[-1])
    ax.plot(ts2, PEs2,'b', label='PE, Step: '+h_names[-1])
    ax.plot(ts2, TEs2, 'g', label='TE, Step: '+h_names[-1])

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_title('Kinetic, Potential, and Total Energies for {}'.format(method1))
    ax.set_xlabel('Time')
    ax.set_ylabel('Energy')
    plt.savefig("Energies_{}.jpg".format(method1), pad_inches=1, overwrite=True) #comment out if you don't want to save them
    plt.show()



 #Plot Orbits and Energies according to Method   
orbit_plotting(h2=step_sizes[-1])
orbit_plotting(h2=step_sizes[-1], method1='RK2')
orbit_plotting(h2=step_sizes[-1], method1='RK4')
    
energy_plotting(h2=step_sizes[-1])
energy_plotting(h2=step_sizes[-1], method1='RK2')
energy_plotting(h2=step_sizes[-1], method1='RK4')






