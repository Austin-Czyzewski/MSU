########################################################################
# Team Rabbits: Austin Czyzewski | Megan Davis | Molly Janasik | David Pham
# AST 304, Fall 2018
# Michigan State University
########################################################################

"""
<Integration of orbits and computation of energies.>
"""

import numpy as np
from numpy.linalg import norm
# ode.py is the routine containing integrators
from ode import fEuler, rk2, rk4

# Use this to identify methods
integration_methods = {
    'Euler':fEuler,
    'RK2':rk2,
    'RK4':rk4
    }
    
# Energies
def kinetic_energy(v):
    """
    Returns kinetic energy per unit mass: KE(v) = 0.5 v*v.
    
    Arguments
        v (array-like)
            velocity vector
    """
    return 0.5*np.dot(v,v)

def potential_energy(x,m):
    """
    Returns potential energy per unit mass: PE(x, m) = -m/norm(r)

    Arguments
        x (array-like)
            position vector
        m (scalar)
            total mass in normalized units
    """
    r = norm(x)
    return -m/r

def total_energy(z,m):
    """
    Returns energy per unit mass: E(z,m) = KE(v) + PE(x,m)

    Arguments
        z (array-like)
            position and velocity vectors
            position = z[0:2]
            velocity = z[2:4]
        m (scalar)
            total mass in normalized units
    """
    # In order to break z into position, velocity vectors, we use array slices:
    # Here, z[n:m] means taking the elements of z with indices n <= j < m
    r = z[0:2]  # Start with index 0 and take two indices: 0 and 1
    v = z[2:4]  # Start with index 2 and take two indices: 2 and 3
    KE = kinetic_energy(v)
    PE = potential_energy(r,m)
    total_E = KE + PE

    # Replace the following two lines
    return total_E

def derivs(t,z,m):
    """
    Computes derivatives of position and velocity for Kepler's problem 
    
    Arguments
        t (scalar)
            ??? Is this the period or the time?
            ??? PLZ help
            should be time?-Meg
        z (array-like)
            position and velocity vectors
            position = z[0:1]
            velocity = z[2:3]
        m (scalar)
            total mass in normalized units
    Returns
        numpy array dzdt with components [ dx/dt, dy/dt, dv_x/dt, dv_y/dt ]
    """
    # Fill in the following steps
    # 1. Split z into position vector and velocity vector (view total_energy as an example)
    #r = z[0:2]  # First two arguments of Z are the position vectors
    v = z[2:4]  # Last two is the velocity vector
    # 2. Compute the norm of position vector, and use it to compute the force
    r = np.sqrt(z[0]**2 + z[1]**2) # The direction of the position vector
    
    drdt = v
    dvdt = [(-m/r**3)*z[0], (-m/r**3)*z[1]] #need cubed on bottom #issues here
    # 3. Compute drdt (array [dx/dt, dy/dt])
    # 4. Compute dvdt (array [dvx/dt, dvy/dt])

    # Join the arrays
    dzdt = np.concatenate((drdt,dvdt))
    return dzdt

def integrate_orbit(z0,m,tend,h,method='RK4'):
    """
    Integrates orbit starting from an initial position and velocity from t = 0 
    to t = tend.
    
    Arguments:
        z0 (array)
            starting values for position and velocity

        m (scalar)
            normalized mass
    
        tend (scalar)
            end of time
            note that tstart is always 0.0
    
        h (scalar)
            length of timestep
    
        method ('Euler', 'RK2', or 'RK4')
            identifies which stepper routine to use (default: 'RK4')

    Returns
        ts, Xs, Ys, KEs, PEs, TEs := arrays of time, x postions, y positions, 
        and energies (kin., pot., total) 
    """

    # Set the initial time and phase space array
    t = 0.0
    z = z0

    # Expected number of steps
    Nsteps = int(tend/h)+1

    # Arrays holding t, x, y, kinetic energy, potential energy, and total energy
    ts = np.zeros(Nsteps)
    Xs = np.zeros(Nsteps)
    Ys = np.zeros(Nsteps)
    KEs = np.zeros(Nsteps)
    PEs = np.zeros(Nsteps)
    TEs = np.zeros(Nsteps)

    # Store the initial point
    ts[0] = t
    Xs[0] = z[0]
    Ys[0] = z[1]
    KEs[0] = kinetic_energy(norm(z[2:4]))
    PEs[0] = potential_energy(z[0:2],m)
    TEs[0] = total_energy(z,m)
    

    # Select the stepping method
    advance_one_step = integration_methods[method]
    # Run through the steps
    for step in range(0,Nsteps):
        z = advance_one_step(derivs,t,z,h,args=m)
        # Insert statement here to increment t by the stepsize h
        
        # Store values
        ts[step] = t
        Xs[step] = z[0]
        Ys[step] = z[1]
        KEs[step] = kinetic_energy(norm(z[2:4]))
        PEs[step] = potential_energy(z[0:2],m)
        TEs[step] = total_energy(z,m)
        
        t+= h
    return ts, Xs, Ys, KEs, PEs, TEs
    
def set_initial_conditions(a, m, e):
    """
    set the initial conditions for the orbit.  The orientation of the orbit is 
    chosen so that y0 = 0.0 and vx0 = 0.0.
    
    Arguments
        a
            semi-major axis in AU
        m
            total mass in Msun
        e
            eccentricity ( x0 = (1+e)*a )
    
    Returns:
    x0, y0, vx0, vy0, eps0, Tperiod := initial position and velocity, energy, 
        and period
    """
        
    # Fill in the following lines with the correct formulae
    # Total energy per unit mass
    # period of motion
    Tperiod = (np.pi/(np.sqrt(2)))*(m)*(abs(e))**(-3/2)

    # Initial position
    # Fill in the following lines with the correct formulae
    x0 = (1.5)*a
    y0 = 0.0

    # Initial velocity is in y-direction; we compute it from the energy
    # Fill in the following lines with the correct formulae
    vx0 = 0.0
    vy0 = np.sqrt(2*((e)+(m/x0)))
    
    return np.array([x0,y0,vx0,vy0]), e, Tperiod
