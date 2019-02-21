########################################################################
# Team Rabbits: Austin Czyzewski | Megan Davis | Molly Janasik | David Pham
# AST 304, Fall 2018
# Michigan State University
########################################################################

"""
This program documents the routines to advance solution by one timestep.
Must define dz/dt = f(t,z), t, z, h
"""

# All routines that take a single step should have the same interface
def fEuler(f,t,z,h,args=()):
    """
<<<<<<< HEAD
    Forward Euler/ First Order method: 
        
=======
    Forward Euler ODE Solver
    
>>>>>>> 609d61bd47b75474e0ef3cf925b4e2ce9c92d824
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
<<<<<<< HEAD
        t and z
            arguments for f(t,z)
        h is the step  
=======
    
        t (scalar)
            time to evaluate f at
         
        z (array)
            initial position and velocity vectors
            
        h (scalar)
            size of timestep
    
>>>>>>> 17da8785e3fb1ed610d49b1d510e0199ba8233b7
        args (tuple, optional)
            additional arguments to pass to f
    
    Returns
        znew = z(t+h)
    """
    
    # The following trick allows us to pass additional parameters to f
    # First we make sure that args is of type tuple; if not, we make it into that form
    if not isinstance(args,tuple):
        args = (args,)
    
    # When we call f, we use *args to pass it as a list of parameters.
    # For example, if elsewhere we define f like def f(t,z,x,y): ...
    # then we would call this routine as znew = fEuler(f,t,z,h,args=(x,y))
    
    return z + h*f(t,z,*args)

# You will need to flesh out the following routines for a second-order
# Runge-Kutta step and a fourth order Runge-Kutta step.

def rk2(f,t,z,h,args=()):
    """
<<<<<<< HEAD
    Second Order Runge-Kutta method:
        
=======
    Runge-Kutta Second Order ODE Solver
    
>>>>>>> 609d61bd47b75474e0ef3cf925b4e2ce9c92d824
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
<<<<<<< HEAD
        t and z are arguments for f(t,z)
        h is the step  
=======
    
        t (scalar)
            time to evaluate f at
         
        z (array)
            initial position and velocity vectors
            
        h (scalar)
            size of timestep
    
>>>>>>> 17da8785e3fb1ed610d49b1d510e0199ba8233b7
        args (tuple, optional)
            additional arguments to pass to f
    Computes
        zm where zm is the step to the midpoint of the interval
    Returns
        z(t+h), now more accurate
    """
    
    if not isinstance(args,tuple):
        args = (args,)
    
    # Delete the line "pass" when you put in the full routine
    zm= z+(h/2)*(f(t,z,*args))
    return z + h*f(t+(h/2), zm, *args)

def rk4(f,t,z,h,args=()):
    """
<<<<<<< HEAD
    Fourth-order Runge-Kutta method:
        
=======
    Runge-Kutta Forth Order ODE Solver
    
>>>>>>> 609d61bd47b75474e0ef3cf925b4e2ce9c92d824
    Arguments
        f(t,z,...)
            function that contains the RHS of the equation dz/dt = f(t,z,...)
<<<<<<< HEAD
        t and z are arguments for f(t,z)
        h is the step  
=======
    
        t (scalar)
            time to evaluate f at
         
        z (array)
            initial position and velocity vectors
            
        h (scalar)
            size of timestep
    
>>>>>>> 17da8785e3fb1ed610d49b1d510e0199ba8233b7
        args (tuple, optional)
            additional arguments to pass to f
    Computes
        
    Returns
        znew = z(t+h), now even MORE accurate
    """
   
    if not isinstance(args,tuple):
        args = (args,)
    
    zm1 = z + (h/2)*f(t,z, *args)
    zm2 = z + (h/2)*f(t+(h/2), zm1, *args)
    zp= z + h*f(t+(h/2),zm2, *args)
    
    k1 = f(t, z, *args)
    k2 = f(t + h/2, zm1, *args)
    k3 = f(t + h/2, zm2, *args)
    k4 = f(t + h, zp, *args)
    
    return z + h*(k1 + 2*k2 + 2*k3 + k4)/6
    
