import numpy as np
from firedrake import *

EgMax = 6

# **************** Example 0 ************************
# **************** Example 0 ************************
# **************** Example 0 ************************

class Example(object):
    def __init__(self, ss):
        self.desc = ss
        return

eg0 = Example("Guess Flow for all zero or perturbation on velocity")

def initial0():
    eg0.PureDirichletBc = True;
    eg0.eta = 1;
    rho = 1.0

    bdy0 = [1,2,3,4]   # Dirichlet data at x=0 and y = 0
    bdy1 = []   # Neumann   data at x=1 and y = 1
 
    penalty = 0.1             # Penatly for |d|^2 = 1
    alpha = [0,0,1,1,0,0]     # Leslie Constants
    kk = [1,1,1,0,0]

    return (rho, bdy0, bdy1, penalty, kk, alpha)

def exact0(t,xx):
    c0 = 0.1
    z = xx[0] - xx[0]
    
#     v   = np.array([z,z])
#     dvt = np.array( [z,z])
#     dvx = np.array([[z,z], [z,z]])
 
    
    v   = np.array([z, c0*(1 - xx[0]*xx[0])])
    dvt = np.array( [z,z])
    dvx = np.array([[z,-2*c0*xx[0]], [z,z]])

    mu = 1.0
    
    p =z; dpt = z; dpx = np.array([z,z])

    d   = np.array([z,z])
    ddt = np.array( [z,z])
    ddx = np.array([[z,z], [z,z]])
    return ([v,p,d], [dvx,dpx,ddx], [dvt,dpt,ddt])

def rhs0(t, xx):
    z = xx[0] - xx[0]
    mu = 1.0
    return ([0,0],[0,0])
   



eg0.initialize = initial0; eg0.exact = exact0; eg0.rhs = rhs0
##################################################
# ****************Symmetric flow and antisymmetric director with eta<0 and eta>0************************

# **************** Example 1 ************************
eg1 = Example("Guess Flow For perturbation on direction with eta<0 and eta>0")

def initial1():
    eg1.PureDirichletBc = False;
    eg1.eta = 1;
    #    eg1.eta = -1;

    rho = 1.0

    bdy0 = [1,2,3,4]   # Dirichlet data at x=0 and y = 0
    bdy1 = []   # Neumann   data at x=1 and y = 1
 
    penalty = 0.1             # Penatly for |d|^2 = 1
    alpha = [0,0,1,1,0,0]     # Leslie Constants
    kk = [0,1,0,0,0]

    return (rho, bdy0, bdy1, penalty, kk, alpha)

def exact1(t,xx):
    z = xx[0] - xx[0]
 
    c0 = 0.2
    theta  = c0 * sin(2*pi*xx[0])
    
    

    v = np.array( [z,z])
    dvt = np.array( [z,z])
    dvx = np.array([[z,z], [z,z]])


    mu = 1.0
    
    p = z
    dpt = z;
    dpx = np.array([z,z])

    d   = np.array([sin(theta), cos(theta)])
    ddt = np.array( [z,z])
    ddx = np.array([[z,z], [z,z]])
    return ([v,p,d], [dvx,dpx,ddx], [dvt,dpt,ddt])

def rhs1(t, xx):
    z = xx[0] - xx[0]

    eta = eg1.eta
    return ([z,z],[z,z])
 
   #return ([-eta*0.100e-1*sin(2*pi*x)*pi*cos(2*pi*x),(-eta*0.5e-3*cos(2*pi*x)**2-eta*0.995e-1)*pi*cos(2*pi*x)/(0.25e-2*cos(2*pi*x)**2+.9975)**(1/2)], [z,z])
   



eg1.initialize = initial1; eg1.exact = exact1; eg1.rhs = rhs1

####################################
# **************** Example 2 ************************
# ****************Symmetric flow and antisymmetric director************************
eg2 = Example("Flow For perturbation on direction with eta<0 ")

def initial2():
    eg2.PureDirichletBc = False;
    eg2.eta = -1;

    rho = 1.0

    bdy0 = [1,2,3,4]   # Dirichlet data at x=0 and y = 0
    bdy1 = []   # Neumann   data at x=1 and y = 1
 
    penalty = 0.1             # Penatly for |d|^2 = 1
    alpha = [0,0,1,1,0,0]     # Leslie Constants
    kk = [0,1,0,0,0]

    return (rho, bdy0, bdy1, penalty, kk, alpha)

def exact2(t,xx):
    z = xx[0] - xx[0]
 
    c0 = 0.2
    theta  = c0 * sin(2*pi*xx[0])
    
    

    v = np.array( [z,z])
    dvt = np.array( [z,z])
    dvx = np.array([[z,z], [z,z]])


    mu = 1.05B
    
    p = z
    dpt = z;
    dpx = np.array([z,z])

    d   = np.array([sin(-theta), cos(-theta)])
    ddt = np.array( [z,z])
    ddx = np.array([[z,z], [z,z]])
    return ([v,p,d], [dvx,dpx,ddx], [dvt,dpt,ddt])

def rhs2(t, xx):
    z = xx[0] - xx[0]

    eta = eg2.eta
    return ([z,z],[z,z])
 



eg2.initialize = initial2; eg2.exact = exact2; eg2.rhs = rhs2



####################################
####################################
####################################
# **************** Example 3 ************************
# ****************Symmetric director and antisymmetric flow
eg3 = Example("Flow For perturbation on direction with eta>0 ")

def initial3():
    eg3.PureDirichletBc = False;
    eg3.eta = 1;

    rho = 1.0

    bdy0 = [1,2,3,4]   # Dirichlet data at x=0 and y = 0
    bdy1 = []   # Neumann   data at x=1 and y = 1
 
    penalty = 0.1             # Penatly for |d|^2 = 1
    alpha = [0,0,1,1,0,0]     # Leslie Constants
    kk = [0,1,0,0,0]

    return (rho, bdy0, bdy1, penalty, kk, alpha)

def exact3(t,xx):
    z = xx[0] - xx[0]
 
    c0 = 0.1
    theta  = c0 * sin(pi*xx[0])
    
    

    v = np.array( [z,z])
    dvt = np.array( [z,z])
    dvx = np.array([[z,z], [z,z]])


    mu = 1.0
    
    p = z
    dpt = z;
    dpx = np.array([z,z])

    d   = np.array([sin(-theta), cos(-theta)])
    ddt = np.array( [z,z])
    ddx = np.array([[z,z], [z,z]])
    return ([v,p,d], [dvx,dpx,ddx], [dvt,dpt,ddt])

def rhs3(t, xx):
    z = xx[0] - xx[0]

    eta = eg3.eta
    return ([z,z],[z,z])
 



eg3.initialize = initial3; eg3.exact = exact3; eg3.rhs = rhs3
#####################################
Examples = [eg0, eg1, eg2, eg3]

def mkWofRel(penalty, kk, alpha):
    "mkWofRel(pentalty for |d|=1, kk=[k1,k2,k3,k4,q], alpha=[alpha1,..,alpha5])"

    import numpy as np

    print("mkWofRel() in 2d. Assuming Parodi Condition")
    print("    mkWofRel() penalty for |d| = 1 ........ {}".format(penalty))
    print("    mkWofRel() chiral constant q .......... {}".format(kk[4]))
    print("    mkWofRel() Frank constants kk ......... {}".format(kk[0:4]))
    print("    mkWofRel() Leslie constants alpha ..... {}".format(alpha[0:5]))

    k1 = kk[0]; k2 = kk[1]; k3 = kk[2]; k4 = kk[3]  # keep usual names
    
    def Wof(dd, gradd):
        "Wof(dd=director, gradd=grad(dd))"
        
        # Maple Generated Code
        
        wof = .25*(-dd[0]**2-dd[1]**2+1)**2/penalty+.5*k1*(gradd[0, 0]+gradd[1, 1])**2+.5*(k2-k4)*(gradd[0, 0]**2+gradd[0, 1]**2+gradd[1, 0]**2+gradd[1, 1]**2-(gradd[0, 0]+gradd[1, 1])**2-(gradd[1, 0]-gradd[0, 1])**2)+.5*k3*(gradd[1, 0]-gradd[0, 1])**2

        dWdd = [-1.00*(-dd[0]**2-dd[1]**2+1)*dd[0]/penalty, -1.00*(-dd[0]**2-dd[1]**2+1)*dd[1]/penalty]

        dWdgd = [[k1*(gradd[0, 0]+gradd[1, 1])-(k2-k4)*gradd[1, 1], (k2-k4)*gradd[1, 0]-k3*(gradd[1, 0]-gradd[0, 1])], [(k2-k4)*gradd[0, 1]+k3*(gradd[1, 0]-gradd[0, 1]), k1*(gradd[0, 0]+gradd[1, 1])-(k2-k4)*gradd[0, 0]]]
      
        return (wof, np.array(dWdd), np.array(dWdgd))

    # END Wof()

    # Use Virga's notation (a4 = mu)

    a1 = alpha[0]; a2 = alpha[1]; a3 = alpha[2]; a4 = alpha[3]; a5 = alpha[4]
    
    g1 = a3 - a2
    g2 = a3 + a2
    g3 = a3 + a2 + 2*a5
    g4 = a1
    mu = a4

    print("    mkWofRel() Raleigan coefficients gamma  {}".format(
        [g1,g2,g3,g4,mu]))

    def Rel(dd, ddot, gradv):
        "Rel(dd=director, ddot=dd_t + (v.grad)dd, gradv=grad(v))"

        rel = .5*g1*((ddot[0]-(.5*gradv[0, 1]-.5*gradv[1, 0])*dd[1])**2+(ddot[1]-(-.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0])**2)+g2*((gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*(ddot[0]-(.5*gradv[0, 1]-.5*gradv[1, 0])*dd[1])+((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*(ddot[1]-(-.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]))+.5*g3*((gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])**2+((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])**2)+.5*g4*((gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*dd[0]+((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*dd[1])**2+.5*mu*(gradv[0, 0]**2+2*(.5*gradv[1,0]+.5*gradv[0,1])**2+gradv[1, 1]**2)
        
        dRddt =[.5*g1*(2*ddot[0]-2*(.5*gradv[0, 1]-.5*gradv[1, 0])*dd[1])+g2*(gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1]), .5*g1*(2*ddot[1]-2*(-.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0])+g2*((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])]
        
        dRdgv = [[g2*dd[0]*(ddot[0]-(.5*gradv[0, 1]-.5*gradv[1, 0])*dd[1])+g3*(gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*dd[0]+g4*((gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*dd[0]+((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*dd[1])*dd[0]**2+mu*gradv[0, 0], .5*g1*(-(ddot[0]-(.5*gradv[0, 1]-.5*gradv[1, 0])*dd[1])*dd[1]+(ddot[1]-(-.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0])*dd[0])+g2*(.5*(ddot[0]-(.5*gradv[0, 1]-.5*gradv[1, 0])*dd[1])*dd[1]-.5*(gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*dd[1]+.5*(ddot[1]-(-.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0])*dd[0]+.5*((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*dd[0])+.5*g3*((gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*dd[1]+((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*dd[0])+1.00*g4*((gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*dd[0]+((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*dd[1])*dd[1]*dd[0]+.5*mu*(gradv[0, 1]+gradv[1,0])], [.5*g1*((ddot[0]-(.5*gradv[0, 1]-.5*gradv[1, 0])*dd[1])*dd[1]-(ddot[1]-(-.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0])*dd[0])+g2*(.5*(ddot[0]-(.5*gradv[0, 1]-.5*gradv[1, 0])*dd[1])*dd[1]+.5*(gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*dd[1]+.5*(ddot[1]-(-.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0])*dd[0]-.5*((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*dd[0])+.5*g3*((gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*dd[1]+((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*dd[0])+1.00*g4*((gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*dd[0]+((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*dd[1])*dd[1]*dd[0]+.5*mu*(gradv[0, 1]+gradv[1,0]), g2*dd[1]*(ddot[1]-(-.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0])+g3*((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*dd[1]+g4*((gradv[0, 0]*dd[0]+(.5*gradv[0, 1]+.5*gradv[1, 0])*dd[1])*dd[0]+((.5*gradv[0, 1]+.5*gradv[1, 0])*dd[0]+gradv[1, 1]*dd[1])*dd[1])*dd[1]**2+mu*gradv[1, 1]]]

        
        return (rel, np.array(dRddt), np.array(dRdgv))

    return (Wof, Rel)

 # Using Sympy doesn't work;
 # subs(.) won't substitute in non-sympy object (e.g. Fenics expressions)

def mkWofRelSymPy(penalty, kk, alpha):  
    "mkWofRel(pentalty for |d|=1, kk=[k1,k2,k3,k4,q], alpha=[alpha1,..,alpha5])"
    import numpy as np
    import sympy as sym

    print("mkWofRel() Assuming Parodi Condition")
    print("    mkWofRel() penalty for |d| = 1 ........ {}".format(penalty))
    print("    mkWofRel() chiral constant q .......... {}".format(kk[4]))
    print("    mkWofRel() Frank constants kk ......... {}".format(kk[0:4]))
    print("    mkWofRel() Leslie constants alpha ..... {}".format(alpha[0:5]))

    # d = director, v = velocity

    d  = np.array(sym.symbols('d0:2'))
    dt = np.array(sym.symbols('dt0:2'))                # d_t + (v.grad) d
    gd = np.reshape(sym.symbols('gd0:2(0:2)'), (2,2))  # grad(d)

    gv = np.reshape(sym.symbols('gv0:2(0:2)'), (2,2))  # grad(v)
    dv = 0.5 * (gv + np.transpose(gv))                 # D(v)
    wv = 0.5 * (gv - np.transpose(gv))                 # W(v)

    # Oseen Frank Energy without chiral term
    
    curld = gd[1,0] - gd[0,1];    # 2d curl = third component
    divd  = gd[0,0] + gd[1,1];

    # penalty = sym.symbols('penalty')
    # k1,k2,k3,k4 = sym.symbols('k1,k2,k3,k4')

    k1 = kk[0]; k2 = kk[1]; k3 = kk[2]; k4 = kk[3]  # keep usual names

    rng = range(0,2)

    ww = (1.0/(4.0*penalty)) * (1-d[0]**2-d[1]**2)**2 \
          + (k1/2.0) * divd**2  \
          + (k2-k4)/2.0 * (np.sum(gd[i,j]**2 for i in rng for j in rng)
                         - divd**2 - curld**2) \
          + (k3/2.0) * curld**2

    dWdd = np.array([sym.diff(ww, d[i]) for i in rng])

    dWdgd = np.array([[sym.diff(ww, gd[i,j]) for j in rng] for i in rng])

    def Wof(dd, gradd):
        "Wof(dd=director, gradd=grad(dd))"
        ss = [(d[i],dd[i]) for i in rng] \
             + [(gd[i,j], gradd[i,j]) for i in rng for j in rng]

        return (ww.subs(ss),
                np.array([dWdd[i].subs(ss) for i in rng]),
                np.array([[dWdgd[i,j].subs(ss) for j in rng] for i in rng]))

    # END Wof()

    # a1,a2,a3,a4,a5,a6 = sym.symbols('a1,a2,a3,a4,a5,a6')
    # g1,g2,g3,g4,mu = sym.symbols('g1,g2,g3,g4,mu')

    # Use Virga's notation (a4 = mu)

    a1 = alpha[0]; a2 = alpha[1]; a3 = alpha[2]; a4 = alpha[3]; a5 = alpha[4]
    
    g1 = a3 - a2
    g2 = a3 + a2
    g3 = a3 + a2 + 2.0*a5
    g4 = a1
    mu = a4

    dcirc = dt - np.dot(wv,d)
    dvd   = np.dot(dv,d)
        
    rr = 0.5*(g1*np.dot(dcirc,dcirc) + g2*np.dot(dvd,dcirc)
              + g3*np.dot(dvd,dvd) + g4*np.dot(dvd,d)**2
              + mu*np.sum(dv[i,j]**2 for i in rng for j in rng) )

    dRddt = np.array([sym.diff(rr, dt[i]) for i in rng])

    dRdgv = np.array([[sym.diff(rr, gv[i,j]) for j in rng] for i in rng])

    def Rel(dd, ddot, gradv):
        "Rel(dd=director, ddot=dd_t + (v.grad)dd, gradv=grad(v))"
        ss = [(d[i], dd[i]) for i in rng] \
             + [(dt[i], ddot[i]) for i in rng] \
             + [(gv[i,j], gradv[i,j]) for i in rng for j in rng]

        return (rr.subs(ss),
                np.array([dRddt[i].subs(ss) for i in rng]),
                np.array([[dRdgv[i,j].subs(ss) for j in rng] for i in rng]))

    # END Rel()


    return (Wof, Rel)
    
    