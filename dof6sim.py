    # 
    # see military handbook for missile flight simulation ch.12 simulation synthesis (205)
    # 
    # state vector xs variables:
    #   0   x
    #   1   y
    #   2   z
    #   3   u
    #   4   v
    #   5   w
    #   6   phi
    #   7   theta
    #   8   psi
    #   9   p
    #   10  q
    #   11  r
    ##
from scipy.integrate import solve_ivp 
import numpy as np

import sys, os, importlib
sys.path.append(os.path.join(os.getcwd(), '..'))

import C4dynamics as c4d
importlib.reload(c4d)

from control_system import control_system
importlib.reload(control_system)

#
# 
##
missile = c4d.rigidbody()
target = c4d.datapoint(x = 4000, y = 1000, z = -3000
                        , vx = -250, vy = 0, vz = 0)
seeker = c4d.seekers.lineofsight(tau1 = 0.01, tau2 = 0.01)
ctrl = control_system()

t = 0
dt = 5e-3
tf = 60
vm = 30
# x = np.zeros(12)

# 
# atmospheric properties up to 2000m
##
pressure = 101325 # pressure pascals
rho = 1.225       # density kg/m^3
vs = 340.29       # speed of sound m/s

# 
# parameters for initial example
##
Gn = 250 

s = 0.0127
d = 0.127
mach = 0.8

cD0 = 0.8
cLa = 39
cMa = -170 
cMd = 250
cMqcMadot = -13000

k = 0.0305

m = 57
xcm = 1.35
xref = 1.35

g = 9.8


ixx = 1
iyy = izz = 61 

#
# init
#
# The initial missile pointing direction and angular rates are calculated in the fire-control block. 
# For the example simulation, a simple algorithm is employed in which the missile is
# pointed directly at the target at the instant of launch, and missile angular rates at launch are assumed to be negligible.
# The unit vector uR in the direction from the missile to the target is calculated by normalizing the range vector R.
## 
rTM = target.pos() - missile.pos()
ucl = rTM / np.linalg.norm(rTM) # center line unit vector 
missile.vx, missile.vy, missile.vz = vm * ucl 
missile.psi = np.arctan(ucl[1] / ucl[0])
missile.theta = np.arctan(-ucl[2] / np.sqrt(ucl[1]**2 + ucl[0]**2))
missile.phi = 0
u, v, w = missile.BI() @ missile.vel()
                
def eqm(t, xs, rb): 
    
    x, y, z, u, v, w, phi, theta, psi, p, q, r = xs





    #
    # calc
    ## 


    BI = c4d.dcm321(phi, theta, psi)
    ucl = BI @ [[1], [0], [0]] # is it? 
    ucl = np.array([[np.cos(theta) * np.cos(psi)]
            , [np.cos(theta) * np.sin(psi)] 
            , [np.sin(-theta)]])
    vm = np.transpose(BI) @ [[u], [v], [w]]
    vm_total = np.sqrt(vm[0]**2 + vm[0]**2 + vm[0]**2)

    alpha = np.arctan2(w, u)
    beta  = np.arctan2(-v, u)
    uvm = vm / vm_total
    alpha_total = 0# np.arccos(uvm @ ucl)

    # 
    # aerodynamic forces
    ##



    # 
    # guidance and control
    ## 
    # d1, d2, d3, d4 = 0, 0, 0, 0
    acmd_yb, acmd_zb = 0, 0
    dpitch = -Gn * acmd_zb / Q
    dyaw = -Gn * acmd_yb / Q




    # lift and drag
    cL = cLa * alpha_total
    L = Q * s * cL 

    cD = cD0 + k * cL**2
    D = Q * s * cD

    # in body frame
    A = D * np.cos(alpha_total) - L * np.sin(alpha_total)
    N = D * np.sin(alpha_total) + L * np.cos(alpha_total)

    fAb = np.array([[-A]
                    , [N * (-v / np.sqrt(v**2 + w**2))]
                    , [N * (-w / np.sqrt(v**2 + w**2))]])

    cNy = fAb[1] / Q / s
    cNz = fAb[2] / Q / s

    # 
    # aerodynamic moments 
    ## 

    cNb = cMa
    cNd = cMd 
    # cNr = cMq 
    # cNbdot = cMadot 
    cNrcNbdot = cMqcMadot
    cMref = cMa * alpha + cMd * dpitch
    cNref = cNb * beta  + cNd * dyaw


    # wrt center of mass
    cM = cMref - cNz * (xcm - xref) / d + d / (2 * v) * cMqcMadot * q
    cN = cNref - cNy * (xcm - xref) / d + d / (2 * v) * cNrcNbdot * r

    lA = 0              # aerodynamic moemnt in roll
    mA = Q * cM * s * d # aerodynamic moment in pitch
    nA = Q * cN * s * d # aerodynamic moment in yaw 


    # 
    # gravity
    ## 
    fGe = [[0], [0], [m * g]]
    fGb = BI @ fGe 



    #
    # translational motion derivatives
    ##
    dx = vm[0]
    dy = vm[1]
    dz = vm[2]

    du = (fAb[0] + fGb[0]) / m - (q * w - r * v)
    dv = (fAb[1] + fGb[1]) / m - r * u
    dw = (fAb[2] + fGb[2]) / m + q * u



    # 
    # euler angles derivatives 
    ## 

    dphi   = (q * np.sin(phi) + r * np.cos(phi)) * np.tan(theta)
    dtheta =  q * np.cos(phi) - r * np.sin(phi)
    dpsi   = (q * np.sin(phi) + r * np.cos(phi)) / np.cos(theta)

    # 
    # angular motion derivatives 
    ## 
    dp     = (lA - q * r * (izz - iyy)) / ixx
    dq     = (mA - p * r * (ixx - izz)) / iyy
    dr     = (nA - p * q * (iyy - ixx)) / izz

    return dx, dy, dz, du, dv, dw, dphi, dtheta, dpsi, dp, dq, dr


while t <= tf:

    vm = missile.V()
    ucl = missile.BI() @ np.array([1, 0, 0]) # unit centerline vector
    
    #
    # atmospheric calculations    
    ##
    h = -missile.z   # missile altitude above sea level, m
    mach = vm / vs # mach number 
    Q = 1 / 2 * rho * vm**2 # dynamic pressure 
    
    # 
    # relative position
    ##
    vTM = target.vel() - missile.vel() # missile-target relative velocity 
    rTM = target.pos() - missile.pos() # relative position 
    uR = rTM / np.linalg.norm(rTM) # unit range vector 
    vc = -uR * vTM # closing velocity 

    # 
    # seeker 
    ## 
    wf = seeker.measure(rTM, vTM)
    
    # 
    # guidance and control 
    ##
    Gs = 4 * vm
    acmd = Gs * np.cross(wf, ucl)
    ab_cmd = missile.BI() @ acmd 
    afp, afy = ctrl.upate(ab_cmd, Q)
    dpitch = afp - aoa 
    dyaw = afy - aos  



    x = missile.x, missile.y, missile.z, u, v, w, missile.phi, missile.theta, missile.psi, missile.p, missile.q, missile.r
    x = solve_ivp(eqm, [t, t + missile._dt], x, args = (missile, )).y[:, -1]
    
    # https://www.mathworks.com/matlabcentral/answers/411616-how-to-use-a-for-loop-to-solve-ode
    # Y = [y0 zeros(length(y0), length(tspan))];
    # for i=1:length(tspan)
    #     ti = tspan(i); yi = Y(:,i);
    #     k1 = f(ti, yi);
    #     k2 = f(ti+dt/2, yi+dt*k1/2);
    #     k3 = f(ti+dt/2, yi+dt*k2/2);
    #     k4 = f(ti+dt  , yi+dt*k3);
    #     dy = 1/6*(k1+2*k2+2*k3+k4);
    #     Y(:,i+1) = yi +dy;
    
    missile.x, missile.y, missile.z, u, v, w, missile.phi, missile.theta, missile.psi, missile.p, missile.q, missile.r = x
    t += dt
    missile.store(t)
    
missile.draw('x')
missile.draw('top')
missile.draw('z')
    
    
    
    
    
    
    
    
    
    
    
    
    




    
    
    
    
    
    
    