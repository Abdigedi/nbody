import math
from matplotlib import pylab
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import mpl_toolkits

G=6.67e-11
M1=1.898e27
M2=2.087e27
M3=2.1827e27
daysec=24*60*60
AU = 1.5e11

epsilon=100000000

#Startposition for masserne
x1 = -0.5*AU
y1 = 0
z1 = 0
x1_liste = []
y1_liste = []
z1_liste = []

x2 = 0.5*AU
y2 = 0
z2 = 0
x2_liste = []
y2_liste = []
z2_liste = []

x3 = 0
y3 = 0.25*AU
z3 = 0.5*AU
x3_liste = []
y3_liste = []
z3_liste = []

#hastigheder

vx1 = 200
vy1 = 1000
vz1 = 1000

vx2 = -1000
vy2 = 500
vz2 = -500

vx3 = -500
vy3 = 200
vz3 = 500

Epot=[]
Ekin_M=[]
Ekin_J=[]
Emek_M=[]
Emek_J=[]
T=[]

t=0
dt=0.1*daysec

def a11(R,r):
    return (G*M2*r)/(R**3)

def a12(R,r):
    return (G*M3*r)/(R**3)

def a21(R, r):
    return (G*M1*r)/(R**3)

def a22(R, r):
    return -(G*M3*r)/(R**3)

def a31(R,r):
    return -(G*M2*r)/(R**3)

def a32(R,r):
    return -(G*M1*r)/(R**3)



def T11(R,r,h):
    k1=a11(R,r)
    k2=a11(R+0.5*h,r+0.5*h*k1)
    k3=a11(R+0.5*h,r+0.5*h*k2)
    k4 = a11(R+h,r+h*k3)
    return (1/6)*(k1+2*k2+2*k3+k4)

def T12(R,r,h):
    k1=a12(R,r)
    k2=a12(R+0.5*h,r+0.5*h*k1)
    k3=a12(R+0.5*h,r+0.5*h*k2)
    k4 = a12(R+h,r+h*k3)
    return (1/6)*(k1+2*k2+2*k3+k4)

def T21(R,r,h):
    k1=a21(R,r)
    k2=a21(R+0.5*h,r+0.5*h*k1)
    k3=a21(R+0.5*h,r+0.5*h*k2)
    k4 = a21(R+h,r+h*k3)
    return (1/6)*(k1+2*k2+2*k3+k4)

def T22(R,r,h):
    k1=a22(R,r)
    k2=a22(R+0.5*h,r+0.5*h*k1)
    k3=a22(R+0.5*h,r+0.5*h*k2)
    k4 = a22(R+h,r+h*k3)
    return (1/6)*(k1+2*k2+2*k3+k4)

def T31(R,r,h):
    k1=a31(R,r)
    k2=a31(R+0.5*h,r+0.5*h*k1)
    k3=a31(R+0.5*h,r+0.5*h*k2)
    k4 = a31(R+h,r+h*k3)
    return (1/6)*(k1+2*k2+2*k3+k4)

def T32(R,r,h):
    k1=a32(R,r)
    k2=a32(R+0.5*h,r+0.5*h*k1)
    k3=a32(R+0.5*h,r+0.5*h*k2)
    k4 = a32(R+h,r+h*k3)
    return (1/6)*(k1+2*k2+2*k3+k4)

#acceleration
while t<5000*daysec:
    r12x=x2-x1
    r12y=y2-y1
    r12z=z2-z1

    r13x=x3-x1
    r13y=y3-y1
    r13z=z3-z1

    r23x=x3-x2
    r23y=y3-y2
    r23z=z3-z2

    R12=math.sqrt((r12x**2) + (r12y**2) + (r12z**2)) + epsilon
    R13=math.sqrt((r13x**2) + (r13y**2) + (r13z**2)) + epsilon
    R23=math.sqrt((r23x**2) + (r23y**2) + (r23z**2)) + epsilon

    V1 = math.sqrt(vx1**2+vy1**2+vz1**2)
    V2 = math.sqrt(vx2**2+vy2**2+vz2**2)
    V3 = math.sqrt(vx3**2+vy3**2+vz3**2)

    a1x = a11(R12, r12x) + a12(R13, r13x)
    a1y = a11(R12, r12y) + a12(R13, r13y)
    a1z = a11(R12, r12z) + a12(R13, r13z)

    a2x = (a21(R12, r12x) + a22(R23, r23x))
    a2y = (a21(R12, r12y) + a22(R23, r23y))
    a2z = (a21(R12, r12z) + a22(R23, r23z))

    a3x = (a31(R13, r13x) + a31(R23, r23x))
    a3y = (a31(R13, r13y) + a31(R23, r23y))
    a3z = (a31(R13, r13z) + a31(R23, r23z))

    h=1
    vx1 = vx1 + T11(R12,r12x,h)*h*dt + T12(R13,r13x,h)*h*dt
    vy1 = vy1 + T11(R12,r12y,h)*h*dt + T12(R13,r13y,h)*h*dt
    vz1 = vz1 + T11(R12,r12z,h)*h*dt + T12(R13,r13z,h)*h*dt

    vx2 = vx2 + T21(R12,r12x,h)*h*dt + T12(R23,r23x,h)*h*dt
    vy2 = vy1 + T21(R12,r12y,h)*h*dt + T12(R23,r23y,h)*h*dt
    vz1 = vz1 + T21(R12,r12z,h)*h*dt + T12(R23,r23z,h)*h*dt

    vx3 = vx3 + T31(R13,r13x,h)*h*dt + T31(R23,r23x,h)*h*dt
    vy3 = vy3 + T31(R13,r13y,h)*h*dt + T31(R23,r23y,h)*h*dt
    vz3 = vz3 + T31(R13,r13z,h)*h*dt + T31(R23,r23z,h)*h*dt

    #position

    x1+=vx1*dt
    y1+=vy1*dt
    z1+=vz1*dt

    x2+=vz2*dt
    y2+=vy2*dt
    z2+=vz2*dt

    x3+=vz3*dt
    y3+=vy3*dt
    z3+=vz3*dt

    x1_liste.append(x1)
    y1_liste.append(y1)
    z1_liste.append(z1)

    x2_liste.append(x2)
    y2_liste.append(y2)
    z2_liste.append(z2)

    x3_liste.append(x3)
    y3_liste.append(y3)
    z3_liste.append(z3)

    # Epot.append(-(G*Mj*Mm)/R)
    # Ekin_J.append(0.5*Mj*(Vj**2))
    # Ekin_M.append(0.5*Mm*Vm**2)
    # Emek_J.append((0.5*Mj*Vj)**2-(G*Mj*Mm)/R)
    # Emek_M.append(0.5*Mm*(Vm**2)-((G*Mj*Mm)/R))


    t+=dt
    T.append(t)

ax=plt.axes(projection='3d')
ax.plot3D(x1_liste, y1_liste, z1_liste, 'r-', linewidth = 2, label="M1")
ax.plot3D(x2_liste, y2_liste, z2_liste, 'g-', linewidth = 2, label="M2")
ax.plot3D(x3_liste, y3_liste, z3_liste, 'b-', linewidth = 2, label="M3")
leg = plt.legend(loc='upper left', frameon=True)
plt.show()

# plt.plot(T,Ekin_M, "-g", label='Ekin')
# plt.plot(T,Epot, "-b", label='Epot')
# plt.plot(T,Emek_M, "-r", label='Emek' )
# leg = plt.legend(loc='upper left', frameon=True)
#
#
# pylab.xlabel("Tid (s)")
# pylab.ylabel("E (J)")
# plt.show()
