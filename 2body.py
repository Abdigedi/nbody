import math
from matplotlib import pylab
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import mpl_toolkits

G=6.67e-11
Mj=5.9e24
Mm=7.3e22
daysec=24*60*60
AU = 1.5e11
#data for jorden
xj = 0
yj = 0
zj = 0

vxj=0
vyj=0
vzj=0

xj_liste=[]
yj_liste=[]
zj_liste=[]

#for månen:
xm = 384000000
ym = 0
zm = 0

vxm=0
vym=1000
vzm=0

xm_liste=[]
ym_liste=[]
zm_liste=[]

Epot=[]
Ekin_M=[]
Ekin_J=[]
Emek_M=[]
Emek_J=[]
T=[]

t=0
dt=0.1*daysec

def am(R,r):
    return -(G*Mj*r)/(R**3)

def aj(R,r):
    return (G*Mm*r)/(R**3)

def Tj(R,r,h):
    k1=aj(R,r)
    k2=aj(R+0.5*h,r+0.5*h*k1)
    k3=aj(R+0.5*h,r+0.5*h*k2)
    k4 = aj(R+h,r+h*k3)
    return (1/6)*(k1+2*k2+2*k3+k4)

def Tm(R,r,h):
    k1=am(R,r)
    k2=am(R+0.5*h,r+0.5*h*k1)
    k3=am(R+0.5*h,r+0.5*h*k2)
    k4 = am(R+h,r+h*k3)
    return (1/6)*(k1+2*k2+2*k3+k4)



#acceleration
while t<150*daysec:
    rx=xm-xj
    ry=ym-yj
    rz=zm-zj

    R=math.sqrt((rx**2) + (ry**2) + (rz**2))
    Vm = math.sqrt(vxm**2+vym**2+vzm**2)
    Vj = math.sqrt(vxj**2+vyj**2+vzj**2)

    amx = am(R,rx)
    amy = am(R,ry)
    amz = am(R,rz)

    ajx = aj(R,rx)
    ajy = aj(R,ry)
    ajz = aj(R,rz)

    h=1
    vxm+=Tm(R,rx,h)*h*dt
    vym+=Tm(R,ry,h)*h*dt
    vzm+-Tm(R,rz,h)*h*dt

    vxj+=Tj(R,rx,h)*h*dt
    vyj+=Tj(R,ry,h)*h*dt
    vzj+=Tj(R,rz,h)*h*dt

    #position

    xm+=vxm*dt
    ym+=vym*dt
    zm+=vzm*dt

    xj+=vzj*dt
    yj+=vyj*dt
    zj+=vzj*dt

    xm_liste.append(xm)
    ym_liste.append(ym)
    zm_liste.append(zm)

    xj_liste.append(xj)
    yj_liste.append(yj)
    zj_liste.append(zj)

    Epot.append(-(G*Mj*Mm)/R)
    Ekin_J.append(0.5*Mj*(Vj**2))
    Ekin_M.append(0.5*Mm*Vm**2)
    Emek_J.append((0.5*Mj*Vj)**2-(G*Mj*Mm)/R)
    Emek_M.append(0.5*Mm*(Vm**2)-((G*Mj*Mm)/R))


    t+=dt
    T.append(t)

ax=plt.axes(projection='3d')
ax.plot3D(xm_liste, ym_liste, zm_liste, linewidth = 2, label="Månen")
ax.plot3D(xj_liste, yj_liste, zj_liste, linewidth = 2, label="Jorden")
plt.show()

plt.plot(T,Ekin_M, "-g", label='Ekin')
plt.plot(T,Epot, "-b", label='Epot')
plt.plot(T,Emek_M, "-r", label='Emek' )
leg = plt.legend(loc='upper left', frameon=True)


pylab.xlabel("Tid (s)")
pylab.ylabel("E (J)")
plt.show()
