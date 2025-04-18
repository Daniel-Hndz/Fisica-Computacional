from numpy import *
import matplotlib.pyplot as plt
import numpy.random as rd

N=16 # numero de particulas
v0=1 # velocidad inicial
L=12 # lado de la caja cuadrada
dt=0.005 #10 ms
d=0.75

def CondIniciales():
     #maxima desviación de su coordenada periodica
    x=arange(1.5,L,3) # distancia entre particulas:3, primer punto a mitad del espacio entre particulas del frontera de la caja:1.5
    r=zeros([N,2],float)

    i=0
    for a in x:
        for b in x:
            r[i,:]=array([a,b], float)
            i+=1
    x=(copy(r[:,0]))+(rd.uniform(-d,d, size=N))

    y=(copy(r[:,1]))+(rd.uniform(-d,d, size=N))

    vx,vy=rd.uniform(-v0,v0, size=N),rd.uniform(-v0,v0, size=N)
    x0,y0=x-vx*dt,y-vy*dt

    return x0,x,vx,y0,y,vy


def Acel(x,y):
    accx=zeros([N,N],float)
    accy=zeros([N,N],float)
    for i in range(N-1):
        for j in range(i+1, N):
            dx=x[j]-x[i]
            dy=y[j]-y[i]
            #Evalución de las distancias con respecto a las condiciones de contorno
            if abs(dx) > L/2: #L/2, proviene de tomar las coordenadas de la particula  como si fueran el centro de la red periodica
                if dx>0:
                    dx-=L
                if dx<0:
                    dx=L-dx
            if abs(dy) > L/2:
                if dx>0:
                    dx-=L
                if dx<0:
                    dx=L-dx
            r=sqrt(dx**2+dy**2)
            if r<=3.5: # 3.5 es la distancia de 'cutoff'
                accx[i,j]=24*((2/(r**13))-(1/(r**7)))*(dx/r)#Calculo de la aceleración por el potencial de Lennard-Jones
                accy[i,j]=24*((2/(r**13))-(1/(r**7)))*(dy/r)
    
    ax=sum(accx,axis=0)
    ay=sum(accy,axis=0)
    return ax,ay


def Update(xp,xc,ax,yp,yc,ay):
    x=2*xc-xp+ax*(dt)**2
    y=2*yc-yp+ay*(dt)**2
    x%=L
    y%=L
    vx=(x-xp)/(2*dt)
    vy=(y-yp)/(2*dt)
    return xc,x,vx,yc,y,vy


T=7 #s
xp,x,vx,yp,y,vy=CondIniciales()
plt.plot(xp,yp,'.', zorder=17)
trax=copy(x)
tray=copy(y)
Rap=sqrt(copy(vx)**2+copy(vy)**2)
for t in range(0,int(T/dt)):
    ax,ay=Acel(x,y)
    xp,x,vx,yp,y,vy=Update(xp,x,ax,yp,y,ay)
    Rap=vstack([Rap,sqrt(vx**2+vy**2)])
    trax=vstack([trax,x])
    tray=vstack([tray,y])



#1-Condiciones iniciales

### Bucle de 0 a t en pasos dt ###
#>  #2-Calcular y almacenar las aceleraciones
#>  #3-Calcular y almacenar nuevas posiciones y velocidades
#>  #4-Evaluar condiciones de contorno
#>  #5-Repetir de 2 a 4 

colores = [ 'g', 'r', 'c', 'm', 'y', 'k', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan', 'lime', 'teal', 'maroon']

for i in range(N):
    plt.plot(trax[:,i], tray[:,i],'.', markersize=1, color=colores[i] )

plt.xlim(0,L)
plt.ylim(0,L)
plt.show()
