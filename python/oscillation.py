## =================================================================== ##
#  this is file bolha.py, created at 03-May-2011                #
#  maintained by Gustavo Rabello dos Anjos                              #
#  e-mail: gustavo.rabello@gmail.com                                    #
## =================================================================== ##

import numpy as np
from pylab import plot,axis,xlabel,ylabel,title,legend,savefig,show

x=0
y=0

n=2
sigma=1
R=0.5

mu=0.001
rho=1
rho_l=1
rho_g=0.001

nu=mu/rho
tao=R/(5*nu)
w=np.sqrt( 24*sigma/( (3*rho_l+2*rho_g)*R*R*R) )

final=100
yt=np.zeros([final],dtype=float)
yt1=np.zeros([final],dtype=float)
yt2=np.zeros([final],dtype=float)
time=np.zeros([final],dtype=float)

a0=0.01
y0=1.000
t=0
dt=0.01
for i in range(0,final):
 yt[i]=y0 + a0 * np.exp(-t/tao) * np.cos(w*t)
 yt1[i]= y0 + a0 * np.exp(-t/tao)
 yt2[i]= y0 - a0 * np.exp(-t/tao)
 time[i]=t
 t=t+dt
 print time[i],yt[i],yt1[i],yt2[i]

#plot(time,yt)
#xlabel('X')
#ylabel('Y')
#title('bubble shape')
#show()



