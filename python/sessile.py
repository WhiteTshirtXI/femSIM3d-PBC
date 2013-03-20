## =================================================================== ##
#  this is file bolha.py, created at 03-May-2011                #
#  maintained by Gustavo Rabello dos Anjos                              #
#  e-mail: gustavo.rabello@gmail.com                                    #
## =================================================================== ##

import numpy as np
from pylab import plot,axes,axis,xlabel,ylabel,title,legend,savefig,show

# r = x
x=0.258
z=0.0

ds=1.0/200
Eo=2*0.9

p=2.75
theta=0+0.0

n=100
X=np.zeros([n],dtype=float)
Z=np.zeros([n],dtype=float)
for i in range(0,n):
 theta = theta + (Eo*(p-z) - np.sin(theta)/x )*ds
 x=x+np.cos(theta)*ds
 z=z+np.sin(theta)*ds
 X[i] = x
 Z[i] = z
 print theta,x,z

#--------------------------------------------------
# # read data
# read_data = np.loadtxt("../dist.dat")
# XA = read_data[:,0]
# YA = read_data[:,1]
# ZA = read_data[:,2]
#-------------------------------------------------- 

#--------------------------------------------------
# plot(X,Z)
# plot(YA,ZA+3.97,marker='o')
# xlabel('X')
# ylabel('Z')
# title('bubble shape')
# axes().set_aspect('equal')
# show()
#-------------------------------------------------- 



# reta y* x kappa* (gnuplot)
# set size ratio 1
# f(x) = 2*2.56-2.2*x
#plot 'dist.dat' using ($3+3.97):4,f(x)

# bolha analitica (gnuplot)
# set size ratio -1
#plot 'sessile3d' using 2:3 with lines,'dist.dat' using 2:($3+3.97)
