## =================================================================== ##
#  this is file bhaga.py, created at 19-Jan-2012                #
#  maintained by Gustavo Rabello dos Anjos                              #
#  e-mail: gustavo.rabello@gmail.com                                    #
## =================================================================== ##


import numpy as np
g = 9.81
d = 0.0261
rho=1350

print "Eo  M   Re_exp  N             R_fem         mu"

# Oblate ellipsoid
Eo = 116
M = 848
Re_exp = 2.47
N = np.sqrt(Eo**3/M)
R_fem = np.sqrt(N)
mu = np.sqrt( (rho**2 * d**3 * g)/N )
v_exp=Re_exp*mu/(rho*d)
print str(Eo) + "  " + str(M) + "  " + str(Re_exp) + "   " + str(N) +"  " +\
      str(R_fem) + "  " + str(mu) + "  " + str(v_exp)


# Oblate ellipsoid
Eo = 116
M = 266
Re_exp = 3.57
N = np.sqrt(Eo**3/M)
R_fem = np.sqrt(N)
mu = np.sqrt( (rho**2 * d**3 * g)/N )
v_exp=Re_exp*mu/(rho*d)
print str(Eo) + "  " + str(M) + "  " + str(Re_exp) + "   " + str(N) +"  " +\
      str(R_fem) + "  " + str(mu) + "  " + str(v_exp)

# Oblate ellipsoid
Eo = 116
M = 41.1
Re_exp = 7.16
N = np.sqrt(Eo**3/M)
R_fem = np.sqrt(N)
mu = np.sqrt( (rho**2 * d**3 * g)/N )
v_exp=Re_exp*mu/(rho*d)
print str(Eo) + "  " + str(M) + "  " + str(Re_exp) + "   " + str(N) +"  " +\
      str(R_fem) + "  " + str(mu) + "  " + str(v_exp)


# Oblate ellipsoid
Eo = 116
M = 5.51
Re_exp = 13.3
N = np.sqrt(Eo**3/M)
R_fem = np.sqrt(N)
mu = np.sqrt( (rho**2 * d**3 * g)/N )
v_exp=Re_exp*mu/(rho*d)
print str(Eo) + "  " + str(M) + "  " + str(Re_exp) + "   " + str(N) +"  " +\
      str(R_fem) + "  " + str(mu) + "  " + str(v_exp)


# Oblate ellipsoid
Eo = 116
M = 1.31
Re_exp = 20.4
N = np.sqrt(Eo**3/M)
R_fem = np.sqrt(N)
mu = np.sqrt( (rho**2 * d**3 * g)/N )
v_exp=Re_exp*mu/(rho*d)
print str(Eo) + "  " + str(M) + "  " + str(Re_exp) + "   " + str(N) +"  " +\
      str(R_fem) + "  " + str(mu) + "  " + str(v_exp)


# Oblate ellipsoid
Eo = 116
M = 0.103
Re_exp = 42.2
N = np.sqrt(Eo**3/M)
R_fem = np.sqrt(N)
mu = np.sqrt( (rho**2 * d**3 * g)/N )
v_exp=Re_exp*mu/(rho*d)
print str(Eo) + "  " + str(M) + "  " + str(Re_exp) + "   " + str(N) +"  " +\
      str(R_fem) + "  " + str(mu) + "  " + str(v_exp)


# spherical cap (closed wake)
Eo = 115
M = 4.63*10**-3
Re_exp = 94.0
N = np.sqrt(Eo**3/M)
R_fem = np.sqrt(N)
mu = np.sqrt( (rho**2 * d**3 * g)/N )
v_exp=Re_exp*mu/(rho*d)
print str(Eo) + "  " + str(M) + "  " + str(Re_exp) + "   " + str(N) +"  " +\
      str(R_fem) + "  " + str(mu) + "  " + str(v_exp)


# spherical cap (open wake)
Eo = 114
M = 8.60*10**-4
Re_exp = 151
N = np.sqrt(Eo**3/M)
R_fem = np.sqrt(N)
mu = np.sqrt( (rho**2 * d**3 * g)/N )
v_exp=Re_exp*mu/(rho*d)
print str(Eo) + "  " + str(M) + "  " + str(Re_exp) + "   " + str(N) +"  " +\
      str(R_fem) + "  " + str(mu) + "  " + str(v_exp)

