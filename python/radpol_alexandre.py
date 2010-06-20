#! /usr/bin/env python
# -*- coding: utf-8 -*-

######### Modules #########################################
from pylab import *
from scipy.special import j0,jn

###########################################################
# Champ d'après 
# A.April, Opt. Lett. 331563 (2008)
lambda0 = 0.8e-6
k = 2.0*pi/lambda0
a = 10.0/k
wo = (sqrt(2.0)/k)*sqrt(sqrt(1.0+(k*a)**2)-1.0);

print
print "*******************"
print "lambda0    = " + str(lambda0)
print "ka         = " + str(k*a)
print "wo/lambda0 = " + str(wo/lambda0)
print "*******************"
print

###########################################################
lim = 10.0
if (k*a <= 100): lim = 6.0
if (k*a <= 70): lim = 4.0
if (k*a <= 30): lim = 3.0
if (k*a <= 10): lim = 2.0
lim *= lambda0
np = 1000
r = linspace(0,lim,np)
dr = lim/np
z = 0

R = sqrt(r**2 + (z + a*1j)**2)
theta = arccos((z+a*1j)/R)
kR = k*R

###########################################################
# Champ électrique du TM01 [Eqs. (6) de A.April, Opt. Lett. 331563 (2008)]
# E0 est défini à partir de la puissance 
# [Eq. (12) de A. April, J. Opt. Soc. Am. A 27, 76-81 (2010)]
eta0 = 120.0*pi
P = 1.0e12
E0 = sqrt(P*exp(2.0*k*a)*(4.0*eta0*k**2)/pi\
     /(2.0*k*a*sinh(2.0*k*a) - cosh(2.0*k*a) + 1.0 - 2.0*(k*a)**2))
print
print "*******************"
print "P = %.3e W" %P
print "E0 = %.3e V/m" %E0
print
     
H0 = E0/eta0
H0p = 0.5*H0*k**3*wo**3
E0p = eta0*H0p

Hphi = H0p*sqrt(pi/kR)*jn(1+0.5,kR)*sin(theta)*exp(-k*a)
Er = -0.5*1j*E0p*jn(2+0.5,kR)*sin(2.0*theta)*sqrt(pi/kR)*exp(-k*a)
Ez = -2.0/3.0*1j*E0p*(jn(0+0.5,kR) + jn(2+0.5,kR)*0.25*(1.0+3.0*cos(2.0*theta)))\
       *sqrt(pi/kR)*exp(-k*a)

###########################################################
# Densité d'énergie
mu0 = 12.566370614e-7
eps0 = 8.854187817e-12
wr = 0.5*eps0*abs(Er)**2
wz = 0.5*eps0*abs(Ez)**2
wphi = 0.5*mu0*abs(Hphi)**2

###########################################################
# Figure
figure(figsize=(16.0,8.5))
tsize = 16
maxi = max(wr + wz)
rscl = sqrt(R**2 - (z + a*1j)**2)/lambda0
title(r"$w_0 = %f \lambda_0,\quad ka = %f$" %(wo,k*a),size = tsize)
plot(rscl,wr/maxi,label=r"$|E_r|^2$")
plot(rscl,wz/maxi,label=r"$|E_z|^2$")
plot(rscl,(wr+wz)/maxi,label=r"$|E_r|^2 + |E_z|^2$")
plot(rscl,wphi/maxi,label=r"$|H_\phi|^2$")
xlabel(r"$r/\lambda$",size=tsize)
ylabel(unicode("Champ au carré"),size=tsize)
ylim(0,1)
xlim(0,lim/lambda0)
grid(True)
legend(loc="best")

###########################################################
# Puissance
# Ref: A. April, J. Opt. Soc. Am. A 27, 76-81 (2010)
# Eq. (12)
P = pi/(4.0*eta0*k**2)*E0**2*exp(-2.0*k*a)\
    *(2.0*k*a*sinh(2.0*k*a) - cosh(2.0*k*a) + 1.0 - 2.0*(k*a)**2)

Pr = Pz = Pphi = 0
for i in range(0,np):
     S = real(Er[i]*Hphi[i].conjugate())*r[i]
     Pr += wr[i]*r[i]
     Pz += wz[i]*r[i]
     Pphi += wphi[i]*r[i]

S*=0.5*dr*2.0*pi
Pr*=dr*2.0*pi
Pz*=dr*2.0*pi
Pphi*=dr*2.0*pi
Ptot = Pr + Pz + Pphi
ratio1 = real(Pr/Pz)
ratio2 = real((Pr+Pz)/Pphi)

print
print "*******************"
print "Pr   = %.3e " %real(Pr)
print "Pz   = %.3e " %real(Pz)
print "Pphi = %.3e " %real(Pphi)
print "Pr/Pz        = " + str(ratio1)
print "(Pr+Pz)/Pphi = " + str(ratio2)
print
print "Ptheo   = %.3e W" %P
print "S       = %.3e W" %S
print "Ptheo/S = %.3e "  %(P/S)
print "Ptot    = %.3e W" %(Pr + Pz + Pphi)
print "S/Ptot  = %.3e " %(S/Ptot)
print
print "*******************"
print

###########################################################
show()

######### Fin #############################################
