#! /usr/bin/env python
# -*- coding: utf-8 -*-

######### Modules #########################################
from pylab import *
from scipy.special import j0,jn

###########################################################
# Champ d'après 
# A.April, Opt. Lett. 331563 (2008)
k = 2.0*pi
a = 1.0/k
wo = (sqrt(2.0)/k)*sqrt(sqrt(1.0+(k*a)**2)-1.0);

print
print "*******************"
print "ka = " + str(k*a)
print "wo = " + str(wo)
print "*******************"
print

###########################################################
lim = 10.0
if (k*a <= 100): lim = 6.0
if (k*a <= 70): lim = 4.0
if (k*a <= 30): lim = 3.0
if (k*a <= 10): lim = 2.0
r = linspace(0,lim,1000)
z = 0

R = sqrt(r**2 + (z + a*1j)**2)
theta = arccos((z+a*1j)/R)
kR = k*R

###########################################################
# Un gaussienne
figure(1)
Ex = jn(0+0.5,kR)*sqrt(pi/kR)
wx = abs(Ex)**2
maxi = max(wx)
rscl = sqrt(R**2 - (z + a*1j)**2)
title(r"$w_0 = %f \lambda_0,\quad ka = %f$" %(wo,k*a),size = 14)
plot(rscl,wx/maxi,label=r"$|E_x|^2$")
plot(rscl,exp(-2.0*r**2/wo**2),label=r"$\exp(-2r^2/w_0^2)$")
ylim(0,1)
xlim(0,lim)
grid(True)
legend(loc="best")

###########################################################
# Champ électrique du TM01 [Eqs. (6) de A.April, Opt. Lett. 331563 (2008)]
Er = 0.5*jn(2+0.5,kR)*sin(2.0*theta)*sqrt(pi/kR)
Ez = 2.0/3.0*(jn(0+0.5,kR) + jn(2+0.5,kR)*0.25*(1.0+3.0*cos(2.0*theta)) )*sqrt(pi/kR)

###########################################################
# Figure
figure()
wr = abs(Er)**2
wz = abs(Ez)**2
maxi = max(wr + wz)
rscl = sqrt(R**2 - (z + a*1j)**2)
title(r"$w_0 = %f \lambda_0,\quad ka = %f$" %(wo,k*a),size = 14)
plot(rscl,wr/maxi,label=r"$|E_r|^2$")
plot(rscl,wz/maxi,label=r"$|E_z|^2$")
plot(rscl,(wr+wz)/maxi,label=r"$|E_r|^2 + |E_z|^2$")
ylim(0,1)
xlim(0,lim)
grid(True)
legend(loc="best")

###########################################################
show()

######### Fin #############################################
