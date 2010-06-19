#! /usr/bin/env python
# -*- coding: utf-8 -*-

######### Modules #########################################
from pylab import *
from scipy.special import j0,jn

###########################################################
# Champ d'après 
# C. J. R. Sheppard and S. Saghaﬁ, Opt. Lett. 24, 1543 (1999)
wo = 0.5
lo = 1.0
k = 2.0*pi/lo
zR = 0.5*k*wo**2
zo = wo*sqrt(1+(0.5*k*wo)**2)

###########################################################
print
print "*******************"
print "wo = " + str(wo)
print "kzo = " + str(k*zo)
print "*******************"
print

###########################################################
lim = 2.0
x = linspace(-lim,lim,1000)
y = 0
z = 0
R = sqrt(x**2 + y**2 + (z - zo*1j)**2)

###########################################################
kR = k*R
f = (jn(0+0.5,kR) + jn(2+0.5,kR))*sqrt(pi/kR)
g = (jn(0+0.5,kR) - 0.5*jn(2+0.5,kR))*sqrt(pi/kR)

Ex = x*(f - g)/R**2*(z - zo*1j)
Ey = y*(f - g)/R**2*(z - zo*1j)
Ez = ( f*(z - zo*1j)**2 + g*(x**2 + y**2) )/R**2

###########################################################
# Densité d'énergie
wx = abs(Ex)**2
wy = abs(Ey)**2
wz = abs(Ez)**2

###########################################################
# Graphiques
maxi = max(array([max(wx),max(wy),max(wz),max(wx + wy + wz)]))
title(r"$w_0 = %f \lambda_0,\quad kz_0 = %f$" %(wo,k*zo),size = 14)
plot(x,wx/maxi,label=r"$|E_x|^2$")
plot(x,wz/maxi,label=r"$|E_z|^2$")
plot(x,(wx+wz)/maxi,label=r"$|E_x|^2+|E_z|^2$")
ylim(0,1)
xlim(0,lim)
legend(loc="best")
grid(True)

###########################################################
show()

######### Fin #############################################