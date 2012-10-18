#! /usr/bin/env python
# -*- coding: utf-8 -*-

######### Modules #########################################
from pylab import *
from scipy.special import j0,jn

mu0 = 12.566370614e-7
eps0 = 8.854187817e-12
eta0 = 120.0*pi
c    = 2.99792458e8

tsize = 16

###########################################################
# Champ d'après 
lambda0 = 0.8e-6
k = 2.0*pi/lambda0
omega = k*c
a = 350.0/k
wo = (sqrt(2.0)/k)*sqrt(sqrt(1.0+(k*a)**2)-1.0);
zR = 0.5*k*wo**2

print
print "*******************"
print "lambda0    = " + str(lambda0)
print "ka         = " + str(k*a)
print "wo/lambda0 = " + str(wo/lambda0)
print "zR/lambda0 = " + str(zR/lambda0)
print "*******************"
print

###########################################################
lim = 20.0
if (k*a <= 400): lim = 14.0
if (k*a <= 300): lim = 10.0
if (k*a <= 100): lim = 6.0
if (k*a <= 70): lim = 4.0
if (k*a <= 30): lim = 3.0
if (k*a <= 10): lim = 2.0
lim *= lambda0
np = 1000
r = linspace(0,lim,np)
dr = lim/np
z = 0
t = 0

w = wo*sqrt(1+z**2/zR**2)
R = z + zR**2/z
phi_g = arctan(z/zR)
phi_c = k*r**2/(2.0*R)
psi = omega*t - k*z + 2.0*phi_g - phi_c

###########################################################
# Champ électrique du TM01 paraxial
P = 1.0e12
E0 = sqrt(4.0*eta0*P/pi*exp(-1.0)/wo**2)
phi_0 = 0.0
Er = E0*exp(0.5)*(sqrt(2.0)*r/wo)*(wo/w)**2*exp(-r**2/w**2)*exp(1j*(psi - phi_0))
#*cos(psi - phi_0)
#phi_0 = -pi/2.0
Ez = E0*exp(0.5)*(2.0*sqrt(2.0)/(k*wo))*(wo/w)**2*exp(-r**2/w**2)\
      *((1-r**2/w**2)*exp(1j*(psi - phi_0)) - k*r**2/(2.0*R)*exp(1j*(psi - phi_0)))
      #*((1-r**2/w**2)*sin(psi - phi_0) - k*r**2/(2.0*R)*cos(psi - phi_0))

###########################################################
# Tracé des composantes
fig1 = figure(figsize=(16.0,8.5))
fig1.canvas.set_window_title(unicode("Distribution du champ électrique (paraxial)","utf-8"))
plot(r/lambda0,Er,label="Er")
plot(r/lambda0,Ez,label="Ez")
#plot(r/lambda0,Er/max(Er),label="Er")
#plot(r/lambda0,Ez/max(Ez),label="Ez")
xlabel(r"$r/\lambda$",size=tsize)
ylabel(unicode("Champ électrique normalisé","utf-8"),size=tsize)
#ylim(-0.2,1.001)
xlim(0,lim/lambda0)
grid(True)
legend(loc="best")

print
print "*******************"
print "P = %.3e W" %P
print "E0 = %.3e V/m" %E0
print

###########################################################
# Densité d'énergie
wr = 0.25*eps0*Er*Er.conjugate()
wz = 0.25*eps0*Ez*Ez.conjugate()
#wr = 0.25*eps0*Er**2
#wz = 0.25*eps0*Ez**2

###########################################################
# Figure
fig2 = figure(figsize=(16.0,8.5))
fig2.canvas.set_window_title(unicode("Densité d'énergie électrique (paraxial)","utf-8"))
maxi = max(wr + wz)
rscl = r/lambda0
title(r"$w_0 = %f \lambda_0,\quad ka = %f$" %(wo/lambda0,k*a),size = tsize)
plot(rscl,wr/maxi,label=r"$|E_r|^2$")
plot(rscl,wz/maxi,label=r"$|E_z|^2$")
plot(rscl,(wr+wz)/maxi,label=r"$|E_r|^2 + |E_z|^2$")
xlabel(r"$r/\lambda$",size=tsize)
ylabel(unicode("Champ au carré","utf-8"),size=tsize)
ylim(0,1)
xlim(0,lim/lambda0)
grid(True)
legend(loc="best")

###########################################################
# Puissance
Ptheo = E0**2/(2.0*eta0)*pi/2.0*wo**2*exp(1.0)

P = 0.0
for i in range(0,np):
     P += Er[i]**2*r[i]
P*=pi/eta0*dr

print
print "*******************"
print "P   = %.3e W" %P
print "Ptheo = %.3e W" %Ptheo
print

###########################################################
show()

######### Fin #############################################
