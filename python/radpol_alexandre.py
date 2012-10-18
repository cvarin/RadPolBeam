#! /usr/bin/python
# -*- coding: utf-8 -*-

######### Modules #########################################
from pylab import *
from scipy.special import j0,jn

eta0 = 120.0*pi
mu0 = 12.566370614e-7
eps0 = 8.854187817e-12
eta0 = 120.0*pi
c    = 2.99792458e8

tsize = 16

###########################################################
# Champ d'après 
# A.April, Opt. Lett. 331563 (2008)
lambda0 = 0.8e-6
k = 2.0*pi/lambda0
#a = 278.581/k
a = 350.0/k
wo = (sqrt(2.0)/k)*sqrt(sqrt(1.0+(k*a)**2)-1.0);

print
print "*******************"
print "lambda0    = " + str(lambda0)
print "ka         = " + str(k*a)
print "wo/lambda0 = " + str(wo/lambda0)
print "wo         = " + str(wo)
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

R = sqrt(r**2 + (z + a*1j)**2)
theta = arccos((z+a*1j)/R)
kR = k*R

###########################################################
# Champ électrique du TM01 [Eqs. (6) de A.April, Opt. Lett. 331563 (2008)]
# E0 est défini à partir de la puissance 
# [Eq. (12) de A. April, J. Opt. Soc. Am. A 27, 76-81 (2010)]
#E0 = 1.0e13
#H0 = E0/eta0
#H0p = 0.5*H0*k**3*wo**3
#E0p = eta0*H0p

###########################################################
# Puissance en régime paraxial : Eq. (13) de
# A. April, J. Opt. Soc. Am. A 27, 76-81 (2010)
P = 1.0e12
E0 = sqrt(P/2.0/(pi*a/(2.0*k))*(4.0*eta0))
E0p = E0
H0p = E0p/eta0

Hphi = -H0p*sqrt(pi/kR/2.0)*jn(1+0.5,kR)*sin(theta)*exp(-k*a)
Er = 0.5*1j*E0p*jn(2+0.5,kR)*sin(2.0*theta)*sqrt(pi/kR/2.0)*exp(-k*a)
Ez = 2.0/3.0*1j*E0p*(jn(0+0.5,kR) + jn(2+0.5,kR)*0.25*(1.0+3.0*cos(2.0*theta)))\
       *sqrt(pi/kR/2.0)*exp(-k*a)

###########################################################
# Puissance
# Ref: A. April, J. Opt. Soc. Am. A 27, 76-81 (2010)
###########################################################
# Eq. (12)
Pana = pi/(4.0*eta0*k**2)*E0**2*exp(-2.0*k*a)\
    *(2.0*k*a*sinh(2.0*k*a) - cosh(2.0*k*a) + 1.0 - 2.0*(k*a)**2)
    
###########################################################
# Eq. (13)
Pparaxial = 2.0*(pi*a/(2.0*k))*E0**2/(4.0*eta0)

###########################################################
# Eq. (3)
Pnum = 0
for i in range(0,np):
     Pnum += real(Er[i]*Hphi[i].conjugate())*r[i]
Pnum*=0.5*2.0*pi*dr

###########################################################
# Normalisation du champ
K00 = sqrt(Pana/Pnum)
Hphi *= K00
Er *= K00
Ez *= K00

###########################################################
# Tracé des composantes
fig1 = figure(figsize=(16.0,8.5))
fig1.canvas.set_window_title(unicode("Distribution du champ électrique (non-paraxial)","utf-8"))
plot(r/lambda0,real(Er),label="Er")
plot(r/lambda0,real(-1j*Ez),label="Ez")
#plot(r/lambda0,real(Er)/real(Er).max(),label="Er")
#plot(r/lambda0,real(-1j*Ez)/real(-1j*Ez).max(),label="Ez")
xlabel(r"$r/\lambda$",size=tsize)
ylabel(unicode("Champ électrique normalisé","utf-8"),size=tsize)
#ylim(-0.2,1.0001)
xlim(0,lim/lambda0)
grid(True)
legend(loc="best")

###########################################################
# Densité d'énergie
wr = 0.25*eps0*abs(Er)**2
wz = 0.25*eps0*abs(Ez)**2
wphi = 0.25*mu0*abs(Hphi)**2

###########################################################
# Figure
fig2 = figure(figsize=(16.0,8.5))
fig2.canvas.set_window_title(unicode("Densité d'énergie électrique (non-paraxial)","utf-8"))
maxi = max(wr + wz)
rscl = sqrt(R**2 - (z + a*1j)**2)/lambda0
title(r"$w_0 = %f \lambda_0,\quad ka = %f$" %(wo/lambda0,k*a),size = tsize)
plot(rscl,wr/maxi,label=r"$|E_r|^2$")
plot(rscl,wz/maxi,label=r"$|E_z|^2$")
plot(rscl,(wr+wz)/maxi,label=r"$|E_r|^2 + |E_z|^2$")
#plot(rscl,wphi/maxi,label=r"$|H_\phi|^2$")
xlabel(r"$r/\lambda$",size=tsize)
ylabel(unicode("Champ au carré","utf-8"),size=tsize)
ylim(0,1)
xlim(0,lim/lambda0)
grid(True)
legend(loc="best")

###########################################################
# Puissance
Pnum = 0
for i in range(0,np):
     Pnum += real(Er[i]*Hphi[i].conjugate())*r[i]
Pnum*=0.5*2.0*pi*dr
print
print "*******************"
print "E0  = %.3e V/m" %E0
print "K00 = %.3e" %(K00)
print
print "Pana  = %.3e W" %Pana
print "Pnum  = %.3e W" %Pnum
print "Ppara = %.3e W" %Pparaxial
print
print "Pana/Ppara = %.3e" %(Pana/Pparaxial)
print "Pnum/Pana  = %.3e" %(Pnum/Pana)
print

###########################################################
# Comparaison de l'énergie dans chaque composante.
Pr = 0
Pz = 0
Pphi = 0
for i in range(0,np):
     Pr += wr[i]*r[i]
     Pz += wz[i]*r[i]
     Pphi += wphi[i]*r[i]
Pr*=2.0*pi*dr*c
Pz*=2.0*pi*dr*c
Pphi*=2.0*pi*dr*c
Ptot = Pr + Pz + Pphi
ratio1 = real(Pr/Pz)
ratio2 = real((Pr+Pz)/Pphi)

print "*******************"
print "Pr   = %.3e " %real(Pr)
print "Pz   = %.3e " %real(Pz)
print "Pphi = %.3e " %real(Pphi)
print "-------------------"
print "Ptot = %.3e " %real(Ptot)
print "-------------------"
print
print "Pr/Pz        = " + str(ratio1)
print "(Pr+Pz)/Pphi = " + str(ratio2)
print

###########################################################
print "*******************"
print
show()

######### Fin #############################################
