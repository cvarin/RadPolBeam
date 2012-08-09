#!/usr/bin/python2
# -*- coding: utf-8 -*-

from pylab import *
from scipy import constants as cst

lambda0 = 0.8e-6
k = 2.0*pi/lambda0
omega = k*cst.c
#a = 278.581/k
a = linspace(0,350)/k

wo = (sqrt(2.0)/k)*sqrt(sqrt(1.0+(k*a)**2)-1.0);
eta0 = 120.0*pi

E0 = omega*cst.c*cst.m_e/cst.e

Pana = pi/(4.0*eta0*k**2)*E0**2*exp(-2.0*k*a)\
    *(2.0*k*a*sinh(2.0*k*a) - cosh(2.0*k*a) + 1.0 - 2.0*(k*a)**2)

plot(wo/1.0e-6,Pana/1.0e12,"r")


a = linspace(0,2500)/k
Pparaxial = 2.0*(pi*a/(2.0*k))*E0**2/(4.0*eta0)
wo = (sqrt(2.0)/k)*sqrt(sqrt(1.0+(k*a)**2)-1.0);
plot(wo/1.0e-6,Pparaxial/1.0e12,"xb")
xscale("log")
yscale("log")
xlim(0,10)
#ylim(1,)
show()
