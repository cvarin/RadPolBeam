#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pylab import *
import sys

kzo = float(sys.argv[1])
k = 2.0*pi
zo = 1.0/k
wo = (sqrt(2.0)/k)*sqrt(sqrt(1.0+(kzo)**2)-1.0);

###############################################################################
# Lecture du fichier
d = loadtxt("intensite_kzo_%.3f.dat" %kzo)
r = d[:,0]
wr = d[:,1]
wz = d[:,2]
wtot = d[:,3]
maxi = max(wtot)

###############################################################################
# Graphique
figure(figsize=(16.0,8.5))
txtsize = 14
title(r"$w_0 = %f \lambda_0,\quad kz_0 = %f$" %(wo,kzo),size = txtsize + 4)
plot(r,wr/maxi,label=r"$|E_r|^2$")
plot(r,wz/maxi,label=r"$|E_z|^2$")
plot(r,wtot/maxi,label=r"$|E_r|^2 + |E_z|^2$")
ylabel(unicode("Densité d'énergie électrique (u. a.)"),size=txtsize)
xlabel(unicode("Position transversale en multiples de longueur d'onde"),size=txtsize)
legend(loc="best")
grid(True)
show()

######### End of file #########################################################