#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pylab import *
import sys

kzo = float(sys.argv[1])

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
plot(r,wr/maxi,label=r"$|E_r|^2$")
plot(r,wz/maxi,label=r"$|E_z|^2$")
plot(r,wtot/maxi,label=r"$|E_r|^2 + |E_z|^2$")
ylabel(unicode("Densité d'énergie électrique (u. a.)"),size=txtsize)
xlabel(unicode("Position radiale en multiples de longueur d'onde"),size=txtsize)
legend(loc="best")
grid(True)
show()

######### End of file #########################################################