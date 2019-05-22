#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 17:06:13 2019

@author: constanzamunoz
"""
import numpy as np
import matplotlib.pyplot as plt


#All the data comes from an adapted table from Astrophysical Quantities 4th ed.
M_sun = 1.9891* 1e30 #kg
R_sun = 0.7 * 1e9 #m
sun_density = 1408.0 #kg /m**3

"""
Classical cepheids Data 
"""
Cph_p = np.array([0.4, 0.6, 0.8, 1.0, 1.2 , 1.4, 1.6 , 1.8]) #log P [days]
Cph_period = 10**Cph_p
Cph_r = np.array([1.41, 1.55, 1.70, 1.85, 2.00, 2.15, 2.29, 2.44]) #log R/R_sun
Cph_radius = (10**(Cph_r))* R_sun
Cph_m = np.array([0.54, 0.61, 0.68, 0.75, 0.82, 0.89, 0.97, 1.04]) #log M/M_sun
Cph_mass = (10**(Cph_m ))* M_sun
Cph_Mv = np.array([-2.4, -2.9, -3.5, -4.1, -4.7, -5.3, -5.8, -6.4]) #Magnitude V
Cph_L = np.array([2.81, 3.05, 3.30, 3.55, 3.80, 4.06, 4.32, 4.58]) #log L/L_sun
Cph_g = np.array([2.2, 1.9, 1.7, 1.5, 1.3, 1.0, 0.8, 0.6]) #Log g
Cph_BV = np.array([0.49, 0.57, 0.66, 0.75, 0.84, 0.93, 1.01, 1.10]) #B-V

"""
dSct/Dwarf Cepheids Data
"""
d_p = np.array([-1.4, -1.0, -0.7]) #log P [days]
d_period = 10**d_p
d_r = np.array([0.07, 0.37, 0.59]) #log R/R_sun
d_radius = (10**(d_r)) * R_sun
d_m= np.array([0.21, 0.25, 0.31]) #log M/M_sun
d_mass = (10**(d_m)) * (M_sun)
d_L = np.array([0.90, 1.22, 1.58]) #log L/L_sun
d_g = np.array([4.5, 3.9, 3.6]) #Log g

"""
RR Lyrae stars Data
"""
rr_p = np.array([-0.5, -0.3, -0.1]) #log P [days]
rr_period = 10**rr_p
rr_r = np.array([0.54, 0.67, 0.76]) #log R/R_sun
rr_radius = (10**(rr_r)) * R_sun
rr_m = np.array([-0.26, -0.26, -0.26]) #log M/M_sun
rr_mass = (10**(rr_m)) * (M_sun)
rr_L = np.array([1.47, 1.57, 1.59]) #log L/L_sun
rr_g = np.array([3.1, 2.8, 2.7]) #Log g

"""
P2. a) obtain Q dor each subclasses
"""

def Q(P, R, M):  
    """
    Q: pulsation constant
    from mass M [KG], period P [days] and radius R[m] return Q 
    """
    Q = P * ((3*M)/ (4 *np.pi * (R**3) * sun_density))**0.5
    return Q

Q_Cph = np.array([])
c_cph = 0
for j in Cph_period:
    value = Q(j, Cph_radius[c_cph], Cph_mass[c_cph])
    Q_Cph = np.append(Q_Cph, value)
    c_cph = c_cph + 1
    
Q_d = np.array([])
c_d = 0
for j in d_period:
    value = Q(j, d_radius[c_d], d_mass[c_d])
    Q_d = np.append(Q_d, value)
    c_d = c_d + 1


Q_rr = np.array([])
c_rr = 0
for j in rr_period:
    value = Q(j, rr_radius[c_rr], rr_mass[c_rr])
    Q_rr = np.append(Q_rr, value)
    c_rr = c_rr + 1

"""
P2. b) P-L, period-radius and period-surface garvity relation 
"""

plt.figure(figsize=(7,5))
plt.title("Period-Luminosity relation (PL) of pulsating variable stars", size= 12)
plt.plot(Cph_p, Cph_L, '-o', color='blue', label ="Classical Cepheids")
plt.plot(d_p, d_L, '-o', color='green', label= r'$\delta$ Sct/Dwarf Cepheids')
plt.plot(rr_p, rr_L, '-o', color='red', label= "RR Lyrae stars")
plt.ylabel(r'log $L /L_{\odot}$', size= 12)
plt.xlabel('log P [days]', size= 12)
plt.legend(prop={'size': 12})
plt.savefig("P-L.png")

plt.figure(figsize=(7,5))
plt.title("Period-Radius relation of pulsating variable stars", size= 12)
plt.plot(Cph_p, Cph_r, '-o', color='blue', label ="Classical Cepheids")
plt.plot(d_p, d_r, '-o', color='green', label= r'$\delta$ Sct/Dwarf Cepheids')
plt.plot(rr_p, rr_r, '-o', color='red', label= "RR Lyrae stars")
plt.ylabel(r'log $R /R_{\odot}$', size= 12)
plt.xlabel('log P [days]', size= 12)
plt.legend(prop={'size': 12})
plt.savefig("P-R.png")

plt.figure(figsize=(7,5))
plt.title("Period-Surface Gravity relation of pulsating variable stars", size= 12)
plt.plot(Cph_p, Cph_g, '-o', color='blue', label ="Classical Cepheids")
plt.plot(d_p, d_g, '-o', color='green', label= r'$\delta$ Sct/Dwarf Cepheids')
plt.plot(rr_p, rr_g, '-o', color='red', label= "RR Lyrae stars")
plt.ylabel('log g', size= 12)
plt.xlabel('log P [days]', size= 12)
plt.legend(prop={'size': 12})
plt.savefig("P-G.png")
    

"""
P2. c) PLC, period-luminosity-color relation for classical cepheids
"""

m_1 = -0.7
m_2 = -0.4
m_3 = 0
mean_Mv_1 = -2.83  -(3.57*Cph_p) +(2.92*Cph_BV) -(0.35*m_1) #<Mv> for metallicity 1
mean_Mv_2 = -2.83  -(3.57*Cph_p) +(2.92*Cph_BV) -(0.35*m_2) #<Mv> for metallicity 2
mean_Mv_3 = -2.83  -(3.57*Cph_p) +(2.92*Cph_BV) -(0.35*m_3) #<Mv> for metallicity 3
plt.figure(figsize=(7,5))
plt.title("PLC: period-luminosity-color relation for classical cepheids ", size= 12)
plt.plot(Cph_p, mean_Mv_1, '-o', color='navy', label =r"$[\frac{Fe}{H}]$ = -0.7")
plt.plot(Cph_p, mean_Mv_2, '-o', color='blue', label =r"$[\frac{Fe}{H}]$ = -0.4")
plt.plot(Cph_p, mean_Mv_3, '-o', color='aqua', label =r"$[\frac{Fe}{H}]$ = 0")
plt.ylabel(r'<$M_{V}$> [mag]', size= 12)
plt.xlabel('log P [days]', size= 12)
plt.legend(prop={'size': 12})
plt.gca().invert_yaxis()
plt.savefig("PLC.png")
    
