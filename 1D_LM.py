#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 13:11:32 2022

@author: Casper Vindahl Jensen
Department of Chemistry
University of Copenhagen

Note about script:
    
This script will take the textfile created by the script "get_E_dipole.py" and solve the vibrational schrödinger equation
from the potential and calculate the oscillator strength of the OH stretch.

This is done by an interpolation of q, E and DM, a renormalization of the q-coordinate to [-1 ; +1].
A Hamiltonian is constructed and a basis of associated Legendre polynomials with m=1 is used.
The Hamiltonian is diagonalised and a wavefunction is constructed from the eigenvalues. 


To use this script you need textfile with displacement q, change in electronic energy E, and the dipole moment in x,y,z
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy import integrate
#%%
data = pd.read_csv("TBHP_1D_wB.txt", header=0, sep="\t", decimal=".", engine="python") #specify the txt file with q,E,mu_xyz
#%%
def constants():
    c = 299792458 #m/s
    h = 6.62607015*10**(-34) #joule*sec 
    h_bar = h/(2*np.pi)
    amu = 1.66053907*10**(-27)
    return c, h, h_bar, amu

[c, h, h_bar, amu] = constants()
    
def makeredmass(m1,m2): 
    realmass_list = pd.Series([1.007825, 2.014102, 12, 13.00335, 14.00307, 15.99491]) #some masses (H, D, C, C13, N, O)
    a = (realmass_list).sub(m1).abs().idxmin() #find mass closest to input
    m1 = realmass_list[a]
    a = (realmass_list).sub(m2).abs().idxmin()
    m2 = realmass_list[a]
    
    redmass = (m1*m2)/(m1+m2)
    return redmass

def makebasis(x, l): #makes Legendre pol
    l+=1 #this is because the range function stops at n-1
    N = np.sqrt(2*(2*l+1)/(2*l*(l+1)* q_range)) #normalization
    P0 = -np.sqrt(1-x**2)
    P1 = -3*x*np.sqrt(1-x**2)
                                    #first two legendre pol are hardcoded. the rest is from recursive formula
    if l == 1:
        return P0*N
    if l == 2:
        return P1*N
    else:
        for a in range(2,l):
            P = ( (2*a + 1)*x*P1 - (a+1)*P0 ) / (a)
            
            P0 = P1*1
            P1 = P*1
        return P*N      #returns the l.th basisfunction



def makedifbasis(x, l): #makes differentiated Legendre pol
    l+=1
    N = np.sqrt((2*l+1)/(2*l*(l+1))) * np.sqrt(2/q_range) * 2/q_range #normalization
    P0 = -np.sqrt(1-x**2)
    P1 = -3*x*np.sqrt(1-x**2)

    if l == 1:
        P_dif = (x[1:-1] / np.sqrt(1-x[1:-1]**2)) 
    
    elif l == 2:
        P_dif =  (6*x[1:-1]**2 - 3)/np.sqrt(1-x[1:-1]**2) 
    
    else:
        for a in range(2,l):
            P = ( (2*a + 1)*x*P1 - (a+1)*P0 ) / (a)
            P0 = P1*1
            P1 = P*1
        P_dif = (l*x[1:-1]*P[1:-1]-(l+1)*P0[1:-1]) / ((x[1:-1])**2 - 1) 


    P_dif = np.insert(P_dif, 0, P_dif[0], axis=0)
    P_dif = np.append(P_dif, P_dif[-1], axis=None)
    return P_dif*N

def Hamilton(N_basis): #makes Hamiltonian
    H = np.zeros((N_basis, N_basis),dtype=float) #hamiltonian 
    V = np.zeros((N_basis, N_basis),dtype=float) #potential energy operator
    T = np.zeros((N_basis, N_basis),dtype=float) #Kinetic energy operator
    
    basis = np.zeros((len(q),N_basis))
    basis_dif = np.zeros((len(q),N_basis))    
    
    for i in range(0,N_basis):
        basis[:,i] = makebasis(x,i)
        basis_dif[:,i] = makedifbasis(x,i)
    
    for i in range(0,N_basis):
        for j in range(i,N_basis):
            V[i,j] = float(integrate.simps(basis[:,i]*basis[:,j]*pot, q))
            V[j,i] = V[i,j]
            
            T[i,j] = ((h_bar**2)/(2*red_mass*h*c*100)) *10**20 *float(integrate.simps(basis_dif[:,i]*basis_dif[:,j], q))
            T[j,i] = T[i,j]
            
    H = V + T
    return H

def makeWF(): #makes Wavefunction
    basis = np.zeros((len(q),Nbasis))  
    wavefunc = np.zeros((len(q),Nvec)) 
    
    for i in range(0,Nbasis):
        basis[:,i] = makebasis(x,i)
        
    for i in range(0,Nvec):
        for j in range(0,Nbasis):
            wavefunc[:,i] += eigvec[j,i]*basis[:,j]
            
    return wavefunc
#%%
#Choose some parameters
m1 = 1 #au
m2 = 16 #au
Nbasis = 50 #basisfunctions
Nvec = 7    #number of eigenvectors wanted. limited by Nbasis.
delta_q = 0.00002 #interpolation of q
#%%
red_mass = makeredmass(m1,m2)
red_mass = red_mass * amu

q_range = data["q"].max() - data["q"].min()
q = np.linspace(data["q"].min(), data["q"].max(), round(q_range/delta_q))

cs = CubicSpline(data["q"], data["E"])
pot = cs(q) #interpolates the electronic energy
#%% interpolates the dipole moments
cs1 = CubicSpline(data["q"], data["mux"])
mux = cs1(q)
cs2 = CubicSpline(data["q"], data["muy"])
muy = cs2(q)
cs3 = CubicSpline(data["q"], data["muz"])
muz = cs3(q)
#%%
x = np.linspace(-1.0,1.0,len(q)) #renormalization of displacement
#%% check if renormalization from q to x is consistent with each other. small variation of 3% is seen in our model
psi = makebasis(x,1)
area_new = float(integrate.simps(psi**2, x))
area_old = float(integrate.simps(psi**2, q))

print(area_new)
print(area_old)
#%%
H = Hamilton(Nbasis) #construct Hamiltonian
eig = np.linalg.eig(H) #pull eigenvalues
sorting = np.argsort(eig[0]) #sorting of eigval from lowest to highest energy
eigE = eig[0][sorting]
eigvec = eig[1][:,sorting] 

WF = makeWF() #make wavefunction from eigenvectors

#%%
plt.plot(q,mux)
plt.plot(q,muy)
plt.plot(q,muz)
#%% calculate oscillator strength of fundamental

v_i = 0
v_f = 1

osc = (4.702*10**(-7)*(eigE[1]-eigE[v_i])*
             (abs(integrate.simps(WF[:,v_i]*mux*WF[:,v_f],q))**2+
              abs(integrate.simps(WF[:,v_i]*muy*WF[:,v_f],q))**2+
              abs(integrate.simps(WF[:,v_i]*muz*WF[:,v_f],q))**2))
#%%
print("transition energy = ", eigE[v_f]-eigE[v_i], "transition intensity = ", osc)

























