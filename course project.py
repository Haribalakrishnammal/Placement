# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 17:44:00 2022

@author: KOMAL
"""
import scipy
import math
import scipy.integrate as scint
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
n = 10 # no of elements
# 0 - COH
# 1 - Cx(OH)2
# 2 - CxO
# 3 - CxO2
# 4 - CxO3
# 5 - CO*
# 6 - cov
# 7 - vac
hash_ = [0 for i in range(n)]
r     = [0 for i in range(n)]
k     = [0 for i in range(n)]
ko    = [0 for i in range(n)]
Ea    = [0 for i in range(n)]
U     = [0 for i in range(n)]
d1    = [0 for i in range(6)]
nj    = [0 for i in range(n)]
alpha_a = [0 for i in range(n)]
alpha_c = [0 for i in range(n)]
R = 8.3144        #J mol-1 K-1
T = 80 + 273.15   #K
F = 96485   #C/mol
g = 3
c_hash =9.55 * (10**-5)   #mol/m2
c_star = 10**-4  #mol/m2
po = 47.41        #kPa at 80C
# po = 19.94        #kPa at 60C
po_ref = 101.325  #kPa
c = 1   #M
c_ref = 1
V = 1.4  #V
S = 100 #m2/g
Nc = 0.29*(15/10) *(10**-9) #mol/m2
MW = 12 

#mol m-2 s-1
ko[1] = 2.35 * (10**-16) * (10**4)
ko[2] = 9.5  * (10**-3)  * (10**4)
ko[3] = 2.19 * (10**-7) * (10**4)
ko[4] = 1.18 * (10**-12) * (10**4)
ko[5] = 4.22 * (10**-18) * (10**4)
ko[6] = 2.35 * (10**-11) * (10**4)
ko[7] = 6.00 * (10**-8) * (10**4)

#kJ mol-1 K-1
Ea[1] = 10
Ea[2] = 110
Ea[3] = 50
Ea[4] = 10
Ea[5] = 10
Ea[6] = 10
Ea[7] = 20

#V
U[1] = 1.00
U[2] = 0.15
U[3] = 0.95
U[4] = 1.00
U[5] = 1.00
U[6] = 1.00
U[7] = 0.57

alpha_a[1] = 0.35
alpha_a[2] = 0.65
alpha_a[3] = 0.48
alpha_a[4] = 0.5
alpha_a[5] = 0.5
alpha_a[6] = 0.5
alpha_a[7] = 0.5

alpha_c[1] = 0
alpha_c[2] = 0
alpha_c[3] = 0.65
alpha_c[4] = 0
alpha_c[5] = 0
alpha_c[6] = 0
alpha_c[7] = 0.5

nj[1] = 1
nj[2] = 2
nj[3] = 2
nj[4] = 4
nj[5] = 0
nj[6] = 4
nj[7] = 2
nr = 0

for i in range (0,n):
    k[i] = ko[i] * math.exp(-Ea[i]/(R))
print ("Ea",Ea, "R",R,"ko",ko)
def model(theta, t, k, alpha_a,T,V,U,po,po_ref,c,c_ref, c_hash, c_star):
    
    theta_cov = theta[0]+3*(theta[1]+theta[2]+theta[3]+theta[4])
    theta_vac = 1 - theta_cov
    
    # r[1] = k[1] * theta_vac * math.exp((((alpha_a[1]*F)/(R*T))*(V-U[1])) - (g*theta_cov))
    # r[2] = k[2] * theta_vac * theta[0] * (po/po_ref) * math.exp(((alpha_a[2]*F)/(R*T))*(V-U[2]))
    # r[3] = k[3] * (((theta_vac**3) * math.exp((((alpha_a[3]*F)/(R*T))*(V-U[3])) - (g*theta_cov))) -
    #                 ( theta[2] * (c/c_ref) * math.exp((((-alpha_c[3]*F)/(R*T))*(V-U[3])) + (g*theta_cov))))
    # r[4] = k[4] * theta_vac * theta[0] * math.exp((((alpha_a[4]*F)/(R*T))*(V-U[4])) - (g*theta_cov))
    # r[5] = -k[5] * theta[2] * theta[5] * math.exp(((-alpha_c[5]*F)/(R*T))*(V-U[5]))
    # r[6] = k[6] * theta_vac * theta[0] * (theta[3]**0.25) * math.exp((((alpha_a[6]*F)/(R*T))*(V-U[6])) - (g*theta_cov))
    # r[7] = k[7] * ((theta[1] * math.exp(((alpha_a[7]*F)/(R*T))*(V-U[7])))  - ( theta[3] * (c/c_ref) * math.exp(((-alpha_c[7]*F)/(R*T))*(V-U[7]) )))
    
    A = theta[0]   #COH
    B = theta[1]  #CxHO2
    C = theta[2]  #CxO
    D = theta[3]#CxO2
    E = theta[4]  #CxO3
    F = theta[5]   #CO*
    r[1] = k[1] * theta_vac * math.exp((((alpha_a[1]*F)/(R*T))*(V-U[1])) - (g*theta_cov))
    r[2] = k[2] * theta_vac * A * (po/po_ref) * math.exp(((alpha_a[2]*F)/(R*T))*(V-U[2]))
    r[3] = k[3] * (((theta_vac**3) * math.exp((((alpha_a[3]*F)/(R*T))*(V-U[3])) - (g*theta_cov))) -
                    ( C * (c/c_ref) * math.exp((((-alpha_c[3]*F)/(R*T))*(V-U[3])) + (g*theta_cov))))
    r[4] = k[4] * theta_vac * A * math.exp((((alpha_a[4]*F)/(R*T))*(V-U[4])) - (g*theta_cov))
    r[5] = -k[5] * C * F * math.exp(((-alpha_c[5]*F)/(R*T))*(V-U[5]))
    r[6] = k[6] * theta_vac * A * (D**0.25) * math.exp((((alpha_a[6]*F)/(R*T))*(V-U[6])) - (g*theta_cov))
    r[7] = k[7] * ((B * math.exp(((alpha_a[7]*F)/(R*T))*(V-U[7])))  - ( D * (c/c_ref) * math.exp(((-alpha_c[7]*F)/(R*T))*(V-U[7]) )))
    
    print ("Vac=",theta_vac, "Cov=",theta_cov)
    
    d1[0] = (r[1]-2*(r[4] + r[6]))/c_hash
    d1[1] = (-r[7])/c_hash
    d1[2] = (r[3]+r[5])/c_hash
    d1[3] = (r[7]-r[5])/c_hash
    d1[4] = (r[4]+r[6])/c_hash
    d1[5] = ((0.5*r[2])+r[5])/c_star
    return d1
#Initial condition
theta_0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
t = np.linspace(0,60,60)
sol = scint.odeint(model,theta_0,t,args=(k, alpha_a,T,V,U,po,po_ref,c,c_ref, c_hash, c_star),)
print (sol)
print (sol[:,2])

# plt.plot(t,sol[:,2],'o')
theta_cov = [0 for i in range (60)]
theta_vac = [0 for i in range (60)]
iCO2 = [0 for i in range (60)]
itotal = [0 for i in range (60)]

for i in range (0,60):
    theta_cov[i] = sol[i,0]+3*(sol[i,1]+sol[i,2]+sol[i,3]+sol[i,4])
    theta_vac[i] = 1 - theta_cov[i]
# plt.plot(t,theta_vac,'g')  
theta = pd.DataFrame({'COH': sol[:,0], 
                          'CxO': sol[:,2],
                         'vac': theta_vac})
datatoexcel = pd.ExcelWriter('theta_data.xlsx')
theta.to_excel(datatoexcel)
datatoexcel.save()
#Total current
for i in range(0,60):
    nr = 0
    theta_cov = sol[i,0]+3*(sol[i,1]+sol[i,2]+sol[i,3]+sol[i,4])
    theta_vac = 1 - theta_cov
    r[1] = k[1] * theta_vac * math.exp((((alpha_a[1]*F)/(R*T))*(V-U[1])) - (g*theta_cov))
    r[2] = k[2] * theta_vac * sol[i,0] * (po/po_ref) * math.exp(((alpha_a[2]*F)/(R*T))*(V-U[2]))
    r[3] = k[3] * (((theta_vac**3) * math.exp((((alpha_a[3]*F)/(R*T))*(V-U[3])) - (g*theta_cov))) -
                    ( sol[i,2] * (c/c_ref) * math.exp((((-alpha_c[3]*F)/(R*T))*(V-U[3])) + (g*theta_cov))))
    r[4] = k[4] * theta_vac * sol[i,0] * math.exp((((alpha_a[4]*F)/(R*T))*(V-U[4])) - (g*theta_cov))
    r[5] = -k[5] * sol[i,2] * sol[i,5] * math.exp(((-alpha_c[5]*F)/(R*T))*(V-U[5]))
    r[6] = k[6] * theta_vac * sol[i,0] * (sol[i,3]**0.25) * math.exp((((alpha_a[6]*F)/(R*T))*(V-U[6])) - (g*theta_cov))
    r[7] = k[7] * ((sol[i,1] * math.exp(((alpha_a[7]*F)/(R*T))*(V-U[7])))  - ( sol[i,3] * (c/c_ref) * math.exp(((-alpha_c[7]*F)/(R*T))*(V-U[7]) )))
    for j in range (0,8):
        nr = nr + (nj[j]*r[j])
    itotal[i] = S*Nc*MW*F*nr
    iCO2[i] = S*Nc*MW*F*nj[2]*r[2]
    
# print ("iCO2=",iCO2,"itotal=",itotal)
print ("iCO2=",iCO2)
plt.plot(t,iCO2,'o')

iCO2 = pd.DataFrame({'iCO2 ':iCO2})
datatoexcel = pd.ExcelWriter('CO2_1.4V.xlsx')
iCO2.to_excel(datatoexcel)
datatoexcel.save()

time = [0.0953795*60,0.113732*60,0.135615*60,0.171609*60,0.204629*60,0.274583*60,0.391007*60,0.494784*60,0.626342*60,0.84078*60]
CO2_current= [0.13895,0.156536,0.176349,0.191029,0.215207,0.242569,0.262894,0.284777,0.296543,0.321309]

plt.plot(time,CO2_current,'g')










