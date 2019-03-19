################################################################################
#
# TL_tlerence.py
#
# Purpose: Investigating how changes in pyhsical parameters effect impedance
# of microstrip line
#
# Written by: Bailey Brookes for the 4B25 RF Systems Coursework
#
################################################################################

import numpy as np
import matplotlib.pyplot as plt
import math

# Plot parameters
plt.rc('xtick',labelsize=21)
plt.rc('ytick',labelsize=21)

# Matching parameters calcuated by AWR/Specified by user, dimensions in mm
H_calc  = 0.22
W_calc  = 0.4071
T_calc  = 0.035
Er_calc = 4.5
Z0_calc = 48.4

# Anaylysis parameters
num_points = 100

# Fucntion defintions ##########################################################
def Impedance(H,W,T,Er):
    # All equations from https://www.allaboutcircuits.com/tools/microstrip-impedance-calculator/
    Weff = W + (T/np.pi)*np.log((4*np.e)/np.sqrt(math.pow(T/H,2) + math.pow(T/(W*np.pi + 1.1*T*np.pi),2))) * (Er+1)/(2*Er)
    X1 = 4 * ((14*Er + 8)/(11*Er)) * (H/Weff)
    X2 = math.sqrt(16 * math.pow(H/Weff, 2) * math.pow((14*Er + 8)/(11*Er),2) + ((Er + 1)/(2*Er)) * math.pow(np.pi, 2))
    Z0 = (eta/(2*np.pi*np.sqrt(2)*np.sqrt(Er+1))) * np.log(1 + 4*(H/Weff)*(X1 + X2))
    return Z0

def percentDiff(array, calc_value):
    min = array.min()
    max = array.max()
    max_change = ((max - calc_value) / calc_value) * 100
    min_change = ((min - calc_value) / calc_value) * 100
    results = [max_change, min_change]
    return results

# Main #########################################################################
Hs  = np.linspace(H_calc*0.9 , H_calc*1.1 , num = num_points) # Height of track above ground plane in mm (0.635)
Ws  = np.linspace(W_calc*0.9 , W_calc*1.1 , num = num_points) # Width of track in mm (1.1907)
Ts  = np.linspace(T_calc*0.9 , T_calc*1.1 , num = num_points) # Thickeness of track in mm (0.01)
Ers  = np.linspace(Er_calc*0.9, Er_calc*1.1, num = num_points)# Dieletric Constant of subtract (4.5)
eta = 376.730313                                              # Impedance of free space
tol = np.linspace(-10,10, num = num_points)                   # Define tolerence to plot against

# Empty imdedance arrays to store results
Z0_H  = np.array([])
Z0_W  = np.array([])
Z0_T  = np.array([])
Z0_Er = np.array([])

# Cacluate impedances
for H in Hs:
    W = W_calc
    T = T_calc
    Er = Er_calc
    Z0_H = np.append(Z0_H, Impedance(H,W,T,Er))

for W in Ws:
    H = H_calc
    T = T_calc
    Er = Er_calc
    Z0_W = np.append(Z0_W, Impedance(H,W,T,Er))

for T in Ts:
    H = H_calc
    W = W_calc
    Er = Er_calc
    Z0_T = np.append(Z0_T, Impedance(H,W,T,Er))

for Er in Ers:
    H = H_calc
    W = W_calc
    T = T_calc
    Z0_Er = np.append(Z0_Er, Impedance(H,W,T,Er))


# Plot results
plt.scatter(tol, Z0_H, label = 'H')
plt.scatter(tol, Z0_W, label = 'W')
plt.scatter(tol, Z0_T, label = 'T')
plt.scatter(tol, Z0_Er, label = '$\epsilon_r$')
plt.legend(fontsize = '21', loc = 'best')
plt.xlabel('Change from calculated value (%)', fontsize = '21')
plt.ylabel('|Impedance| ($\Omega$)', fontsize = '21')
plt.show()

# Show % change in values
print(percentDiff(Z0_H, Z0_calc))
print(percentDiff(Z0_W, Z0_calc))
print(percentDiff(Z0_T, Z0_calc))
print(percentDiff(Z0_Er, Z0_calc))
