"""
==> @ Crreated by J.M. Abril-Hernández, 2025.

This code is based on TERESA_plots.py (see previous reference), but it computes the theoretical profiles only 
for the set of solvers stored in the file 'Cloud_ages.txt'.

It stores the results in the output file 'Plot_ages.txt', and also generates the same 'Solution.txt' file,
useful if you are only concerned with chronology, allowing you to skip both TERESA_clouds.py and TERESA_plots.py.

The outputs can be used to plot (with gnuplot or other graphic packages) the chronological line derived from the solver 
at the absolute minimum, along with the cloud of solvers corresponding to a 68.3% confidence interval.

==> If you keep the default names for the input and output files, you can run this code without any further action.

"""
import numpy as np
# from scipy.stats import norm
# import matplotlib.pyplot as plt
# import random
import sys


m_i = []  # mass depth at the bottom of the slice i
Aexp_i = []  # 210Pb_exc (Bq/kg) in slice i
sgAexp_i = [] # Error in 210Pb_exc (Bq/kg) in slice i

with open('Core_C1.txt', 'r') as file:
    for line in file:
        if line.strip():
            parts = line.strip().split()
            m_i.append(float(parts[0]))  
            Aexp_i.append(float(parts[1]))
            sgAexp_i.append(float(parts[2]))


# Reading the canonical representative sample of randomly sorted values following a normal typified distribution
N = len(m_i)
file2 = f"./aleat/aleat_{N}.txt"
print(file2)
z_1 = []
z_2 = []
with open(file2, 'r') as file_aleat:
    for line in file_aleat:
        if line.strip():
            parts = line.strip().split()
            z_1.append(float(parts[0]))  
            z_2.append(float(parts[1]))  

z_1 = np.array(z_1)
z_2 = np.array(z_2)

ldPb= 0.03118   # 1/yr. decay constant for 210Pb

# Computing the mass depth scale referred to the mid-point of each sediment slice
mi_m = []
for k in range(N):
    if k == 0:
        dm2=m_i[0]
        mi_m.append(dm2/2)
    else: 
        dm2 = m_i[k]-m_i[k-1]
        mi_m.append(m_i[k-1]+dm2/2)
# print(mi_m)


def profile(Am0ij, wm0ij,sgA0ij,sgw0ij):
    Time_sol = []  # Stores solution for ages
    SolA = []  # Stores solution for initial activity concentrations
    Solw = []   # Stores solution for mean SAR
    Sol_Th = []  # Stores solution for the 210-Pb_exc profile that compares against the empirical one
    Sol_flux = [] # Stores solution for the 210-Pb_exc flux captured by each sediment slice.

    slices = N
    chif = 0.
    tant = 0.
    Am_ini = Am0ij *(1.0 + sgA0ij * z_1)
    wm_ini = wm0ij *(1.0 + sgw0ij * z_2)

    for k in range(N):
        chiant = 1.0E20
        jmin = 0
        if k == 0:
            dm2=m_i[0]
        else: 
            dm2 = m_i[k]-m_i[k-1]
        for s in range(slices):
            Dt = dm2/wm_ini[s]
            Act = Am_ini[s]*np.exp(-ldPb*tant)*(1-np.exp(-ldPb*Dt))/(ldPb*Dt)
            chi = (Act-Aexp_i[k])**2
            if chi < chiant:
                jmin = s
                chiant = chi
        SolA.append(Am_ini[jmin])
        Solw.append(wm_ini[jmin])
        Time_sol.append(tant+dm2/wm_ini[jmin])
        tant = tant+dm2/wm_ini[jmin]
        Am_ini = np.delete(Am_ini, jmin)
        wm_ini = np.delete(wm_ini, jmin)
        slices -= 1
        if k == 0:
            Dt = Time_sol[0]
            tanterior = 0
        else:
            Dt = Time_sol[k]-Time_sol[k-1]
            tanterior = Time_sol[k-1]
        Act = SolA[k]*np.exp(-ldPb*tanterior)*(1-np.exp(-ldPb*Dt))/(ldPb*Dt)
        chif += ((Act-Aexp_i[k])/sgAexp_i[k])**2
        Sol_Th.append(Act)
        Sol_flux.append(SolA[k]*Solw[k]*10)
    
    #   The X_gl value nor the OBJECTIVE function needs to be recalculated here
    return SolA, Solw, Time_sol, Sol_Th, Sol_flux



with open('Plot_ages.txt', 'w', encoding= 'utf-8') as file_out1:
    with open('Cloud_ages.txt', 'r') as file1:
        for line in file1:
            if line.strip():
                parts = line.strip().split()
                Am0ij = float(parts[0])  
                wm0ij = float(parts[1]) 
                sgA0ij = float(parts[2])
                sgw0ij = float(parts[3])
            SolA, Solw, Time_sol, Sol_Th, Sol_flux = profile(Am0ij, wm0ij, sgA0ij, sgw0ij)
            for x1, x2, x3, x4, x5, x6, x7 in zip(m_i, mi_m, SolA, Solw, Time_sol, Sol_Th, Sol_flux):
                file_out1.write(f"{x1:.3f}\t{x2:.3f}\t{x3:.2f}\t{x4:.4f}\t{x5:.2f}\t{x6:.2f}\t{x7:.2f}\n")
            file_out1.write(f"\n")

with open('Solution.txt', 'w', encoding= 'utf-8') as file_out4:
    with open('Absolute_min.txt', 'r') as file4:
        for line in file4:
            if line.strip():
                parts = line.strip().split()
                Am0ij = float(parts[0])  
                wm0ij = float(parts[1]) 
                sgA0ij = float(parts[2])
                sgw0ij = float(parts[3])
                # Additional information in this file is not needed here
            SolA, Solw, Time_sol, Sol_Th, Sol_flux = profile(Am0ij, wm0ij, sgA0ij, sgw0ij)
            for x1, x2, x3, x4, x5, x6, x7 in zip(m_i, mi_m, SolA, Solw, Time_sol, Sol_Th, Sol_flux):
                file_out4.write(f"{x1:.3f}\t{x2:.3f}\t{x3:.2f}\t{x4:.4f}\t{x5:.2f}\t{x6:.2f}\t{x7:.2f}\n")
            break
sys.exit()
