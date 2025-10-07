"""
==> @ Crreated by J.M. Abril-Hernández, 2025.

This code reads a file containing experimental data on 210Pb_exc in sediment core layers.
You must specify the name and path of the file. It also reads a file from the 'aleat/' folder
with the same number of entries.

Required input parameters:
- An initial estimate of the mean value and relative deviation for activity concentrations and SARs.
- Definition of scanning intervals around these values to generate solvers.
- The number NR of divisions per interval. By default, the same NR is used for all four parameters,
  resulting in NR**4 solvers.
- Optional: time-mark information to define an objective function.

The code defines a sorting function to find the best arrangement of (A0i, wi) pairs for each solver.
This means the sorting that minimizes the quadratic distance to the experimental profile.

The adjusted chi-value (X_gl) is computed for the best sorting of each solver.

Outputs:
- `Mapa4D.txt`: contains the four input values and corresponding X_gl for each solver.
- `Absolute_min.txt`: contains the absolute minimum found across all solvers.

These outputs can be used by other scripts for plotting and further analysis.

The parts of the code that need your action are bounded by lines ++++++++++++++

"""
import numpy as np
# from scipy.stats import norm
# import matplotlib.pyplot as plt
import random
import sys

# ++++++++++++++++++++++++  INFORMATION FOR USING A TIME-MARK ++++++++++++++++++++++
         # Please, ignore this part if you are not using a time-mark.
kr = 9  # index for the sediment slice that contains the time mark (note that in Python index stars at zero) 
         # It must be an intger >= 0 and < N (the number of slices in the core)
Tmr = 43 # ' The known age of the slice with the time mark (yr)
sgt = 1.5 # 1-sigma error of Tmr. By default use 1.0
peso = 0* 1 / 2   # Weight to be used in the Objective function. "0" means "no time mark"
                  # 1/2 is a recommended value that you can tune.

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

m_i = []  # mass depth at the bottom of the layer i
Aexp_i = []  # 210Pb_exc (Bq/kg) in layer i
sgAexp_i = [] # Error in 210Pb_exc (Bq/kg) in layer i

# ++++++++++++++++++++++++ UPDATE HERE ONLY THE NAME OF THE TEXT FILE CONTAINING EXPERIMETAL DATA ++++++
with open('Core_C1.txt', 'r') as file:
    for line in file:
        if line.strip():
            parts = line.strip().split()
            m_i.append(float(parts[0]))  
            Aexp_i.append(float(parts[1]))  
            sgAexp_i.append(float(parts[2]))  


# Reading the randomly sorted canonical representative sample of size N of a Normal typified distribution
# This does not require action if the library has been properly built (see file Randon_generator).
N = len(m_i)
print(f"Number of sediment slices : {N} \n")

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

ldPb= 0.03118   # 1/yr. decay constant for 210Pb [http://www.lnhb.fr/home/nuclear-data/nuclear-data-table/]


# This is the main part of the code, where the sorting function is defined. 
# For this code only the chi-value is necessary in return.


def sorting(Am0ij, wm0ij,sgA0ij,sgw0ij):
    Time_sol = []  # To store solution for the chronology
    SolA = []  # To store the solution for the sorted initial activity concentrations
    Solw = []   # As above, but for SARs
    Sol_Th = [] # To store for the theoretical profile (the solution including radioactive decay)

    slices = N
    chif = 0.
    tant = 0.
    Am_ini = Am0ij *(1.0 + sgA0ij * z_1)
    wm_ini = wm0ij *(1.0 + sgw0ij * z_2)
    wm_ini = np.clip(wm_ini, 0.2*wm0, None)  # This is to protect against non physically meaningful values
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
    
    # This closes the k-loop
    #   DEFINITION OF THE OBJETIVE FUNCTION
    chif = chif + peso * N * ((Tmr - Time_sol[kr])/sgt)**2 # see more in publications

    chi_gl=(chif/(N-4))**(0.5)
    return chi_gl

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#           This requires action: INPUT PARAMETERS
# Initial activity concentration Am0 can be estimated from the value of the first slices
# and/or you can use CFCS model (exponential fit) to estimate this value and mean SAR (wm0)
# For relative deviations (sgA0 and sgw0) you can start with recommended values in the range
# 0.1 - 0.3, depending on how "noisy" the experimental profile is. 
# For complex cases requiring values higher than 0.3 -0.4 better to use log_normal distributions
# (read more in publications). 
# If the LN[A(m)] plot shows discontinuities, think about using Multimodal-TERESA (see publications)
# Recommendation: if not sure, star with sgA0 = 0.25, sgw0= 0.25

Am0 = 563.0  # Bq/kg
wm0 = 0.0385  # g/(cm^2 yr).
sgA0 = 0.19
sgw0 = 0.21   

# A factora defines the interval from Am0*(1-factora) to Am0*(1+factora)
# For example, for Am0 =100 and factora = 0.5, the interval is [50, 150]
# The same for the other factors
# It is recommened running this code at least twice. In the first run you can
# define wide intervals, particularly for those magnitudes with more uncertainty
# The model output identifies the region for the absolute minimum.
# In a second run you can use the above values at the minimum as input parameters
# and using factors now defining narrower intervals with a higher resolution

factora = 0.2 
factorw = 0.25 # For very low wm0 you can try alternative definitions of the interval
factorsa = 0.25
factorsw = 0.30 

# The number NR of divisions of each interval NR
# Note that NR = 10 leads to 10**4 solvers, NR = 50 leads 6.25*10^6
# The CPU time increases with NR. You can find a reasonable compromise between resolution
# and CPU time by using NR < 50 and repeated runs

NR = 30
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

dam = factora * Am0 *2/NR
dwm = factorw * wm0 *2/NR
das = factorsa * sgA0 *2/NR
dws = factorsw * sgw0 *2/NR
Chimin= 2.0E20

with open('Mapa4D.txt', 'w', encoding= 'utf-8') as file:
    for ia in range(NR):
        Am0ij = Am0 * (1-factora) + ia * dam
        print(f"Counter ia from 0 to NR = {ia}")
        for iw in range(NR):
            wm0ij = wm0 * (1-factorw) + iw * dwm
            wm0ij = np.clip(wm0ij, 0.2*wm0, None) # Protection
            for isa in range (NR):
                sgA0ij = sgA0 * (1-factorsa) + isa * das
                for isw in range(NR):
                    sgw0ij = sgw0 * (1-factorsw) + isw * dws
                    valor= sorting(Am0ij, wm0ij,sgA0ij,sgw0ij)
                    if valor < Chimin:
                        Chimin = valor
                        minA = Am0ij
                        minw = wm0ij
                        minsa = sgA0ij
                        minsw = sgw0ij               
                    file.write(f"{Am0ij:.2f}\t{wm0ij:.4f}\t{sgA0ij:.4f}\t{sgw0ij:.4f}\t{valor:.3f}\n")
                file.write(f"\n")
    
print(f"Minimum  A= {minA:.2f}, w = {minw:.4f}, sa={minsa:.4f}, sw ={minsw:.4f}, chi_gl ={Chimin:.3f}, slices= {N}\n")
with open('Absolute_min.txt', 'w', encoding= 'utf-8') as file4:
    file4.write(f"{minA:.2f}\t{minw:.4f}\t{minsa:.4f}\t{minsw:.4f}\t{Chimin:.3f}\t{N}\n")
    file4.write(f"{dam/2:.2f}\t{dwm/2:.5f}\t{das/2:.5f}\t{dws/2:.5f}\n")
    file4.write(f"{Am0:.2f}\t{wm0:.5f}\t{sgA0:.5f}\t{sgw0:.5f}\n")
    file4.write(f"{NR}\t{factora:.3f}\t{factorw:.3f}\t{factorsa:.3f}\t{factorsw:.3f}\n")
    file4.write(f"{kr}\t{Tmr:.1f}\t{sgt:.2f}\t{peso:.3f}\n")
    file4.write(f"# Line 1: Aomin, wmin, sAmin, swmin, X_gl, N slices \n")
    file4.write(f"# Line 2: half resolution intervals for Ao, w, sA, sw \n") 
    file4.write(f"# Line 3: Input parameters Am0, wm0, sgA0, sgw0 \n")
    file4.write(f"# Line 4: NR, factora, factorw, factorsa, factorsw \n")
    file4.write(f"# Line 5: Information for time mark: kr, Tmr, sgt, peso \n")
# print(f"resolution dam = {dam/2:.2f}, dwm= {dwm/2:.5f}, das= {das/2:.5f}, dws ={dws/2:.5f}\n")
print("Bye!!")
sys.exit()

 