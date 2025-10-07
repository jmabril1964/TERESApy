"""
==> @ Crreated by J.M. Abril-Hernández, 2025.

This code is based on TERESA_clouds.py, but defines a different threshold for X_gl using the concept of 
"interesting parameters" (see G. Cowan, 1998: *Statistical Data Analysis*). 

The interesting parameters are wm0 and sgw0, as they strictly determine the sequence of SARs and, 
consequently, the chronology. 

The method keeps the values of the other two parameters fixed at their absolute minimum: minA and minw. 

The code reads from 'Mapa4D.txt' only those solvers that contain minA and minw (their total number is NR²), 
and stores in a single output file ('Cloud_ages.txt') those that fall within the new threshold value of X_gl, 
corresponding to the 68.3% confidence level (equivalent to 1σ) under these settings.

==> If you keep the default names for the input and output files, you can run this code without any further action. 

"""
import sys


with open('Absolute_min.txt', 'r') as file:
    for line in file:    
        if line.strip():
            parts = line.strip().split()
            minA = float(parts[0])  
            minw= float(parts[1])  
            minsa= float(parts[2])
            minsw= float(parts[3])
            minChi= float(parts[4])
            N = int(parts[5])
            break
# Read only the first non-empty line, which contains all required parameters.

minChiT = (minChi**2 + 2.30/(N-4))**(0.5)  # # Threshold for 68.3% confidence level (1σ) for two interesting parameters (Δχ² = 2.30)
count1 = count2 = 0

with open('Cloud_ages.txt', 'w', encoding='utf-8') as file1:
    with open('Mapa4D.txt', 'r') as file0:
        for line in file0:
            if line.strip():
                parts = line.strip().split()
                if len(parts) != 5:
                    print("Fallo != 5\n")
                    break 
                Aini = float(parts[0])
                saini = float(parts[2])
                vChi = float(parts[4])

                # Only when the double-conditions is true
                if Aini == minA and saini == minsa:
                    if vChi <= minChiT:
                        file1.write(line)
                        count1 += 1
                    else:
                        count2 += 1
                else:
                    pass  # Nothing to do

print(f"Counts (1) {count1} , (2) {count2} \n")
suma = count1 + count2 
print(f"Total solvers processed {suma}\n")
sys.exit()