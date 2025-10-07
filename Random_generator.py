"""
==> @ Crreated by J.M. Abril-Hernández, 2025.

This code generates the so-called 'canonical representative sample' of size N
for a Normal typified distribution (mean = 0, sigma = 1) denoted as z_0. 
This set of number is copied twice in lists z_1 and z_2 and randomly rearranged.
The result is stored in a 2-column text file named "aleat_N.txt"
The process is repeated for N in range (5,100) [You can increase the range if necesary]
The output files will be created in the same folder where the code actually is placed.
It is suggested to create then a sub-folder named "aleat" and sotre all the aleat... files into it.
This will be your library where the other codes will read the necessary information.
Note that your library is "unique" since, although all the users will share the same z_0,
the random rearrangement can be different from a user to another.
Using a library ensure the repeatability of computations, and separating the effect of
random sorting from the variability in the distributions of physical magnitudes Ao and SAR.

"""
import numpy as np
from scipy.stats import norm
import random

def generar_normal_canonica(N):
    probabilidades = [(i + 0.5) / N for i in range(N)]
    valores = [norm.ppf(p) for p in probabilidades]
    
    return np.array(valores)

# Generating your library of random samples of Normal typified distributions

for N in range(5,100):
    file = f"aleat_{N}.txt"
    z_0 = generar_normal_canonica(N)
    z_1 = z_0.copy()
    random.shuffle(z_1)
    z_2 = z_0.copy()
    random.shuffle(z_2)
    with open(file, 'w', encoding='utf-8') as archivo:
        for x,y in zip(z_1, z_2):
            archivo.write(f"{x:.4f}\t{y:.4f}\n")
print("Bye!!")       