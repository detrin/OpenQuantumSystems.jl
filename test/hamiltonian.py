
import numpy as np

import quantarhei as qr
from quantarhei import Molecule
from quantarhei import Mode
from quantarhei import Aggregate

Nmol = 2
Nvib = 3
omega = 300.
HR = 0.01

mol_l = []
with qr.energy_units("1/cm"):
    mol = Molecule([12500., 12750.])

    mod = Mode(omega)
    mol.add_Mode(mod)
    mod.set_nmax(0, Nvib)
    mod.set_nmax(1, Nvib)
    mod.set_HR(1, HR)
    mol_l.append(mol)

    mol = Molecule([12500., 12800.])

    mod = Mode(omega)
    mol.add_Mode(mod)
    mod.set_nmax(0, Nvib)
    mod.set_nmax(1, Nvib)
    mod.set_HR(1, HR)
    mol_l.append(mol)

    agg = Aggregate(molecules=mol_l)
    agg.init_coupling_matrix()
    for i in range(Nmol - 1):
        agg.set_resonance_coupling(i, i+1, 100)

    agg.build(mult=1)
    Ham_py = agg.get_Hamiltonian().data

Ham_jl = np.genfromtxt("hamiltonian.csv", delimiter=',')

# quantarhei does not work with h*omega*(n + 1/2), but with h*omega*n
for i in range(len(Ham_jl)):
    Ham_jl[i, i] -= omega

Ham_diff = np.abs(Ham_py - Ham_jl)

diff = np.sum(Ham_diff)

'''
import matplotlib.pyplot as plt

plt.imshow(Ham_diff)
plt.colorbar()
plt.show()
'''
print(diff)
