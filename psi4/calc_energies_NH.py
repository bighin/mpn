import time
import numpy as np
import psi4

def print_array(ar):
        print(' '.join(map(str, ar)))

# Set memory
psi4.set_memory('2 GB', quiet=True)
psi4.core.set_output_file('output.dat', False)

numpy_memory = 2

mol = psi4.geometry("""
N
H 1 1.04
symmetry c1
""")

psi4.set_options({'basis': '6-31G',
                  'scf_type': 'pk',
                  'mp2_type': 'conv',
                  'freeze_core': 'false',
                  'e_convergence': 1e-8,
                  'd_convergence': 1e-8})

# First compute RHF energy using Psi4
scf_e, wfn = psi4.energy('SCF', return_wfn=True)

# Grab data from
C = wfn.Ca()
ndocc = wfn.doccpi()[0]
nmo = wfn.nmo()
SCF_E = wfn.energy()
eps = np.asarray(wfn.epsilon_a())

# Compute size of SO-ERI tensor in GB
nmo = wfn.nmo()
ERI_Size = (nmo  ** 4) * 128e-9
print('# Size of the SO ERI tensor will be %4.2f GB.' % ERI_Size)
memory_footprint = ERI_Size * 5.2
if memory_footprint > numpy_memory:
    clean()
    raise Exception("Estimated memory utilization (%4.2f GB) exceeds numpy_memory \
                     limit of %4.2f GB." % (memory_footprint, numpy_memory))

# Integral generation from Psi4's MintsHelper
t = time.time()
mints = psi4.core.MintsHelper(wfn.basisset())
print('# Total time taken for ERI integrals: %.3f seconds.' % (time.time() - t))

# Make spin-orbital MO antisymmetrized integrals
I_mo = np.asarray(mints.mo_spin_eri(C, C))

# Update nocc and nvirt
nso = nmo * 2
nocc = ndocc * 2
nvir = nso - nocc

# Build epsilon tensor (from EP2_SO.py)
eps = np.repeat(eps, 2)
eocc = eps[:nocc]
evir = eps[nocc:]

print("nso", end = " ")
print(int(nso))

print("nocc", end = " ")
print(int(nocc))

print("nvirt", end = " ")
print(int(nvir))

print("eocc", end = " ")
print_array(eps[:nocc])

print("evirt", end = " ")
print_array(eps[nocc:])

print ("hfe %4.11f" % SCF_E)
print("enuc",mol.nuclear_repulsion_energy())

H=np.asarray(mints.ao_kinetic()) + np.asarray(mints.ao_potential())
C_nparray = np.asarray(C)
mo_H = np.dot(np.dot(C_nparray.T, H), C_nparray)
mo_Hdiag = mo_H.diagonal()
mo_Hdiag = np.repeat(mo_Hdiag,2)

print("hdiag", end = " ")
print_array(mo_Hdiag[:nocc])

hfenergy=psi4.energy('hf')
mp2energy=psi4.energy('mp2')
mp3energy=psi4.energy('mp3')
mp4energy=psi4.energy('mp4')
mp5energy=psi4.energy('mp5')
mp6energy=psi4.energy('mp6')

order0=hfenergy
order2=mp2energy-hfenergy
order3=mp3energy-mp2energy
order4=mp4energy-mp3energy
order5=mp5energy-mp4energy
order6=mp6energy-mp5energy

print("\nHartree-Fock energy: ",order0)
print("MP2 contribution: ",order2)
print("MP3 contribution: ",order3)
print("MP4 contribution: ",order4)
print("MP5 contribution: ",order5)
print("MP6 contribution: ",order6)

print("\nHF/MP2: ",order0/order2)
print("MP2/MP3: ",order2/order3)
print("MP3/MP4: ",order3/order4)
print("MP4/MP5: ",order4/order5)
print("MP5/MP6: ",order5/order6)

print("")
print("MP2/MP3: ",order2/order3)
print("MP2/MP4: ",order2/order4)
print("MP2/MP5: ",order2/order5)
print("MP2/MP6: ",order2/order6)

