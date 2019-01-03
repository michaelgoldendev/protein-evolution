import numpy as np
from Bio.PDB import *

import json

with open("../data/families/5_3_exonuclease.ali_2.fam") as f:
    data = json.load(f)

omega = data[0]["omega"]
phi_psi = data[0]["phi_psi"]
sequence = data[0]["sequence"]

def norm(v):
	return np.sqrt((np.sum(v**2)))
  
def angle(v1, v2):
    """Return angle between two vectors."""
    n1 = norm(v1)
    n2 = norm(v2)
    c = np.dot(v1, v2) / (n1 * n2)
    # Take care of roundoff errors
    c = min(c, 1)
    c = max(-1, c)
    return np.arccos(c)

def calc_angle(v1, v2, v3):
    v1 = v1 - v2
    v3 = v3 - v2
    return angle(v1, v3)

def get_omega_list(polypeptide):
        omegalist = []
        lng = len(polypeptide)
        for i in range(0, lng):
            res = polypeptide[i]
            try:
                n = res['N'].get_vector()
                ca = res['CA'].get_vector()
            except Exception:
                # Some atoms are missing
                omegalist.append(-1000.0)
                continue
            # Phi
            if i > 0:
                rp = polypeptide[i - 1]
                try:
                    cap  = rp['CA'].get_vector()
                    cp = rp['C'].get_vector()
                    omega = calc_dihedral(cap, cp, n, ca)
                    omegalist.append(omega)
                except Exception:
                    omegalist.append(-1000.0)
            else:
                omegalist.append(-1000.0)
        return omegalist

p = PDBParser()
structure = p.get_structure('X', "../data/homstrad_with_PDB_2018_Nov_1/5_3_exonuclease/5_3_exonuclease-sup.pdb")
chain = structure[0]["A"]
ppb=PPBuilder()
peptidein = ppb.build_peptides(chain)[0]

aa = 16
bondangle = calc_angle(peptidein[aa-1]['CA'].coord, peptidein[aa-1]['C'].coord, peptidein[aa]['N'].coord)
print((bondangle)/np.pi*180.0)

bondangle = calc_angle(peptidein[aa-1]['C'].coord, peptidein[aa]['N'].coord, peptidein[aa]['CA'].coord)
print((bondangle)/np.pi*180.0)

bondangle = calc_angle(peptidein[aa]['N'].coord, peptidein[aa]['CA'].coord, peptidein[aa]['C'].coord)
print((bondangle)/np.pi*180.0)


#ATOM      1  N   MET A   1      -4.308  -7.299   1.075  1.00 52.25
#ATOM      2  CA  MET A   1      -3.695  -7.053   2.378  1.00 57.02
#ATOM      3  C   MET A   1      -4.351  -5.868   3.097  1.00 55.93


#N1 = np.array([0.0, 0.0, 0.0]) 
#CA1 = np.array([0.61300015, 0.24599981, 1.303]) # 1.47
#C1 = np.array([-0.04299974, 1.4309998, 2.0219998] ) # 1.53

N1 = np.array([-4.308, -7.299, 1.075]) 
CA1 = np.array([-3.695, -7.053, 2.378]) # 1.47
C1 = np.array([-4.351, -5.868, 3.097] ) # 1.53

builder = StructureBuilder.StructureBuilder()
structure = builder.init_structure("structure")
builder.init_model(0)
builder.init_chain("A")
builder.init_seg("")

res = builder.init_residue(Polypeptide.one_to_three(sequence[0]), " ", 1, " ")
builder.init_atom("N", N1, 1.0, 1.0, " ", " N ")
builder.init_atom("CA", CA1, 1.0, 1.0, " ", " CA ")
builder.init_atom("C", C1, 1.0, 1.0, " ", " C ")

numaminoacids = len(phi_psi)
startaa = 1
for aa in range(startaa, numaminoacids):
	res = builder.init_residue(Polypeptide.one_to_three(sequence[aa]), " ", aa+1, " ")
	#R = 1.32 # bond length
	R = norm(peptidein[aa]['N'].coord - peptidein[aa-1]['C'].coord)
	torsionangle = phi_psi[aa-1][1]
	bondangle =  calc_angle(peptidein[aa-1]['CA'].coord, peptidein[aa-1]['C'].coord, peptidein[aa]['N'].coord)
	print(torsionangle,calc_dihedral(peptidein[aa-1]['N'].get_vector(), peptidein[aa-1]['CA'].get_vector(), peptidein[aa-1]['C'].get_vector(), peptidein[aa]['N'].get_vector()))
	#origbondangle = bondangle
	#for bondangle in np.linspace(-2.0*np.pi,2.0*np.pi,num=360):
	D = np.array([R*np.cos(bondangle), R*np.cos(torsionangle)*np.sin(bondangle), R*np.sin(torsionangle)*np.sin(bondangle)])
	ab = (CA1-N1)
	ac = (C1-N1)
	bc = (C1-CA1)
	bchat = bc/norm(bc)
	n = np.cross(ab, bchat)/norm(np.cross(ab, bchat))
	M = np.column_stack((bchat,np.cross(n,bchat),n))
	N1 = np.matmul(M,D) + C1
	builder.init_atom("N", N1, 1.0, 1.0, " ", " N ")
	#print(origbondangle,bondangle,np.matmul(M,D) + C1)
	#print("CA1,C,N1", (calc_angle(CA1,C1,N1))/np.pi*180.0)
	#print("N1-C1", norm(N1-C1))

	#R = 1.47 # bond length
	R = norm(peptidein[aa]['CA'].coord - peptidein[aa]['N'].coord)	
	torsionangle = omega[aa]
	#bondangle = phi_psi[aa][0]
	print(torsionangle,calc_dihedral(peptidein[aa-1]['CA'].get_vector(), peptidein[aa-1]['C'].get_vector(), peptidein[aa]['N'].get_vector(), peptidein[aa]['CA'].get_vector()))
	bondangle = calc_angle(peptidein[aa-1]['C'].coord, peptidein[aa]['N'].coord, peptidein[aa]['CA'].coord) 
	D = np.array([R*np.cos(bondangle), R*np.cos(torsionangle)*np.sin(bondangle), R*np.sin(torsionangle)*np.sin(bondangle)])
	ab = (C1-CA1)
	ac = (N1-CA1)
	bc = (N1-C1)
	bchat = bc/norm(bc)
	n = np.cross(ab, bchat)/norm(np.cross(ab, bchat))
	#n = np.cross(ac, bchat)/norm(np.cross(ac, bchat))
	M = np.column_stack((bchat,np.cross(n,bchat),n))
	CA1 = np.matmul(M,D)  + N1
	builder.init_atom("CA", CA1, 1.0, 1.0, " ", " CA ")
	print("C,N1,CA1", (calc_angle(C1,N1,CA1))/np.pi*180.0)
	#print("CA1-N1", norm(CA1-N1))

	if phi_psi[aa][1] > -100.0:
		#R = 1.53 # bond length
		R = norm(peptidein[aa]['C'].coord - peptidein[aa]['CA'].coord)
		torsionangle = phi_psi[aa][0]
		print(torsionangle,calc_dihedral(peptidein[aa-1]['C'].get_vector(), peptidein[aa]['N'].get_vector(), peptidein[aa]['CA'].get_vector(), peptidein[aa]['C'].get_vector()))
		#bondangle = phi_psi[aa][1]
		bondangle = calc_angle(peptidein[aa]['N'].coord, peptidein[aa]['CA'].coord, peptidein[aa]['C'].coord)
		D = np.array([R*np.cos(bondangle), R*np.cos(torsionangle)*np.sin(bondangle), R*np.sin(torsionangle)*np.sin(bondangle)])
		ab = (N1-C1)
		ac = (CA1-C1)
		bc = (CA1-N1)
		bchat = bc/norm(bc)
		n = np.cross(ab, bchat)/norm(np.cross(ab, bchat))
		#n = np.cross(ac, bchat)/norm(np.cross(ac, bchat))
		M = np.column_stack((bchat,np.cross(n,bchat),n))
		C1 = np.matmul(M,D)  + CA1
		builder.init_atom("C", C1, 1.0, 1.0, " ", " C ")
		print("N1,CA1,C", (calc_angle(N1,CA1,C1))/np.pi*180.0)
		#print("C1-CA1", norm(C1-CA1))
io = PDBIO()
io.set_structure(builder.get_structure())
io.save("out.pdb")

ppb=PPBuilder()
chain = builder.get_structure()[0]["A"]
print(chain)

for (old,new) in zip(phi_psi[0:10],ppb.build_peptides(chain)[0].get_phi_psi_list()[0:10]):
	#print(old,new,old[0]+np.pi)
	print(old,new)
for (old,new) in zip(omega[0:10], get_omega_list(ppb.build_peptides(chain)[0])[0:10]):
	#print(old,new,old[0]+np.pi)
	print(old,new)
"""
from Bio.PDB import *


builder = StructureBuilder.StructureBuilder()
structure = builder.init_structure("structure")
builder.init_model(0)
builder.init_chain("A")
builder.init_seg("")
res = builder.init_residue("MET", " ", 1, " ")
builder.init_atom("N", np.array([0.0,0.0,0.0], "f"), 1.0, 1.0, " ", " N ")
builder.init_atom("CA", np.array([0.0,0.0,0.0], "f"), 1.0, 1.0, " ", " CA ")
builder.init_atom("C", np.array([0.0,0.0,0.0], "f"), 1.0, 1.0, " ", " C ")
builder.init_atom("O", np.array([0.0,0.0,0.0], "f"), 1.0, 1.0, " ", " O ")
#res = builder.init_residue("ASN", "H", 1, "1")
#builder.init_atom("N", [1.0,1.0,1.0], 0.0, 0.5, "2", " N ", "N")
print(builder.chain.get_list())


io = PDBIO()
io.set_structure(builder.get_structure())
io.save("out.pdb")
"""