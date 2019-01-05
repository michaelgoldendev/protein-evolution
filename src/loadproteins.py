from Bio.PDB import *
from Bio import AlignIO
import numpy as np
import os
import json
import subprocess

def fasttree_protein(fastafile):
    infile = os.path.abspath(fastafile)
    executable = os.path.abspath("../binaries/FastTree.exe")
    cmd = "%s -lg %s > %s.nwk" % (executable, infile, infile)
    subprocess.run(cmd, shell=True)
    fin = open("%s.nwk" % infile, "r")
    newick = fin.readlines()[0].strip()
    return newick

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

def get_bond_angles_list(polypeptide):
        bond_angles = [(-1000.0, -1000.0, calc_angle(polypeptide[0]['N'].coord, polypeptide[0]['CA'].coord, polypeptide[0]['C'].coord).item())]

        for aa in range(1,len(polypeptide)):
        	b1 = calc_angle(polypeptide[aa-1]['CA'].coord, polypeptide[aa-1]['C'].coord, polypeptide[aa]['N'].coord)
        	b2 = calc_angle(polypeptide[aa-1]['C'].coord, polypeptide[aa]['N'].coord, polypeptide[aa]['CA'].coord)
        	b3 = calc_angle(polypeptide[aa]['N'].coord, polypeptide[aa]['CA'].coord, polypeptide[aa]['C'].coord)
        	bond_angles.append((b1.item(),b2.item(),b3.item()))
        return bond_angles

def get_bond_lengths(polypeptide):
        bond_lengths = [(-1000.0, norm(polypeptide[0]['CA'].coord - polypeptide[0]['N'].coord).item(), norm(polypeptide[0]['C'].coord - polypeptide[0]['CA'].coord).item())]

        for aa in range(1,len(polypeptide)):
        	b1 = norm(polypeptide[aa]['N'].coord - polypeptide[aa-1]['C'].coord)
        	b2 = norm(polypeptide[aa]['CA'].coord - polypeptide[aa]['N'].coord)	
        	b3 = norm(polypeptide[aa]['C'].coord - polypeptide[aa]['CA'].coord)
        	bond_lengths.append((b1.item(),b2.item(),b3.item()))
        return bond_lengths        

def insertgaps(alignmentseq, vec, gapvalue=-1000.0):
    index = 0
    ret = []
    for c in alignmentseq:
        if c == '-':
            ret.append(gapvalue)
        else:
            ret.append(vec[index])
            index += 1
    return ret

def countnongaps(alignedsequences):
    nongapcounts = np.zeros(len(alignedsequences[0]))
    for seq in alignedsequences:
        for col,c in enumerate(seq):
            if c != '-':
                nongapcounts[col] += 1.0
    return nongapcounts

def removeonlygaps(seq, alignedsequences):
    nongapcounts = countnongaps(alignedsequences)
    ret = ""
    for (count,c) in zip(nongapcounts,seq):
        if count > 0.0:
            ret += c
    return ret

class AlignedProtein:
    def __init__(self):
        self.name = ""
        self.sequence = ""
        self.aligned_sequence = ""
        self.phi_psi = []
        self.aligned_phi_psi = []
        self.omega = []
        self.aligned_omega = []
        self.bond_angles = []
        self.aligned_bond_angles = []
        self.bond_lengths = []
        self.aligned_bond_lengths = []

    def toJSON(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

pdbl = PDBList()

alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"


homstradpath = "../data/homstrad_with_PDB_2018_Nov_1/"
for homstradfolder in os.listdir(homstradpath):
    homstradfamilypath = os.path.join(homstradpath,homstradfolder)
    print(homstradfamilypath)
    if os.path.isdir(homstradfamilypath):
        for alignmentfile in os.listdir(homstradfamilypath):
            if alignmentfile.endswith(".ali"):
                alignment = AlignIO.read(open(os.path.join(homstradfamilypath, alignmentfile)), "pir")
                family = []
                for recordindex, record in enumerate(alignment):                    
                    try:
                        pdbfile = os.path.join(homstradfamilypath, alignmentfile)[:-4] + "-sup.pdb"
                        #print(homstradpdbfile)
                        #pdbfile = pdbl.retrieve_pdb_file(record.id[:4], file_format="pdb", pdir="pdbs/")
                        if os.path.exists(pdbfile):
                            chainid = str(alphabet[recordindex])

                            p = PDBParser()
                            structure = p.get_structure('X', pdbfile)
                            peptideindex = 0
                            for model in structure:
                                for chain in model:
                                    if chainid == None or chain.id.upper() == chainid:
                                        #print(record.description)
                                        ppb=PPBuilder()
                                        alignedprotein = AlignedProtein()
                                        for peptide in ppb.build_peptides(chain):
                                            alignedprotein.sequence += str(peptide.get_sequence())
                                            alignedprotein.phi_psi.extend(peptide.get_phi_psi_list())
                                            alignedprotein.omega.extend(get_omega_list(peptide))
                                            alignedprotein.bond_angles.extend(get_bond_angles_list(peptide))
                                            alignedprotein.bond_lengths.extend(get_bond_lengths(peptide))
                                        if str(record.seq).replace("-","") == alignedprotein.sequence:
                                            alignedprotein.name = str(record.id)
                                            alignedprotein.aligned_sequence = str(record.seq)
                                            for pos in range(len(alignedprotein.phi_psi)):
                                                if alignedprotein.phi_psi[pos][0] == None and alignedprotein.phi_psi[pos][1] == None:
                                                    alignedprotein.phi_psi[pos] = (-1000.0, 1000.0)
                                                elif alignedprotein.phi_psi[pos][0] == None:
                                                    alignedprotein.phi_psi[pos] = (-1000.0, alignedprotein.phi_psi[pos][1])
                                                elif alignedprotein.phi_psi[pos][1] == None:
                                                    alignedprotein.phi_psi[pos] = (alignedprotein.phi_psi[pos][0], -1000.0)
                                            alignedprotein.aligned_phi_psi = insertgaps(record.seq, alignedprotein.phi_psi, gapvalue=(-1000.0, -1000.0))
                                            alignedprotein.aligned_omega = insertgaps(record.seq, alignedprotein.omega, gapvalue=-1000.0)
                                            family.append(alignedprotein)
                                            
                                        else:
                                            print("BAD")
                    except Exception as e:
                        print("err",e)
                alignedseqs = [protein.aligned_sequence for protein in family]
                if len(family) > 0:
                    print(len(family))
                    print(alignedseqs)
                    print(countnongaps(alignedseqs))
                    fastafile = "../data/families/%s_%d.fas" % (alignmentfile, len(family))
                    fastaout = open(fastafile, "w")
                    for alignedprotein in family:
                        fastaout.write(">%s\n" % alignedprotein.name)
                        alignedprotein.aligned_sequence = removeonlygaps(alignedprotein.aligned_sequence, alignedseqs)
                        alignedprotein.aligned_phi_psi = insertgaps(alignedprotein.aligned_sequence, alignedprotein.phi_psi, gapvalue=(-1000.0, -1000.0))
                        alignedprotein.aligned_omega = insertgaps(alignedprotein.aligned_sequence, alignedprotein.omega, gapvalue=-1000.0)
                        alignedprotein.aligned_bond_angles = insertgaps(alignedprotein.aligned_sequence, alignedprotein.bond_angles, gapvalue=(-1000.0, -1000.0, -1000.0))
                        alignedprotein.aligned_bond_lengths = insertgaps(alignedprotein.aligned_sequence, alignedprotein.bond_lengths, gapvalue=(-1000.0, -1000.0, -1000.0))
                        fastaout.write(alignedprotein.aligned_sequence+"\n")                        
                    fastaout.close()
                    newick_tree = fasttree_protein(fastafile)
                    fout = open("../data/families/%s_%d.fam" % (alignmentfile, len(family)), "w")
                    fout.write("{\n")
                    fout.write("\"newick_tree\": \"%s\",\n" % newick_tree)                    
                    fout.write("\"proteins\": [")
                    fout.write(",".join([alignedprotein.toJSON() for alignedprotein in family]))
                    fout.write("]\n")

                    fout.write("}")
                    fout.flush()
