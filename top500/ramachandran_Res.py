##
# Compute density plot for each aminoacid usin top500.tsv proteins
#
#      CLeandro 2014
#
# Please see accompanying webpage:
# 
# http://www.warwick.ac.uk/go/peter_cock/python/ramachandran/calculate/
#
# Parte of this code was relied on Thomas Hamelryck's Bio.PDB module in BioPython:
#
# http://www.biopython.org
#
# It assumes the input files pdb is in the current directory,
# and generates an output file tsv in the current directory.
import matplotlib.pyplot as mpl
import math
import numpy as np
import scipy
import matplotlib.pyplot as plt

input_file = open("top500.tsv","r")
lines = input_file.readlines()
input_file.close()

x=[]
y=[]
#RES={'PHE':0, 'ASP':0, 'THR':0, 'ARG':0, 'TRP':0, 'VAL':0, 'CYS':0, 'SER':0, 'ALA':0, 'GLY':0, 'MET':0, 'TYR':0, 'ASN':0, 'PRO':0, 'LYS':0, 'HIS':0, 'GLN':0, 'ILE':0, 'LEU':0, 'GLU':0}
RES={'PHE':0, 'THR':0, 'ARG':0, 'TRP':0,  'SER':0, 'ALA':0, 'GLY':0, 'MET':0, 'TYR':0, 'ASN':0, 'PRO':0, 'LYS':0, 'HIS':0, 'ILE':0, 'LEU':0, 'GLU':0}
print "About to save angles to file..."
for res in RES.keys():
    output_file = open("biopython_%s.tsv"%(res),"w")
    for line in lines:
        (pdb_code, chains) = line.split("\t")
        print "About to load Bio.PDB and the PDB file", pdb_code
        import Bio.PDB
        structure = Bio.PDB.PDBParser().get_structure(pdb_code, "%s.pdb" % pdb_code)
        print "Done"
        for model in structure :
            for chain in model :
                print "Chain %s"%(str(chain.id))
                polypeptides = Bio.PDB.CaPPBuilder().build_peptides(chain)
                for poly_index, poly in enumerate(polypeptides) :
                    phi_psi = poly.get_phi_psi_list()
                    resslist= enumerate(poly)
                    for res_index, residue in  resslist:
                            phi, psi = phi_psi[res_index]
                            if residue.get_resname()== res and phi and psi :
                                RES[res]=RES[res]+1
                                x.append(degrees(phi))
                                y.append(degrees(psi))
                                #Don't write output when missing an angle
                                output_file.write("%f,%f\n"% (phi, psi))
    output_file.close()
    #xmin = min(x)
    #xmax = max(x)
    #ymin = min(y)
    #ymax = max(y)

    #plt.subplots_adjust(hspace=0.5)
    #plt.subplot(121)
    #plt.hexbin(x,y, cmap=plt.cm.YlOrRd_r)
    #plt.axis([xmin, xmax, ymin, ymax])
    #plt.title("Ramachandran")
    #cb = plt.colorbar()
    #cb.set_label('counts')

    #plt.subplot(122)
    #plt.hexbin(x,y,bins='log', cmap=plt.cm.YlOrRd_r)
    #plt.axis([xmin, xmax, ymin, ymax])
    #plt.title("With a log color scale")
    #cb = plt.colorbar()
    #cb.set_label('log10(N)')
    #plt.savefig('../ramachandran_500%s.eps'%(res))
    #plt.show()

print "Done"
