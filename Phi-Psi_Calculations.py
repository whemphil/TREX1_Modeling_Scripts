# Obvious function

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Dihedral
from MDAnalysis.analysis.dihedrals import Ramachandran

u = mda.Universe('mT1.mol2')
r=u.select_atoms('resid 1-462')
R=Ramachandran(r).run()

data=R.angles

Phi=data[:,:,0]
Psi=data[:,:,1]

import csv

with open("phi.csv","w+") as my_csv:
    csvWriter = csv.writer(my_csv,delimiter=';')
    csvWriter.writerows(Phi)
with open("psi.csv","w+") as my_csv2:
    csvWriter = csv.writer(my_csv2,delimiter=';')
    csvWriter.writerows(Psi)
