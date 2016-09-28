./sim sams.run
grep ' E ' dump.pdb | sort -n -k6 > sort.pdb
~/util/sap121 2gisP.pdb sort.pdb | grep Weighted
#rasmol -script super.ras super.pdb
