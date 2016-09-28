scp mull-login:simrna/simRNA/flat$argv[1]grem/batch$argv[2]/models/model$argv[3].cas best.cas
~/util/sap121 true.pdb best.cas | grep Wei
rasmol -script ../super.ras
./scoreWrms links.dat best.cas 8 4 50 | tail -2
~/util/sap121 smoothB.pdb smoothA.pdb | grep Wei
rasmol -script ../super.ras

