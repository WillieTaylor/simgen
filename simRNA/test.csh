# 1 = A or B starting form

make
	./sim test$argv[1].run > sim.log
	tail -1 sim.log
	grep ' E ' dump.pdb | sort -n -k6 > sort.pdb
	tcsh ~/util/renumpdb.csh sort.pdb temp.pdb
	./scoreWrms links.dat temp.pdb 8 4 50 > score.log
	set bum = `head -1 score.log | awk '{print $3}'`
	set rms = `tail -1 score.log | awk '{print $6}'`
	set sco = `tail -2 score.log | head -1 | awk '{print $12}'`
	echo $rms $sco $bum 
	rasmol -script chain.ras temp.pdb
