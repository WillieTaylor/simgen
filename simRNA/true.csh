make
if ( -e models ) then
	rm models/*
else
	mkdir models
endif
if ( -e rms.dat ) mv rms.dat rms.save.dat
@ n = 0
@ best = 9999
@ high = -999
while ( 1 )
	@ n++
	./sim true.run > sim.log
	grep ' E ' dump.pdb | sort -n -k6 > sort.pdb
	tcsh ~/util/renumpdb.csh sort.pdb temp.pdb
	./scoreWrms links.dat temp.pdb 8 4 50 > score.log
	set bum = `head -1 score.log | awk '{print $3}'`
	set rms = `tail -1 score.log | awk '{print $6}'`
	set sco = `tail -2 score.log | head -1 | awk '{print $12}'`
	echo $n $rms $sco $bum >> rms.dat
	echo Model $n  rms = $rms   score = $sco  bumps = $bum
	echo REMARK   rms = $rms   score = $sco  bumps = $bum > models/model$n.cas
	cat temp.pdb >>  models/model$n.cas
	# keep best rms model
	@ rmsi = `echo $rms | awk '{printf("%d\n", $1*100)}'`
	if ( $rmsi < $best ) then
		echo new lowest rmsd = $rms  score = $sco
		cp temp.pdb best.pdb
		 ~/util/sap121 smoothA.pdb smoothB.pdb | grep 'over all'
		rasmol -script super.ras super.pdb &
		@ best = $rmsi
	endif
	# keep highest scoring
	@ scoi = `echo $sco | awk '{printf("%d\n", $1)}'`
	if ( $scoi > $high ) then
		echo new highest score = $sco  rms = $rms
		cp temp.pdb high.pdb
		tcsh main/view.csh links.dat temp.pdb 50 5
		@ high = $scoi
	endif
	rm temp.pdb
	rm sort.pdb
end

#gnuplot> plot [][3000:] 'rms.dat' u 2:($3*$3-($4*10)) not, 3333 not, 'knot.dat' u 2:($3*$3-($4*10)) not
#gnuplot> plot [7:][50:] 'rms.dat' u 2:($3*$3/($4+10)) not,  110 not, 'knot.dat' u 2:($3*$3/($4+10)) not

#wtaylor@wt:~/simgen/simRNA$ cat rms.dat | awk '{print $0, $3*$3-($4*10)}' | sort -nr -k5 | head -20
#wtaylor@wt:~/simgen/simRNA$ cat rms.dat | awk '{print $0, $3*$3/($4+10)}' | sort -nr -k5 | head -20

