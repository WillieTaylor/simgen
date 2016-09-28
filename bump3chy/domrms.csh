rm domrms.dat
@ lens = `cat final.out | wc -l`
@ len = 129
@ i = 0
@ m = 0
@ bad = -999
grep -v TER final.out > final.pdb
grep -v TER start.out > start.pdb
while ( $m < $lens )
	@ i++
	@ m += $len
	head -$m final.pdb | tail -$len > protA.cas
	@ k = $len - 5
	head -$k protA.cas > tmp.pdb
	@ k = $len - 10
	tail -$k tmp.pdb > protA.cas
	@ j = 0
	@ n = 0
	while ( $n < $lens )
		@ j++
		@ n += $len
		head -$n start.pdb | tail -$len > protB.cas
		@ k = $len - 5
		head -$k protB.cas > tmp.pdb
		@ k = $len - 10
		tail -$k tmp.pdb > protB.cas
		set rms = `~/util/sap121 protA.cas protB.cas | tail -1 | awk '{print $5}'`
		@ got = `echo $rms | awk '{printf("%d",$1*100)}'`
		if ( $got > $bad ) then
			@ bad = $got
			cp super.pdb worst.pdb
		endif
		echo $i $j $rms >> domrms.dat
	end
end 
set ave = `awk '{n++; s+=$3; print s/n}' domrms.dat | tail -1`
#set min = `sort -n -k3 domrms.dat | head -1`
#set max = `sort -n -k3 domrms.dat | tail -1`
set min = `sort -n -k3 domrms.dat | head -1 | awk '{print $3}'`
set max = `sort -n -k3 domrms.dat | tail -1 | awk '{print $3}'`
echo "	RMS: " min = $min  ave = $ave  max = $max
